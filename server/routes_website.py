import logging
import os
import traceback
import requests
from threading import Thread
from pprint import pformat

from flask import jsonify, request, Blueprint

from pypka.main import getTitrableSites

from pka2pI import pkas_2_pdb, clean_pdb
from database import DB_SESSION
from models import Job, Protein, Input, Results
from const import PKPDB_PARAMS, DIR_PATH, CONFIG
from routes_direct import pkpdb_query, save_pdb
from utils import submit_pypka_job, get_subID, run_pKAI

routes_bp = Blueprint("routes", __name__)


@routes_bp.route("/getTitrableSitesNumber", methods=["POST"])
def getNumberOfTitratableSites():
    pdbfile = request.json["PDB"]

    subID = get_subID(request)
    newfilename = save_pdb(pdbfile, subID)

    try:
        out_sites, chains_res = getTitrableSites(newfilename, ser_thr_titration=False)
    except Exception as e:
        tb = traceback.format_exc()
        logging.error(tb)
        response_dict = {"error": str(e)}
        response = jsonify(response_dict)
        response.headers.add("Access-Control-Allow-Origin", "*")
        return response

    os.remove(newfilename)

    nchains = len(chains_res.keys())
    nres = 0
    for chain, sites in chains_res.items():
        nres += len(sites)

    response_dict = {"nchains": nchains, "nsites": nres}

    response = jsonify(response_dict)
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


@routes_bp.route("/getSubID", methods=["POST"])
def getSubID():
    subID = get_subID(request)

    response = jsonify({"subID": subID})
    response.headers.add("Access-Control-Allow-Origin", "*")

    return response


@routes_bp.route("/submitSim", methods=["POST"])
def submitCalculation():
    # Get IP, country, and city information
    ip = request.remote_addr
    location_info = requests.get(f"https://ipinfo.io/{ip}/json").json()
    country = location_info.get("country", "Unknown")
    city = location_info.get("city", "Unknown")

    pdbfile = request.json["pdbfile"]
    pdbid = request.json["pdbcode"]

    input_naming_scheme = request.json["inputFileNamingScheme"]

    pHmin = request.json["pHmin"]
    pHmax = request.json["pHmax"]
    pHstep = request.json["pHstep"]

    epsin = request.json["proteinDielectric"]
    epsout = request.json["solventDielectric"]
    ionic = request.json["ionicStrength"]

    nsites = request.json["nsites"]
    nchains = request.json["nchains"]

    outputpKs = request.json["outputpKValues"]
    outputfile = request.json["outputPDBFile"]
    outputfilenaming = request.json["outputFileNamingScheme"]
    outputfilepH = request.json["outputFilepH"]

    on_PKPDB = request.json["onPKPDB"]
    default_parameters = request.json["defaultParams"]

    outputemail = request.json["email"]

    model = request.json["model"]

    subID = get_subID(request)

    newfilename = save_pdb(pdbfile, subID)

    if outputpKs:
        pH = f"{pHmin},{pHmax}"
    else:
        pH = str(outputfilepH)

    if model in ("pKAI", "pKAI+"):
        tmp_f = f"{DIR_PATH}/pdbs/{subID}_cleaned.pdb"
        clean_pdb(newfilename, tmp_f)
        newfilename = tmp_f

        results = run_pKAI(newfilename, model)

        if outputfile:
            pkas_2_pdb(
                subID,
                newfilename,
                f"{DIR_PATH}/pdbs_out/out_{subID}.pdb",
                outputfilepH,
                outputfilenaming,
                results,
            )
            results["pdb_out"] = True
        response = jsonify({"subID": subID, **results, "params": pformat(PKPDB_PARAMS)})
        response.headers.add("Access-Control-Allow-Origin", "*")
        return response

    elif on_PKPDB and default_parameters and not outputfile:
        response = pkpdb_query(pdbid)

        return response
    else:
        parameters = {
            "structure": newfilename,
            "pH": pH,
            "pHstep": pHstep,
            "epsin": epsin,
            "epssol": epsout,
            "ionicstr": ionic,
            "ffinput": input_naming_scheme,
            "scaleM": 2,
            "convergence": 0.1,
            "pbc_dimensions": 0,
            "ncpus": CONFIG["SLURM_JOB_NCORES"],
            "clean": True,
            "keep_ions": True,
            "ser_thr_titration": False,
            "titration_output": f"{DIR_PATH}/titrations/titration_{subID}.out",
            "output": f"{DIR_PATH}/pkas/pKas_{subID}.out",
        }

        if outputfile:
            parameters["structure_output"] = (
                f"{DIR_PATH}/pdbs_out/out_{subID}.pdb",
                outputfilepH,
                outputfilenaming,
            )

        sub_parameters = (
            pdbid,
            pdbfile,
            outputemail,
            nsites,
            nchains,
            ip,
            country,
            city,
        )

        Thread(
            target=submit_pypka_job,
            args=(
                parameters,
                sub_parameters,
                subID,
            ),
        ).start()

        response = jsonify({"subID": subID})
        response.headers.add("Access-Control-Allow-Origin", "*")

        return response


@routes_bp.route("/getSubmissions", methods=["POST", "GET"])
def get_submission():
    submission_IDS = (
        DB_SESSION.query(Job.job_id, Job.dat_time, Protein.pdb_code, Job.sub_id)
        .join(Input, Input.protein_id == Protein.protein_id)
        .filter(Input.job_id == Job.job_id)
        .order_by(Job.sub_id.desc())
        .limit(25)
        .all()
    )
    sub_ids = [i[:] for i in submission_IDS]

    # queued = []
    # sbrun = subprocess.run(f"s-id pedror | grep slurm_", shell=True, capture_output=True)
    # slurm_queue = sbrun.stdout.decode("utf-8").strip()
    # for line in slurm_queue.splitlines():
    #    subID = line.split()[-1].replace('.py', '').replace('slurm_', '')
    #    queued.append(subID)

    finished = DB_SESSION.query(Results.job_id)

    inprogress_IDS = (
        DB_SESSION.query(Job.job_id, Job.dat_time, Job.sub_id)
        .filter(Job.job_id.not_in(finished))
        .all()
    )
    queued_ids = [i[:2] + ("QUEUED", i[2]) for i in inprogress_IDS][::-1]

    response = jsonify(queued_ids + sub_ids)
    response.headers.add("Access-Control-Allow-Origin", "*")

    return response
