from flask import Flask, jsonify, request
from flask_cors import CORS
import requests
from pprint import pformat
import datetime
import os
import smtplib
from dotenv import dotenv_values
import logging
from threading import Thread
import traceback

import json

from pypka.main import getTitrableSites

from database import db_session, Base, engine
from models import Job, Protein, Input, UsageStats, Results

from slurm import create_slurm_file
from pka2pI import pkas_2_titcurve, titcurve_2_pI, exclude_cys, pkas_2_pdb, clean_pdb
from pka2pI import pkas_2_titcurve, titcurve_2_pI, exclude_cys, pkas_2_pdb, clean_pdb

from pkai.pKAI import pKAI

# from pkpdb import retrieve_from_pkpdb, retrieve_pkpdb_titcurve, retrieve_pkpdb_pis

from flask_limiter import Limiter
from flask_limiter.util import get_remote_address


import pandas as pd

df_pkpdb_pkas = pd.read_csv(
    "static/pkas.csv", header=0, sep=";"
)  # , compression="gzip"
df_pkpdb_pI = pd.read_csv(
    "static/isoelectric.csv", header=0, sep=";" # , compression="gzip"
)
df_pkpdb_titcurves = pd.read_csv(
    "static/titrationcurves.csv", header=0, sep=";" # , compression="gzip"
)

Base.metadata.create_all(bind=engine)

STATUS = "live"

logging.basicConfig(filename="server.log", level=logging.DEBUG)

root = logging.getLogger("werkzeug")
root.handlers = logging.getLogger().handlers
root.setLevel(logging.INFO)

app = Flask(__name__, static_url_path="/static")
app.config["SECRET_KEY"] = str(datetime.datetime.today())
app.config["CORS_HEADERS"] = "Content-Type"

CORS(app, resources={r"/*": {"origins": "*"}})

dir_path = os.path.dirname(os.path.realpath(__file__))
config = dotenv_values(f"{dir_path}/../.env")

limiter = Limiter(
    get_remote_address,
    app=app,
    default_limits=["500 per hour"],
    storage_uri="memory://",
)

limiter = Limiter(
    get_remote_address,
    app=app,
    default_limits=["500 per hour"],
    storage_uri="memory://",
)


@app.teardown_appcontext
def shutdown_session(exception=None):
    db_session.remove()
    # print(engine.pool.status())


def plus_one_pkpdb_downloads():
    usage = db_session.query(UsageStats).first()
    usage.pkpdb_downloads += 1
    db_session.commit()


def plus_one_pkpdb_queries():
    usage = db_session.query(UsageStats).first()
    usage.pkpdb_queries += 1
    db_session.commit()


def plus_one_pypka_subs():
    usage = db_session.query(UsageStats).first()
    usage.pypka_subs += 1
    db_session.commit()


@app.route("/stats")
def get_stats():
    stats = db_session.query(
        UsageStats.pkpdb_queries,
        UsageStats.pkpdb_downloads,
        UsageStats.pypka_subs,
        UsageStats.pKAI_subs,
        UsageStats.pKAI_subs,
    ).first()

    response_dict = {
        "pKPDB Queries": stats[0],
        "pKPDB Downloads": stats[1],
        "PypKa Jobs": stats[2],
        "pKAI Jobs": stats[3],
    }

    response = jsonify(response_dict)
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


@app.route("/isoelectric.csv")
@limiter.limit("10 per hour")
@app.route("/isoelectric.csv")
@limiter.limit("10 per hour")
def export_pis():
    plus_one_pkpdb_downloads()
    return app.send_static_file("isoelectric.csv")
    return app.send_static_file("isoelectric.csv")


@app.route("/pkas.csv")
@limiter.limit("10 per hour")
@app.route("/pkas.csv")
@limiter.limit("10 per hour")
def export_pkas():
    plus_one_pkpdb_downloads()
    return app.send_static_file("pkas.csv")
    return app.send_static_file("pkas.csv")


@app.route("/similarity090.csv")
@limiter.limit("10 per hour")
@app.route("/similarity090.csv")
@limiter.limit("10 per hour")
def export_clusters():
    plus_one_pkpdb_downloads()
    return app.send_static_file("similarity090.csv")


@app.route("/test")
@limiter.limit("2 per hour")
def test():
    return jsonify({"status": "live"})


@app.route("/test")
@limiter.limit("2 per hour")
def test():
    return jsonify({"status": "live"})

@app.route("/")
@limiter.limit("100 per hour")
@limiter.limit("100 per hour")
def hello():
    status = "live"
    if STATUS:
        status = STATUS
    return jsonify(
        {
            "status": status,
            "endpoints": """
        /pkas/<idcode>     GET/POST    Retrieves the results from the pKPDB if available, and runs a pKAI calculation otherwise
                           Example: https://api.pypka.org/pkas/4LZT
        
        /pkpdb/<idcode>    GET/POST    Queries the pKPDB database
                           Example: https://api.pypka.org/pkpdb/4LZT

        /pKAI/<idcode>     GET/POST    Runs a pKAI calculation
                           Example: https://api.pypka.org/pkpdb/4LZT
        """,
        }
    )


@app.route("/query/<idcode>")
@limiter.limit("100 per hour")
def pkpdb_query(idcode):
    if idcode == "CRON_JOB":
        idcode = "4lzt"
    else:
        plus_one_pkpdb_queries()

    idcode = idcode.lower()

    df_protein_results = df_pkpdb_pkas.query(f"idcode == '{idcode}'").sort_values(
        ["chain", "residue_number"]
    )
    if len(df_protein_results) > 0:
        results_list = df_protein_results[
            ["chain", "residue_name", "residue_number", "pk"]
        ].values.tolist()

        tit_curve = json.loads(
            df_pkpdb_titcurves.query(f"idcode == '{idcode}'").values.tolist()[0][2]
        )
        tit_x, tit_y = list(tit_curve.keys()), list(tit_curve.values())

        response_dict = {
            "tit_x": tit_x,
            "tit_y": tit_y,
            "pI": df_pkpdb_pI.query(f"idcode == '{idcode}'").values.tolist()[0][2],
            "pKas": results_list,
            "params": pformat(PKPDB_PARAMS),
        }

    else:
        response_dict = {}

    response = jsonify(response_dict)
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


@app.route("/pkpdb/<idcode>", methods=["GET", "POST"])
@limiter.limit("100 per hour")
def exists_on_pkpdb(idcode):
    results = df_pkpdb_pkas.query(f"idcode == '{idcode.lower()}'")
    if len(results) > 0:
        results = True
    else:
        results = False
    response = jsonify(results)
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


def run_pKAI(pdb, model):
    results = pKAI(pdb, model_name=model, device="cpu", threads=2)
    results = [[i[0], i[2], i[1], i[3]] for i in results]
    results = exclude_cys(pdb, results)
    tit_x, tit_y = pkas_2_titcurve(pdb, results)
    pI = titcurve_2_pI(tit_x, tit_y)
    return {"pKas": results, "tit_x": tit_x, "tit_y": tit_y, "pI": round(pI, 2)}


@app.route("/pKAI/<idcode>", methods=["GET", "POST"])
@limiter.limit("100 per hour")
def run_pKAI_idcode(idcode):
    subID = get_subID(request)

    r = requests.get(f"https://files.rcsb.org/download/{idcode}.pdb")
    pdb_content = r.content.decode("utf-8")
    if "The requested URL was not found on this server." in pdb_content:
        results = f"Error: {idcode} not found"

    else:
        newfilename = save_pdb(pdb_content, subID)

        results = run_pKAI(newfilename, "pKAI")

    response = jsonify(results)
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


@app.route("/pkas/<idcode>", methods=["GET", "POST"])
@limiter.limit("100 per hour")
def get_pkas_from_idcode(idcode):
    idcode = idcode.lower()
    subID = get_subID(request)

    df_protein_results = df_pkpdb_pkas.query(f"idcode == '{idcode}'").sort_values(
        ["chain", "residue_number"]
    )
    if len(df_protein_results) > 0:
        results_list = df_protein_results[
            ["chain", "residue_name", "residue_number", "pk"]
        ].values.tolist()
        response_dict = {
            "idcode": idcode.upper(),
            "method": "PypKa (pKPDB)",
            "pdb": f"https://files.rcsb.org/download/{idcode}.pdb",
            "pI": df_pkpdb_pI.query(f"idcode == '{idcode}'").values.tolist()[0][2],
            "pKas": results_list,
            "tit_curve": {
                key: round(value, 2)
                for key, value in json.loads(
                    df_pkpdb_titcurves.query(f"idcode == '{idcode}'").values.tolist()[
                        0
                    ][2]
                ).items()
            },
        }
        return jsonify(response_dict)

    # PDB IDCODE
    if len(idcode) == 4:
        r = requests.get(f"https://files.rcsb.org/download/{idcode}.pdb")
        pdb_content = r.content.decode("utf-8")
    else:
        idcode_upper = idcode.upper()
        r = requests.get(
            f"https://alphafold.ebi.ac.uk/files/AF-{idcode_upper}-F1-model_v4.pdb"
        )
        pdb_content = r.content.decode("utf-8")

    if "ATOM " not in pdb_content:
        return f"Error: {idcode} not found"

    newfilename = save_pdb(pdb_content, subID)

    try:
        # Run pKAI
        results = run_pKAI(newfilename, "pKAI")
        response_dict = {
            "idcode": idcode.upper(),
            "method": "pKAI 1.2.0",
            "pdb": f"https://alphafold.ebi.ac.uk/files/AF-{idcode_upper}-F1-model_v4.pdb",
            **results,
            "tit_curve": {x: y for x, y in zip(results["tit_x"], results["tit_y"])},
        }
        del response_dict["tit_x"]
        del response_dict["tit_y"]
    except Exception as e:
        response_dict = {"idcode": idcode, "Error": str(e)}
        return jsonify(response_dict)

    response = jsonify(response_dict)
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


# @app.route("/run/<idcode>")
# def run(idcode):
#
#     subID = get_subID(request)
#
#     r = requests.get(f"https://files.rcsb.org/download/{idcode}.pdb")
#     pdb_content = r.content.decode("utf-8")
#     if "The requested URL was not found on this server." in pdb_content:
#         return f"Error: {idcode} not found"
#
#     newfilename = save_pdb(pdb_content, subID)
#
#     parameters = {
#         "structure": newfilename,
#         "epsin": 15,
#         "ionicstr": 0.1,
#         "convergence": 0.1,
#         "pbc_dimensions": 0,
#         "ncpus": 1,
#         "clean": True,
#         "ser_thr_titration": False,
#         "titration_output": f"{dir_path}/titrations/titration_{subID}.out",
#         "output": f"{dir_path}/pkas/pKas_{subID}.out",
#     }
#
#     try:
#         response_dict = run_pypka(parameters, subID)
#     except Exception as e:
#         response_dict = {"PDBID": idcode, "Error": str(e)}
#         return jsonify(response_dict)
#
#     response_dict = {"pKs": response_dict["pKas"]}
#     response_dict["PDBID"] = idcode
#     response = jsonify(response_dict)
#     response.headers.add("Access-Control-Allow-Origin", "*")
#     return response


@app.route("/queue-size")
@limiter.limit("5000 per hour")
def get_queue_size():
    # response = jsonify(len(job_queue))
    response = jsonify(0)
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


def send_email(outputemail):
    gmail_user = config["EMAIL"]
    gmail_password = config["PASS"]

    sent_from = gmail_user
    to = outputemail
    subject = "PypKa"
    body = "Hey, whats up?\n\n- You"

    email_text = """\
    From: %s
    To: %s
    Subject: %s

    %s
    """ % (
        sent_from,
        ", ".join(to),
        subject,
        body,
    )

    try:
        server = smtplib.SMTP_SSL("smtp.gmail.com", 465)
        server.ehlo()
        server.login(gmail_user, gmail_password)
        server.sendmail(sent_from, to, email_text)
        server.close()

        logging.info("Email sent!")
    except:
        logging.info("Something went wrong...")


def save_pdb(pdbfile, subID):
    def new_name(subID):
        return f"{dir_path}/pdbs/{subID}.pdb"

    newfilename = new_name(subID)
    while os.path.isfile(newfilename):
        subID += "_"
        newfilename = new_name(subID)

    with open(newfilename, "w") as f_new:
        f_new.write(pdbfile)
    return newfilename


def get_subID(request):
    subID = "{}{}".format(datetime.datetime.today(), hash(str(request.headers)))
    return subID.replace("-", "").replace(":", "").replace(" ", "").replace(".", "")


@app.route("/getTitrableSitesNumber", methods=["POST"])
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


@app.route("/getSubID", methods=["POST"])
def getSubID():
    subID = get_subID(request)

    response = jsonify({"subID": subID})
    response.headers.add("Access-Control-Allow-Origin", "*")

    return response


@app.route("/submitSim", methods=["POST"])
def submitCalculation():
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
        tmp_f = f"{dir_path}/pdbs/{subID}_cleaned.pdb"
        clean_pdb(newfilename, tmp_f)
        newfilename = tmp_f

        results = run_pKAI(newfilename, model)

        if outputfile:
            pkas_2_pdb(
                subID,
                newfilename,
                f"{dir_path}/pdbs_out/out_{subID}.pdb",
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
            "ncpus": config["SLURM_JOB_NCORES"],
            "clean": True,
            "keep_ions": True,
            "ser_thr_titration": False,
            "titration_output": f"{dir_path}/titrations/titration_{subID}.out",
            "output": f"{dir_path}/pkas/pKas_{subID}.out",
        }

        if outputfile:
            parameters["structure_output"] = (
                f"{dir_path}/pdbs_out/out_{subID}.pdb",
                outputfilepH,
                outputfilenaming,
            )

        sub_parameters = (pdbid, pdbfile, outputemail, nsites, nchains)

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


def submit_pypka_job(job_params, sub_params, subID):
    plus_one_pypka_subs()

    pdbid, pdbfile, outputemail, nsites, nchains = sub_params
    cur_date = datetime.datetime.today()
    new_job = Job(dat_time=cur_date, email=outputemail, sub_id=subID)
    db_session.add(new_job)
    db_session.commit()

    pid = None
    if pdbid:
        pid = (
            db_session.query(Protein.protein_id)
            .filter(Protein.pdb_code == pdbid)
            .first()
        )
    else:
        pdbid = f"#{new_job.job_id}"
    if pid:
        pid = pid[0]
    else:
        new_protein = Protein(
            pdb_code=pdbid, pdb_file=pdbfile, nsites=nsites, nchains=nchains
        )
        db_session.add(new_protein)
        db_session.commit()
        pid = new_protein.protein_id

    job_id = new_job.job_id
    db_session.close()
    create_slurm_file(job_params, subID, job_id, pid)


@app.route("/getSubmissions", methods=["POST", "GET"])
def get_submission():
    submission_IDS = (
        db_session.query(Job.job_id, Job.dat_time, Protein.pdb_code, Job.sub_id)
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

    finished = db_session.query(Results.job_id)

    inprogress_IDS = (
        db_session.query(Job.job_id, Job.dat_time, Job.sub_id)
        .filter(Job.job_id.not_in(finished))
        .all()
    )
    queued_ids = [i[:2] + ("QUEUED", i[2]) for i in inprogress_IDS][::-1]

    response = jsonify(queued_ids + sub_ids)
    response.headers.add("Access-Control-Allow-Origin", "*")

    return response


@app.route("/getFile", methods=["GET", "POST"])
@limiter.limit("50 per hour")
def get_file():
    subID = request.json["subID"]
    ftype = request.json["file_type"]

    if ftype == "original_pdb":
        fname = f"{dir_path}/pdbs/{subID}.pdb"
    elif ftype == "pdb_out":
        fname = f"{dir_path}/pdbs_out/out_{subID}.pdb"

    with open(fname) as f:
        content = f.read()

    response = jsonify(content)
    response.headers.add("Access-Control-Allow-Origin", "*")

    return response


PKPDB_PARAMS = {
    "CpHMD_mode": False,
    "LIPIDS": {},
    "bndcon": 3,
    "clean_pdb": True,
    "couple_min": 2.0,
    "cutoff": -1,
    "epsin": 15.0,
    "epssol": 80.0,
    "eqsteps": 1000,
    "ffID": "G54A7",
    "ff_family": "GROMOS",
    "ffinput": "GROMOS",
    "gsize": 81,
    "ionicstr": 0.1,
    "keep_ions": False,
    "maxc": 0.1,
    "mcsteps": 200000,
    "nanoshaper": -1,
    "nlit": 500,
    "nonit": 0,
    "pHmax": 12.0,
    "pHmin": 0.0,
    "pHstep": 0.25,
    "pbc_dim": 0,
    "pbx": False,
    "pby": False,
    "perfil": 0.9,
    "precision": "single",
    "relfac": 0.75,
    "relpar": 0.75,
    "scaleM": 2.0,
    "scaleP": 1,
    "seed": 1234567,
    "ser_thr_titration": False,
    "slice": 0.05,
    "version": "2.1.0",
}

if __name__ == "__main__":
    # app.run(host="127.0.0.1", debug=True)
    app.run(host="0.0.0.0", port="5000")
