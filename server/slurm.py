import os
import logging
import traceback

from pypka import Titration
from pypka import __version__ as pypka_version
from models import Residue, Results, Pk, Input
from database import db_session

dir_path = os.path.dirname(os.path.realpath(__file__))

logging.basicConfig(
    filename="/home/pedror/PypKa-Server-Back/server/server.log", level=logging.DEBUG
)

root = logging.getLogger("werkzeug")
root.handlers = logging.getLogger().handlers
root.setLevel(logging.INFO)


def run_pypka(parameters, subID, get_params=False):
    """"""
    logging.info(f"{subID} started")

    new_dir = f"/tmp/tmp_{subID}"
    os.mkdir(new_dir)
    os.chdir(new_dir)
    try:
        # with open(f"LOG_{subID}", "a") as f, stdout_redirected(f):
        tit = Titration(parameters)

        tit_x = []
        tit_y = []
        for pH, prot in tit.getTitrationCurve().items():
            tit_x.append(pH)
            tit_y.append(prot)

        pKs = []
        for site in tit:
            pK = site.pK
            if pK:
                pK = round(site.pK, 2)
            else:
                pK = "-"

            res_number = site.getResNumber()

            pKs.append([site.molecule.chain, site.res_name, res_number, pK])

        params = tit.getParameters()
        params_dicts = tit.getParametersDict()
        pI = tit.getIsoelectricPoint()

        if isinstance(pI, tuple):
            pI, limit, charge = pI

    except Exception as e:
        raise e

    logging.info(f"{subID} finished succesfully")

    response_dict = {
        "titration": [tit_x, tit_y],
        "pKas": pKs,
        "parameters": params,
        "pI": round(pI, 2),
    }

    if get_params:
        return response_dict, params_dicts
    else:
        return response_dict


def run_pypka_job(job_params, subID, job_id, pid):
    logging.info(f"{subID} added to the Queue")

    try:
        response_dict, final_params = run_pypka(job_params, subID, get_params=True)
        response_dict["error"] = None
    except Exception as e:
        tb = traceback.format_exc()
        logging.error(tb)
        response_dict = {"error": str(e)}

    pdb_out = None
    if "structure_output" in job_params and not response_dict["error"]:
        with open(f"{dir_path}/pdbs_out/out_{subID}.pdb") as f:
            pdb_out = f.read()
    response_dict["pdb_out"] = pdb_out

    logging.warning(response_dict["error"])

    if response_dict["error"]:
        logging.info(f"LOGGING ERROR of jobid #{job_id}")
        new_results = Results(
            job_id=job_id,
            error=response_dict["error"],
        )
        db_session.add(new_results)
        db_session.commit()
        new_input = Input(job_id=job_id, protein_id=pid)
        db_session.add(new_input)
        db_session.commit()
        return

    logging.info(f"Handling parameters")

    to_keep = [
        "CpHMD_mode",
        "ffID",
        "ff_family",
        "ffinput",
        "clean_pdb",
        "LIPIDS",
        "keep_ions",
        "ser_thr_titration",
        "cutoff",
        "slice",
    ]
    pypka_params, delphi_params, mc_params = final_params
    pypka_params = {key: value for key, value in pypka_params.items() if key in to_keep}
    pypka_params["version"] = pypka_version
    if isinstance(mc_params["pH_values"], float):
        mc_params["pH_values"] = [mc_params["pH_values"]]
    else:
        mc_params["pH_values"] = list(mc_params["pH_values"])

    logging.info(f"Inserting into Input")

    new_input = Input(
        job_id=job_id,
        protein_id=pid,
        pypka_set=pypka_params,
        pb_set=delphi_params,
        mc_set=mc_params,
    )
    db_session.add(new_input)
    db_session.commit()

    logging.info(f'inserting {response_dict["pKas"]} pKa values')

    for (chain, res_name, res_number, pK) in response_dict["pKas"]:
        resid = (
            db_session.query(Residue.res_id)
            .filter(Residue.protein_id == pid)
            .filter(Residue.res_name == res_name)
            .filter(Residue.chain == chain)
            .filter(Residue.res_number == res_number)
            .first()
        )

        if resid:
            resid = resid[0]
        else:
            new_residue = Residue(
                protein_id=pid,
                res_name=res_name,
                chain=chain,
                res_number=res_number,
            )
            db_session.add(new_residue)
            db_session.commit()
            resid = new_residue.res_id

        if pK == "-":
            pK = None

        new_pk = Pk(job_id=job_id, res_id=resid, pk=pK)
        db_session.add(new_pk)
        db_session.commit()

    logging.info(f"inserting Results")
    pdb_out_ph = None
    if "structure_output" in job_params:
        logging.error(f'structure_out -> {job_params["structure_output"]}')
        pdb_out_ph = job_params["structure_output"][1]

    new_results = Results(
        job_id=job_id,
        tit_curve=response_dict["titration"],
        isoelectric_point=response_dict["pI"],
        pdb_out=response_dict["pdb_out"],
        pdb_out_ph=pdb_out_ph,
    )
    db_session.add(new_results)
    db_session.commit()

    # if outputemail:
    #    send_email(outputemail)

    logging.info(f"{subID} exiting")


def create_slurm_file(params, subID, job_id, pid):
    slurm_f = f"{dir_path}/submissions/slurm_{subID}.py"
    with open(slurm_f, "w") as f:
        f.write(
            f"""#! /usr/bin/python
import sys
sys.path.append("{dir_path}")
from slurm import run_pypka_job

run_pypka_job({params}, {subID}, {job_id}, {pid})
"""
        )

    submit_job(subID, slurm_f)


def submit_job(job_name, job_script, ncores=16, partitions="debug"):
    cmd = f"sbatch -p {partitions} -N 1 -n {ncores} -o {job_name}.out -e {job_name}.err {job_script}"
    os.system(cmd)
