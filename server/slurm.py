import os
import logging
import traceback
import subprocess
from time import sleep
import datetime
from pypka import Titration
from pypka import __version__ as pypka_version
from models import Residue, Results, Pk, Input, Job
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


def report_error(job_id, error_msg, pid):
    logging.info(f"LOGGING ERROR of jobid #{job_id}")
    with db_session() as session:
        new_results = Results(
            job_id=job_id,
            error=error_msg,
        )
        session.add(new_results)
        session.commit()
        new_input = Input(job_id=job_id, protein_id=pid)
        session.add(new_input)
        session.commit()


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
        report_error(job_id, response_dict["error"], pid)
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

    job = db_session.query(Job).filter(Job.job_id == job_id).first()
    job.dat_time_finish = datetime.datetime.today()
    db_session.commit()

    # if outputemail:
    #    send_email(outputemail)

    db_session.close()

    logging.info(f"{subID} exiting")


def create_slurm_file(params, subID, job_id, pid):
    slurm_f = f"{dir_path}/submissions/slurm_{subID}.py"
    with open(slurm_f, "w") as f:
        f.write(
            f"""#! /home/pedror/miniconda3/envs/pypka/bin/python3
import sys
sys.path.append("{dir_path}")
from slurm import run_pypka_job

run_pypka_job({params}, {subID}, {job_id}, {pid})
"""
        )

    submit_job(subID, slurm_f)

    with db_session() as session:
        has_reported = (
            session.query(Results.job_id).filter(Results.job_id == job_id).all()
        )
        logging.info(f"{subID} has reported: {has_reported}")
        if not has_reported:
            error_msg = f"Job cancelled due to time limit"
            report_error(job_id, error_msg, pid)


def check_idle_machines():
    sbrun = subprocess.run(
        "clues status | grep idle | wc -l",
        shell=True,
        capture_output=True,
    )
    n_machines_idling = int(sbrun.stdout.decode("utf-8").strip())
    return n_machines_idling


def submit_job(job_name, job_script, ncores=16, partitions="debug"):
    n_machines_idling = check_idle_machines()
    logging.info(f"MACHINES IDLE: {n_machines_idling}")

    sbrun = subprocess.run(
        """clues status | grep "off    enabled" | awk '{print $1}'""",
        shell=True,
        capture_output=True,
    )
    machines_offline = sbrun.stdout.decode("utf-8").strip().splitlines()

    logging.info(f"MACHINES OFF: {machines_offline}")

    if n_machines_idling == 0 and machines_offline:
        sbrun = subprocess.run(
            f"clues poweron {machines_offline[0]}",
            shell=True,
            capture_output=True,
        )
        logging.info(f"POWERING ON: {machines_offline[0]}")

        while n_machines_idling == 0:
            sleep(5)
            n_machines_idling = check_idle_machines()
        logging.info(f"MACHINES IDLE: {n_machines_idling}")

    cmd = f"sbatch -p {partitions} -N 1 -n {ncores} -t 240 -o {dir_path}/submissions/{job_name}.out -e {dir_path}/submissions/{job_name}.out {job_script}"
    logging.info(cmd)

    # os.system(cmd)
    sbrun = subprocess.run(
        cmd,
        shell=True,
        capture_output=True,
    )
    logging.info(f'SLURM SUBMISSION -> {sbrun.stdout.decode("utf-8").strip()}')

    in_queue = True
    while in_queue:
        sleep(5)
        sbrun = subprocess.run(
            f"/home/pedror/bin/s-id pedror | grep {job_name} | wc -l",
            shell=True,
            capture_output=True,
        )
        in_queue = int(sbrun.stdout.decode("utf-8").strip())

    if not os.path.isfile(f"{dir_path}/submissions/{job_name}.out"):
        logging.info(f"RESUBMITTING {job_name}")
        submit_job(job_name, job_script, ncores=16, partitions="debug")
