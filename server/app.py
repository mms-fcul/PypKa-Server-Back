from flask import Flask, jsonify, request, send_file, Response, make_response
from flask_cors import CORS, cross_origin

import hashlib
from pprint import pprint, pformat
import datetime
import os
import smtplib
import requests
from dotenv import dotenv_values
import time
import sys
from contextlib import contextmanager
import logging
from multiprocessing import Process, Queue
from threading import Thread
import traceback
import json


from pypka import Titration
from pypka import __version__ as pypka_version
from pypka.main import getTitrableSites

import db

logging.basicConfig(
    filename="/home/pedror/PypKa-Server-Back/server/server.log", level=logging.DEBUG
)

app = Flask(__name__, static_url_path="/static")
app.config["SECRET_KEY"] = str(datetime.datetime.today())
app.config["CORS_HEADERS"] = "Content-Type"

CORS(app, resources={r"/*": {"origins": "*"}})


def _build_cors_prelight_response():
    response = make_response()
    response.headers.add("Access-Control-Allow-Origin", "*")
    response.headers.add("Access-Control-Allow-Headers", "*")
    response.headers.add("Access-Control-Allow-Methods", "*")
    return response


def _corsify_actual_response(response):
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


dir_path = os.path.dirname(os.path.realpath(__file__))
config = dotenv_values(f"{dir_path}/../.env")

import importlib.util

spec = importlib.util.spec_from_file_location("pkpdb", config["PKPDB"])
pkpdb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pkpdb)

session = pkpdb.session
Protein = pkpdb.Protein
Residue = pkpdb.Residue
Pk_sim = pkpdb.Pk_sim
Pk = pkpdb.Pk

from job_queue import Jobqueue

job_queue = Jobqueue()


@app.route("/isoelectric.csv.gz")
def export_pis():
    return app.send_static_file("isoelectric.csv.gz")


@app.route("/pkas.csv.gz")
def export_pkas():
    return app.send_static_file("pkas.csv.gz")


@app.route("/similarity090.csv.gz")
def export_clusters():
    return app.send_static_file("similarity090.csv.gz")


@app.route("/")
def hello():
    return "It's alive!"


def retrieve_from_pkpdb(idcode):
    pid = session.query(Protein.pid).filter(Protein.idcode == idcode).first()
    results = False
    if pid:
        pks = (
            session.query(
                Residue.chain, Residue.residue_type, Residue.residue_number, Pk.pk
            )
            .join(Protein, Protein.pid == Residue.pid)
            .filter(Residue.resid == Pk.resid)
            .filter(Protein.idcode == idcode)
            .all()
        )
        if pks:
            results = []
            for pk in pks:
                chain, resname, resnumb, pka = pk
                if pka:
                    pka = round(pka, 2)
                else:
                    pka = "-"
                results.append((chain, resname, resnumb, pka))

    return results


def retrieve_pkpdb_titcurve(idcode):
    pid = session.query(Protein.pid).filter(Protein.idcode == idcode).first()
    tit_x = []
    tit_y = []
    if pid:
        tit_curve = session.query(Pk_sim.tit_curve).filter(Pk_sim.pid == pid[0]).all()
        if tit_curve:
            for pH, prot in tit_curve[0][0].items():
                tit_x.append(pH)
                tit_y.append(round(prot, 2))

    return tit_x, tit_y


def retrieve_pkpdb_pis(idcode):
    pid = session.query(Protein.pid).filter(Protein.idcode == idcode).first()
    pi = None
    if pid:
        pis = session.query(Pk_sim.isoelectric_point).filter(Pk_sim.pid == pid[0]).all()
        if pis:
            pi = round(pis[0][0], 2)

    return pi


@app.route("/query/<idcode>")
def pkpdb_query(idcode):

    tit_x, tit_y = retrieve_pkpdb_titcurve(idcode)

    response_dict = {
        "tit_x": tit_x,
        "tit_y": tit_y,
        "pKas": retrieve_from_pkpdb(idcode),
        "pI": retrieve_pkpdb_pis(idcode),
        "params": "DEFAULT PARAMS",
    }

    response = jsonify(response_dict)
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


@app.route("/pkpdb/<idcode>", methods=["GET", "POST"])
def exists_on_pkpdb(idcode):
    results = retrieve_from_pkpdb(idcode)
    if results and len(results) > 0:
        results = True
    response = jsonify(results)
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


@app.route("/run/<idcode>")
def run(idcode):

    subID = get_subID(request)

    r = requests.get(f"https://files.rcsb.org/download/{idcode}.pdb")
    pdb_content = r.content.decode("utf-8")
    if "The requested URL was not found on this server." in pdb_content:
        return f"Error: {idcode} not found"

    newfilename = save_pdb(pdb_content, subID)

    parameters = {
        "structure": newfilename,
        "epsin": 15,
        "ionicstr": 0.1,
        "convergence": 0.1,
        "pbc_dimensions": 0,
        "ncpus": 8,
        "clean": True,
        "ser_thr_titration": False,
        "titration_output": f"{dir_path}/titrations/titration_{subID}.out",
        "output": f"{dir_path}/pkas/pKas_{subID}.out",
    }

    try:
        response_dict = run_pypka(parameters, subID)
    except Exception as e:
        response_dict = {"PDBID": idcode, "Error": str(e)}
        return jsonify(response_dict)

    response_dict = {"pKs": response_dict["pKas"]}
    response_dict["PDBID"] = idcode
    response = jsonify(response_dict)
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


def fileno(file_or_fd):
    fd = getattr(file_or_fd, "fileno", lambda: file_or_fd)()
    if not isinstance(fd, int):
        raise ValueError("Expected a file (`.fileno()`) or a file descriptor")
    return fd


@contextmanager
def stdout_redirected(to=os.devnull, stdout=None):
    if stdout is None:
        stdout = sys.stdout

    stdout_fd = fileno(stdout)
    with os.fdopen(os.dup(stdout_fd), "wb") as copied:
        stdout.flush()
        try:
            os.dup2(fileno(to), stdout_fd)
        except ValueError:
            with open(to, "wb") as to_file:
                os.dup2(to_file.fileno(), stdout_fd)
        try:
            yield stdout
        finally:
            stdout.flush()
            os.dup2(copied.fileno(), stdout_fd)


@app.route("/queue-size")
def get_queue_size():
    response = jsonify(len(job_queue))
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


def launch_pypka_process(subID, params, queue):
    new_dir = f"/tmp/tmp_{subID}"
    os.mkdir(new_dir)
    os.chdir(new_dir)
    try:
        with open(f"LOG_{subID}", "a") as f, stdout_redirected(f):
            tit = Titration(params)

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

        if type(pI) == tuple:
            pI, limit, charge = pI

        queue.put((tit_x, tit_y, pKs, params, params_dicts, pI))

    except Exception as e:
        queue.put(e)


def run_pypka(parameters, subID, get_params=False):
    """"""
    logging.info(f"{subID} added to the Queue")
    job_queue.add(subID)

    while job_queue.top() != subID:
        time.sleep(2)
        logging.info(f"{job_queue.top()} {subID}")

    logging.info(f"{subID} started")

    q = Queue()
    p = Process(
        target=launch_pypka_process,
        args=(
            subID,
            parameters,
            q,
        ),
    )
    p.start()
    p.join()
    results = q.get()

    logging.info(f"{subID} finished succesfully: {isinstance(results, Exception)}")
    logging.info(f"{subID} removed from the Queue")
    job_queue.pop()

    if isinstance(results, Exception):
        raise results
    else:
        tit_x, tit_y, pKs, params, params_dicts, pI = results

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

    subID = get_subID(request)

    newfilename = save_pdb(pdbfile, subID)

    if outputpKs:
        pH = f"{pHmin},{pHmax}"
    else:
        pH = str(outputfilepH)

    if on_PKPDB and default_parameters and not outputfile:
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
            "ncpus": 16,
            "clean": True,
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


def submit_pypka_job(pypka_params, sub_params, subID):
    pdbid, pdbfile, outputemail, nsites, nchains = sub_params
    try:
        response_dict, final_params = run_pypka(pypka_params, subID, get_params=True)
        response_dict["error"] = None
    except Exception as e:
        tb = traceback.format_exc()
        logging.error(tb)
        response_dict = {"error": str(e)}

    pdb_out = None
    if "structure_output" in pypka_params and not response_dict["error"]:
        with open(f"{dir_path}/pdbs_out/out_{subID}.pdb") as f:
            pdb_out = f.read()
    response_dict["pdb_out"] = pdb_out

    cur_date = datetime.datetime.today()

    session = db.session
    Job = db.Job
    Protein = db.Protein
    Input = db.Input
    Residue = db.Residue
    Results = db.Results
    Pk = db.Pk

    new_job = Job(dat_time=cur_date, email=outputemail, sub_id=subID)
    session.add(new_job)
    session.commit()

    pid = session.query(Protein.protein_id).filter(Protein.pdb_code == pdbid).first()
    if pid:
        pid = pid[0]
    else:
        new_protein = Protein(
            pdb_code=pdbid, pdb_file=pdbfile, nsites=nsites, nchains=nchains
        )
        session.add(new_protein)
        session.commit()
        pid = new_protein.protein_id

    if response_dict["error"]:
        new_results = Results(
            job_id=new_job.job_id,
            error=response_dict["error"],
        )
        session.add(new_results)
        session.commit()

        new_input = Input(job_id=new_job.job_id, protein_id=pid)
        session.add(new_input)
        session.commit()

        return

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
    mc_params["pH_values"] = list(mc_params["pH_values"])

    new_input = Input(
        job_id=new_job.job_id,
        protein_id=pid,
        pypka_set=pypka_params,
        pb_set=delphi_params,
        mc_set=mc_params,
    )
    session.add(new_input)
    session.commit()

    for (chain, res_name, res_number, pK) in response_dict["pKas"]:
        resid = (
            session.query(Residue.res_id)
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
            session.add(new_residue)
            session.commit()
            resid = new_residue.res_id

        new_pk = Pk(job_id=new_job.job_id, res_id=resid, pk=pK)
        session.add(new_pk)
        session.commit()

    new_results = Results(
        job_id=new_job.job_id,
        tit_curve=response_dict["titration"],
        isoelectric_point=response_dict["pI"],
        pdb_out=response_dict["pdb_out"],
    )
    session.add(new_results)
    session.commit()

    # if outputemail:
    #    send_email(outputemail)

    logging.info(f"{subID} exiting")


@app.route("/getSubmissions", methods=["POST", "GET"])
def get_submission():

    session = db.session
    Job = db.Job
    Protein = db.Protein
    Input = db.Input

    submission_IDS = (
        session.query(Job.job_id, Job.dat_time, Protein.pdb_code, Job.sub_id)
        .join(Input, Input.protein_id == Protein.protein_id)
        .filter(Input.job_id == Job.job_id)
        .order_by(Job.sub_id.desc())
        .limit(25)
        .all()
    )

    response = jsonify([i[:] for i in submission_IDS])
    response.headers.add("Access-Control-Allow-Origin", "*")

    # TODO: get job_id, protein_name, submission_datetime, pdb_out exists? bolean
    # TODO: fix 1qyv insertion error

    return response


@app.route("/getFile", methods=["POST"])
def get_file(path):
    subID = request.json["subID"]
    ftype = request.json["ftype"]

    # TODO: get correct variable from the database

    # TODO: format the variable accordingly to download as file
    if ftype == "titration":
        pass
    elif ftype == "parameters":
        pass
    elif ftype == "pKas":
        pass
    elif ftype == "mdpdb":
        pass

    fname = "test_file"
    return send_file(fname, cache_timeout=36000)


if __name__ == "__main__":
    app.run(host="127.0.0.1", debug=True)
