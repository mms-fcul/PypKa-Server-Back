from flask import Flask, jsonify, request, send_file
from flask_cors import CORS

from pprint import pformat
import datetime
import os
import smtplib
import requests
from dotenv import dotenv_values
import logging
from threading import Thread
import traceback
from contextlib import contextmanager

from pypka.main import getTitrableSites

from database import db_session, Base, engine
from models import Job, Protein, Input, UsageStats

from slurm import create_slurm_file

# from pkpdb import retrieve_from_pkpdb, retrieve_pkpdb_titcurve, retrieve_pkpdb_pis


Base.metadata.create_all(bind=engine)

STATUS = "live"

logging.basicConfig(
    filename="/home/pedror/PypKa-Server-Back/server/server.log", level=logging.DEBUG
)

root = logging.getLogger("werkzeug")
root.handlers = logging.getLogger().handlers
root.setLevel(logging.INFO)

app = Flask(__name__, static_url_path="/static")
app.config["SECRET_KEY"] = str(datetime.datetime.today())
app.config["CORS_HEADERS"] = "Content-Type"

CORS(app, resources={r"/*": {"origins": "*"}})

dir_path = os.path.dirname(os.path.realpath(__file__))
config = dotenv_values(f"{dir_path}/../.env")


@app.teardown_appcontext
def shutdown_session(exception=None):
    db_session.remove()


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
        UsageStats.pkai_subs,
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


@app.route("/isoelectric.csv.gz")
def export_pis():
    plus_one_pkpdb_downloads()
    return app.send_static_file("isoelectric.csv.gz")


@app.route("/pkas.csv.gz")
def export_pkas():
    plus_one_pkpdb_downloads()
    return app.send_static_file("pkas.csv.gz")


@app.route("/similarity090.csv.gz")
def export_clusters():
    plus_one_pkpdb_downloads()
    return app.send_static_file("similarity090.csv.gz")


@app.route("/")
def hello():
    status = "live"
    if STATUS:
        status = STATUS
    return jsonify({"status": status})


@app.route("/query/<idcode>")
def pkpdb_query(idcode):
    if idcode == "CRON_JOB":
        idcode = "4lzt"
    else:
        plus_one_pkpdb_queries()

    # tit_x, tit_y = retrieve_pkpdb_titcurve(idcode)

    response_dict = {
        # "tit_x": tit_x,
        # "tit_y": tit_y,
        # "pKas": retrieve_from_pkpdb(idcode),
        # "pI": retrieve_pkpdb_pis(idcode),
        "params": pformat(PKPDB_PARAMS),
    }

    response = jsonify(response_dict)
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


@app.route("/pkpdb/<idcode>", methods=["GET", "POST"])
def exists_on_pkpdb(idcode):
    # results = retrieve_from_pkpdb(idcode)
    results = False
    if results and len(results) > 0:
        results = True
    response = jsonify(results)
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
            "ncpus": 32,
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

    # TODO: SUBMIT TO SLURM
    create_slurm_file(job_params, subID, new_job.job_id, pid)

    logging.info(f"{subID} submitted to SLURM")


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

    # TODO: GET SLURM QUEUED
    # queued = job_queue.all()
    # logging.info(queued)
    #
    # inprogress_IDS = (
    #     db_session.query(Job.job_id, Job.dat_time, Job.sub_id)
    #     .filter(Job.sub_id.in_(queued))
    #     .all()
    # )
    # logging.info(inprogress_IDS)
    # queued_ids = [i[:2] + ("QUEUED", i[2]) for i in inprogress_IDS][::-1]

    # response = jsonify(queued_ids + sub_ids)
    response = jsonify(sub_ids)
    response.headers.add("Access-Control-Allow-Origin", "*")

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
    app.run(host="127.0.0.1", debug=True)
