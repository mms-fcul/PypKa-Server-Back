import datetime
import os
import pandas as pd

from pkai.pKAI import pKAI

from pka2pI import pkas_2_titcurve, titcurve_2_pI, exclude_cys
from database import DB_SESSION
from models import Job, Protein, UsageStats
from slurm import create_slurm_file
from const import DIR_PATH


def plus_one_pkpdb_downloads():
    usage = DB_SESSION.query(UsageStats).first()
    usage.pkpdb_downloads += 1
    DB_SESSION.commit()


def plus_one_pkpdb_queries():
    usage = DB_SESSION.query(UsageStats).first()
    usage.pkpdb_queries += 1
    DB_SESSION.commit()


def plus_one_pypka_subs():
    usage = DB_SESSION.query(UsageStats).first()
    usage.pypka_subs += 1
    DB_SESSION.commit()


def plus_one_pkai_subs():
    usage = DB_SESSION.query(UsageStats).first()
    usage.pkai_subs += 1
    DB_SESSION.commit()


def load_static_files():
    df_pkpdb_pkas = pd.read_csv("static/pkas.csv", header=0, sep=";")
    df_pkpdb_pI = pd.read_csv("static/isoelectric.csv", header=0, sep=";")
    df_pkpdb_titcurves = pd.read_csv("static/titrationcurves.csv", header=0, sep=";")
    return df_pkpdb_pkas, df_pkpdb_pI, df_pkpdb_titcurves


def submit_pypka_job(job_params, sub_params, subID):
    plus_one_pypka_subs()

    pdbid, pdbfile, outputemail, nsites, nchains, ip, country, city = sub_params
    cur_date = datetime.datetime.today()
    new_job = Job(
        dat_time=cur_date,
        email=outputemail,
        sub_id=subID,
        ip=ip,
        country=country,
        city=city,
    )
    DB_SESSION.add(new_job)
    DB_SESSION.commit()

    pid = None
    if pdbid:
        pid = (
            DB_SESSION.query(Protein.protein_id)
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
        DB_SESSION.add(new_protein)
        DB_SESSION.commit()
        pid = new_protein.protein_id

    job_id = new_job.job_id
    DB_SESSION.close()
    create_slurm_file(job_params, subID, job_id, pid)


def run_pKAI(pdb, model):
    plus_one_pkai_subs()
    results = pKAI(pdb, model_name=model, device="cpu", threads=2)
    results = [[i[0], i[2], i[1], i[3]] for i in results]
    results = exclude_cys(pdb, results)
    tit_x, tit_y = pkas_2_titcurve(pdb, results)
    pI = titcurve_2_pI(tit_x, tit_y)
    return {"pKas": results, "tit_x": tit_x, "tit_y": tit_y, "pI": round(pI, 2)}


def get_subID(request):
    subID = "{}{}".format(datetime.datetime.today(), hash(str(request.headers)))
    return subID.replace("-", "").replace(":", "").replace(" ", "").replace(".", "")


def save_pdb(pdbfile, subID):
    def new_name(subID):
        return f"{DIR_PATH}/pdbs/{subID}.pdb"

    newfilename = new_name(subID)
    while os.path.isfile(newfilename):
        subID += "_"
        newfilename = new_name(subID)

    with open(newfilename, "w") as f_new:
        f_new.write(pdbfile)
    return newfilename


df_pkpdb_pkas, df_pkpdb_pI, df_pkpdb_titcurves = load_static_files()
