import os
from dotenv import dotenv_values
import importlib.util

dir_path = os.path.dirname(os.path.realpath(__file__))
config = dotenv_values(f"{dir_path}/../.env")

spec = importlib.util.spec_from_file_location("pkpdb", config["PKPDB"])
pkpdb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pkpdb)

pkpdbsession = pkpdb.Session
Protein = pkpdb.Protein
Residue = pkpdb.Residue
Pk = pkpdb.Pk
Pk_sim = pkpdb.Pk_sim


def retrieve_from_pkpdb(idcode):
    with pkpdbsession() as session:
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
    with pkpdbsession() as session:
        pid = session.query(Protein.pid).filter(Protein.idcode == idcode).first()
        tit_x = []
        tit_y = []
        if pid:
            tit_curve = (
                session.query(Pk_sim.tit_curve).filter(Pk_sim.pid == pid[0]).all()
            )
            if tit_curve:
                for pH, prot in tit_curve[0][0].items():
                    tit_x.append(pH)
                    tit_y.append(round(prot, 2))

    return tit_x, tit_y


def retrieve_pkpdb_pis(idcode):
    with pkpdbsession() as session:
        pid = session.query(Protein.pid).filter(Protein.idcode == idcode).first()
        pi = None
        if pid:
            pis = (
                session.query(Pk_sim.isoelectric_point)
                .filter(Pk_sim.pid == pid[0])
                .all()
            )
            if pis:
                pi = round(pis[0][0], 2)

    return pi