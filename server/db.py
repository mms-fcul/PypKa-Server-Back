from sqlalchemy import create_engine, ForeignKeyConstraint, UniqueConstraint
from sqlalchemy import Column, Integer, Date, CHAR, JSON, Time, Text, VARCHAR, REAL
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
import psycopg2
from dotenv import dotenv_values
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
config = dotenv_values(f"{dir_path}/../.env")

user = config["user"]
password = config["password"]
ip = config["host"]
port = config["port"]
database = config["database"]

db_path = f"postgresql://{user}:{password}@{ip}:{port}/{database}"

db = create_engine(db_path)
Base = declarative_base()

Session = sessionmaker(db)
session = Session()


class Input(Base):
    __tablename__ = "input"

    input_id = Column(Integer, primary_key=True)
    job_id = Column(Integer, nullable=False)
    protein_id = Column(Integer, nullable=False)
    pb_set = Column(JSON, nullable=False)
    mc_set = Column(JSON, nullable=False)
    pypka_set = Column(JSON, nullable=False)
    ForeignKeyConstraint(["job_id"], ["job.job_id"])
    ForeignKeyConstraint(["protein_id"], ["protein.protein_id"])


class Job(Base):
    __tablename__ = "job"

    job_id = Column(Integer, primary_key=True)
    sub_id = Column(Text, nullable=False)
    dat_time = Column(Date, nullable=False)
    email = Column(Text)


class Pk(Base):
    __tablename__ = "pk"

    pk_id = Column(Integer, primary_key=True)
    res_id = Column(Integer)
    job_id = Column(Integer)
    pk = Column(REAL)
    ForeignKeyConstraint(["job_id"], ["job.job_id"])
    ForeignKeyConstraint(["res_id"], ["residue.res_id"])


class Protein(Base):
    __tablename__ = "protein"

    protein_id = Column(Integer, primary_key=True)
    nchains = Column(Integer(), nullable=False)
    nsites = Column(Integer, nullable=False)
    pdb_code = Column(CHAR, unique=True)
    pdb_file = Column(JSON, nullable=False)


class Residue(Base):
    __tablename__ = "residue"

    res_id = Column(Integer, primary_key=True)
    protein_id = Column(Integer, nullable=False)
    res_name = Column(CHAR, nullable=False)
    chain = Column(CHAR, nullable=False)
    res_number = Column(Integer, nullable=False)
    ForeignKeyConstraint(["protein_id"], ["protein.protein_id"])


class Results(Base):
    __tablename__ = "results"

    results_id = Column(Integer, primary_key=True)
    job_id = Column(Integer, nullable=False)
    tit_curve = Column(JSON)
    isoelectric_point = Column(REAL)
    pdb_out = Column(JSON)
    pdb_out_ph = Column(REAL)
    error = Column(Text)
    ForeignKeyConstraint(["job_id"], ["job.job_id"])


# import psycopg2 as pg
# from decouple import config
# import json
# import numpy as np
#
# def db_connect():
#    conn = pg.connect(
#        host=config("host"),
#        database=config("database"),
#        user=config("user"),
#        password=config("password"),
#    )
#    cur = conn.cursor()
#
#    return conn, cur
#
#
# def db_close():
#    conn.close()
#    cur.close()
#
#
# def executeSingleSQLstatement(conn, cur, sql, fetchall=False, commit=True):
#    cur.execute(sql)
#    if commit:
#        conn.commit()
#    if fetchall:
#        return cur.fetchall()
#
#
# def insert_new_submission(
#    conn, cur, cur_date, to_insert, pdbid, params, pdbfile, outputemail, error
# ):
#
#    pKas = to_insert["pKas"]
#    pdb_out = to_insert["pdb_out"]
#    titration = to_insert["titration"]
#
#    pypka_set, pb_set, mc_set = params
#
#    mc_set["pH_values"] = mc_set["pH_values"].tolist()
#
#    res_name = [item[1] for item in pKas]
#    chain = [item[0] for item in pKas]
#    res_number = [item[2] for item in pKas]
#    res_name = np.array(res_name).tolist()
#    chain = np.array(chain).tolist()
#    res_number = np.array(res_number).tolist()
#
#    pK = [x.pop(3) for x in pKas]
#    pK = np.array(pK).tolist()
#
#    protein_id = insert_new_protein(conn, cur, pdbid, pdbfile)
#    job_id = insert_new_job(conn, cur, cur_date, outputemail)
#    input_id = insert_new_input(
#        conn, cur, job_id, protein_id, pypka_set, pb_set, mc_set
#    )
#    results_id = insert_new_results(conn, cur, job_id, titration, pdb_out, error)
#    res_id, pk_id = insert_new_residue_pk(
#        conn, cur, protein_id, res_name, chain, res_number, pK, job_id
#    )
#
#
# def clean_json(old_dic):
#    return json.dumps(old_dic).replace("'", "")
#
#
# def insert_new_protein(conn, cur, pdbid, pdbfile):
#    sql_query = f"""
#    SELECT protein_id
#    FROM protein
#    WHERE pdb_code = '{pdbid}';
#    """
#    pid = executeSingleSQLstatement(conn, cur, sql_query, fetchall=True)
#    if pid:
#        return pid[0][0]
#    else:
#        sql_protein = f"""
#        INSERT INTO protein (pdb_code, pdb_file)
#        VALUES ('{pdbid}', '{clean_json(pdbfile)}')
#        RETURNING protein_id;
#        """
#        protein_id = executeSingleSQLstatement(
#            conn, cur, sql_protein, commit=True, fetchall=True
#        )[0][0]
#
#    return protein_id
#
#
# def insert_new_job(conn, cur, dat_time, outputemail):
#    sql_job = f"""
#    INSERT INTO job (dat_time, email)
#    VALUES('{dat_time}', '{outputemail}')
#    RETURNING job_id;
#    """
#    job_id = executeSingleSQLstatement(conn, cur, sql_job, commit=True, fetchall=True)[
#        0
#    ][0]
#    return job_id
#
#
# def insert_new_input(conn, cur, job_id, protein_id, pypka_set, pb_set, mc_set):
#    sql_input = f"""
#    INSERT INTO input (job_id, protein_id, pb_set, mc_set, pypka_set)
#    VALUES('{job_id}','{protein_id}', '{clean_json(pb_set)}','{clean_json(mc_set)}','{clean_json(pypka_set)}')
#    RETURNING input_id;
#    """
#    input_id = executeSingleSQLstatement(
#        conn, cur, sql_input, commit=True, fetchall=True
#    )[0][0]
#    return input_id
#
#
# def insert_new_results(conn, cur, job_id, tit_curve, pdb_out, error):
#    sql_results = f"""
#    INSERT INTO results (job_id, tit_curve, pdb_out, error)
#    VALUES ('{job_id}','{clean_json(tit_curve)}','{clean_json(pdb_out)}', '{error}')
#    RETURNING results_id;
#    """
#    results_id = executeSingleSQLstatement(
#        conn, cur, sql_results, commit=True, fetchall=True
#    )[0][0]
#    return results_id
#
#
# def insert_new_residue_pk(
#    conn, cur, protein_id, res_name, chain, res_number, pK, job_id
# ):
#
#    for name, cha, number, pk in zip(res_name, chain, res_number, pK):
#
#        sql_query = f"""
#        SELECT res_id
#        FROM Residue
#        WHERE protein_id = '{protein_id}' and res_name = '{name}' and chain = '{cha}' and res_number = '{number}'
#        """
#        res_id = executeSingleSQLstatement(conn, cur, sql_query, fetchall=True)
#        if res_id:
#            res_id = res_id[0][0]
#
#        else:
#            sql_residue = f"""
#            INSERT INTO residue (protein_id, res_name, chain, res_number)
#            VALUES ('{protein_id}','{name}','{cha}','{number}')
#            RETURNING res_id;
#            """
#            res_id = executeSingleSQLstatement(
#                conn, cur, sql_residue, commit=True, fetchall=True
#            )[0][0]
#
#        sql_pk = f"""
#        INSERT INTO pk (job_id, res_id, pk)
#        VALUES ('{job_id}','{res_id}', '{pk}')
#        RETURNING pk_id;
#        """
#        pk_id = executeSingleSQLstatement(
#            conn, cur, sql_pk, commit=True, fetchall=True
#        )[0][0]
#    return res_id, pk_id
#