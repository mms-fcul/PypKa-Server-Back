import psycopg2 as pg

def db_connect():
    conn = pg.connect(host="localhost",database="pypkaserver", user="logger", password="t3lmo", port=63333)
    cur = conn.cursor()

    return conn, cur

def db_close():
    conn.close()
    cur.close()

def executeSingleSQLstatement(sql, fetchall=False, commit=False):
    if not fetchall:
        cur.execute(sql)
        results = None
    else:
        results = cur.execute(sql).fetchall()
    if commit:
        conn.commit()
    return results

def insert_new_submission(conn, cur, cur_date, to_insert):
    # unpack to_insert to variables email
    sql_query = f"""
INSERT INTO Job (dat_time, email)
VALUES ({cur_date}, {email})
RETURNING job_id;
    """
    jobID = executeSingleSQLstatement(sql_query)[0]

    sql_query = f"""
INSERT INTO Protein(job_id, pdb_code, pdb_file, Tit_curve, pdb_out)
VALUES ({jobID},{pdb_code},{pdb_file},{titration}, {pdb_out})
RETURNING protein_id;
"""

    sql_query = f"""
INSERT INTO Residues(protein_id, pka, res_name)
VALUES ({proteinID}, {pKs}, {res_name})
"""

    sql_query = f"""
INSERT INTO Set(job_id, pb_set, mc_set, pypka_set)
VALUES ({jobID}, {pb_set}, {mc_set}, {pypka_set});
"""
#perguntar ao pedro onde estão mc_set  as variaveis

