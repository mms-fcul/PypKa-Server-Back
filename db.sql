
CREATE TABLE protein (
    protein_id serial,
    pdb_code character,
    pdb_file json NOT NULL,
    nchains integer NOT NULL,
    nsites integer NOT NULL,
    PRIMARY KEY (protein_id),
    UNIQUE (pdb_code)
);

CREATE TABLE job (
    job_id serial,
    dat_time date NOT NULL,
    email text,
    sub_id text NOT NULL,
    PRIMARY KEY (job_id)
);

CREATE TABLE input (    
    input_id serial,
    job_id integer,
    protein_id integer NOT NULL,
    pb_set json,
    mc_set json,
    pypka_set json,
    PRIMARY KEY (job_id),
    FOREIGN KEY (job_id) REFERENCES job(job_id),
    FOREIGN KEY (protein_id) REFERENCES protein(protein_id)
);

CREATE TABLE results (
    results_id serial,
    job_id integer,
    tit_curve json,
    pdb_out json,
    error text,
    isoelectric_point real,
    pdb_out_ph real,
    PRIMARY KEY (job_id),
    FOREIGN KEY (job_id) REFERENCES job(job_id)
);

CREATE TABLE residue (
    res_id serial,
    protein_id integer NOT NULL,
    res_name character(4) NOT NULL,
    chain character(1) NOT NULL,
    res_number int NOT NULL,
    PRIMARY KEY (res_id),
    UNIQUE (protein_id, chain, res_name, res_number),
    FOREIGN KEY(protein_id) REFERENCES protein(protein_id)
);

CREATE TABLE pk (
    pk_id serial,    
    res_id integer,
    job_id integer,    
    pk real,
    PRIMARY KEY (pk_id),
    UNIQUE (res_id, job_id),
    FOREIGN KEY(res_id) REFERENCES residue(res_id),
    FOREIGN KEY(job_id) REFERENCES job(job_id)
);

create table usage_stats(
    usid serial,
    pkpdb_queries   integer,
    pkpdb_downloads integer,
    pypka_subs      integer,
    pkai_subs       integer,
    PRIMARY KEY (usid)
);
