from sqlalchemy import Column, Integer, Date, CHAR, JSON, Time, Text, VARCHAR, REAL
from sqlalchemy import ForeignKeyConstraint, UniqueConstraint
from database import Base


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
    dat_time_finish = Column(Date, nullable=False)
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


class UsageStats(Base):
    __tablename__ = "usage_stats"

    usid = Column(Integer, primary_key=True)
    pkpdb_queries = Column(Integer)
    pkpdb_downloads = Column(Integer)
    pypka_subs = Column(Integer)
    pkai_subs = Column(Integer)
