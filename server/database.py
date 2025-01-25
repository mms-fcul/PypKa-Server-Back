from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from contextlib import contextmanager
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

engine = create_engine(db_path, pool_size=3, max_overflow=0, pool_recycle=30)
db_session = scoped_session(
    sessionmaker(autocommit=False, autoflush=False, bind=engine)
)

Base = declarative_base()
