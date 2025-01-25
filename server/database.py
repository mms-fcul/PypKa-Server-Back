from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy.ext.declarative import declarative_base

from const import CONFIG

user = CONFIG["user"]
password = CONFIG["password"]
ip = CONFIG["host"]
port = CONFIG["port"]
database = CONFIG["database"]

db_path = f"postgresql://{user}:{password}@{ip}:{port}/{database}"

ENGINE = create_engine(db_path, pool_size=3, max_overflow=0, pool_recycle=30)
DB_SESSION = scoped_session(
    sessionmaker(autocommit=False, autoflush=False, bind=ENGINE)
)

BASE = declarative_base()
