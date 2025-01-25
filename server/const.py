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

STATUS = "live"

import os

DIR_PATH = os.path.dirname(os.path.realpath(__file__))

from dotenv import dotenv_values


CONFIG = dotenv_values(f"{DIR_PATH}/../.env")
