import json
import requests
from pprint import pformat

from flask import jsonify, request, Blueprint

from database import DB_SESSION
from models import UsageStats
from const import PKPDB_PARAMS, STATUS, DIR_PATH
from utils import (
    plus_one_pkpdb_queries,
    plus_one_pkai_subs,
    run_pKAI,
    get_subID,
    save_pdb,
    df_pkpdb_pkas,
    df_pkpdb_pI,
    df_pkpdb_titcurves,
)

direct_routes_bp = Blueprint("direct_routes", __name__)
direct_routes_severelimit_bp = Blueprint("direct_routes_severelimit", __name__)


@direct_routes_bp.route("/stats")
def get_stats():
    stats = DB_SESSION.query(
        UsageStats.pkpdb_queries,
        UsageStats.pkpdb_downloads,
        UsageStats.pypka_subs,
        UsageStats.pKAI_subs,
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


@direct_routes_bp.route("/test")
def test():
    return jsonify({"status": "live"})


@direct_routes_bp.route("/")
def hello():
    status = "live"
    if STATUS:
        status = STATUS
    return jsonify(
        {
            "status": status,
            "endpoints": """
        /pkas/<idcode>     GET/POST    Retrieves the results from the pKPDB if available, and runs a pKAI calculation otherwise
                           Example: https://api.pypka.org/pkas/4LZT
        
        /pkpdb/<idcode>    GET/POST    Queries the pKPDB database
                           Example: https://api.pypka.org/pkpdb/4LZT

        /pKAI/<idcode>     GET/POST    Runs a pKAI calculation
                           Example: https://api.pypka.org/pkpdb/4LZT
        """,
        }
    )


@direct_routes_severelimit_bp.route("/query/<idcode>")
def pkpdb_query(idcode):
    if idcode == "CRON_JOB":
        idcode = "4lzt"
    else:
        plus_one_pkpdb_queries()

    idcode = idcode.lower()

    df_protein_results = df_pkpdb_pkas.query(f"idcode == '{idcode}'").sort_values(
        ["chain", "residue_number"]
    )
    if len(df_protein_results) > 0:
        results_list = df_protein_results[
            ["chain", "residue_name", "residue_number", "pk"]
        ].values.tolist()

        tit_curve = json.loads(
            df_pkpdb_titcurves.query(f"idcode == '{idcode}'").values.tolist()[0][2]
        )
        tit_x, tit_y = list(tit_curve.keys()), list(tit_curve.values())

        response_dict = {
            "tit_x": tit_x,
            "tit_y": tit_y,
            "pI": df_pkpdb_pI.query(f"idcode == '{idcode}'").values.tolist()[0][2],
            "pKas": results_list,
            "params": pformat(PKPDB_PARAMS),
        }

    else:
        response_dict = {}

    response = jsonify(response_dict)
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


@direct_routes_severelimit_bp.route("/pkpdb/<idcode>", methods=["GET", "POST"])
def exists_on_pkpdb(idcode):
    results = df_pkpdb_pkas.query(f"idcode == '{idcode.lower()}'")
    if len(results) > 0:
        results = True
    else:
        results = False
    response = jsonify(results)
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


@direct_routes_severelimit_bp.route("/pKAI/<idcode>", methods=["GET", "POST"])
def run_pKAI_idcode(idcode):
    plus_one_pkai_subs()
    subID = get_subID(request)

    r = requests.get(f"https://files.rcsb.org/download/{idcode}.pdb")
    pdb_content = r.content.decode("utf-8")
    if "The requested URL was not found on this server." in pdb_content:
        results = f"Error: {idcode} not found"

    else:
        newfilename = save_pdb(pdb_content, subID)

        results = run_pKAI(newfilename, "pKAI")

    response = jsonify(results)
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


@direct_routes_severelimit_bp.route("/pkas/<idcode>", methods=["GET", "POST"])
def get_pkas_from_idcode(idcode):
    idcode = idcode.lower()
    subID = get_subID(request)

    df_protein_results = df_pkpdb_pkas.query(f"idcode == '{idcode}'").sort_values(
        ["chain", "residue_number"]
    )
    if len(df_protein_results) > 0:
        results_list = df_protein_results[
            ["chain", "residue_name", "residue_number", "pk"]
        ].values.tolist()
        response_dict = {
            "idcode": idcode.upper(),
            "method": "PypKa (pKPDB)",
            "pdb": f"https://files.rcsb.org/download/{idcode}.pdb",
            "pI": df_pkpdb_pI.query(f"idcode == '{idcode}'").values.tolist()[0][2],
            "pKas": results_list,
            "tit_curve": {
                key: round(value, 2)
                for key, value in json.loads(
                    df_pkpdb_titcurves.query(f"idcode == '{idcode}'").values.tolist()[
                        0
                    ][2]
                ).items()
            },
        }
        return jsonify(response_dict)

    # PDB IDCODE
    if len(idcode) == 4:
        r = requests.get(f"https://files.rcsb.org/download/{idcode}.pdb")
        pdb_content = r.content.decode("utf-8")
    else:
        idcode_upper = idcode.upper()
        r = requests.get(
            f"https://alphafold.ebi.ac.uk/files/AF-{idcode_upper}-F1-model_v4.pdb"
        )
        pdb_content = r.content.decode("utf-8")

    if "ATOM " not in pdb_content:
        return f"Error: {idcode} not found"

    newfilename = save_pdb(pdb_content, subID)

    try:
        # Run pKAI
        results = run_pKAI(newfilename, "pKAI")
        response_dict = {
            "idcode": idcode.upper(),
            "method": "pKAI 1.2.0",
            "pdb": f"https://alphafold.ebi.ac.uk/files/AF-{idcode_upper}-F1-model_v4.pdb",
            **results,
            "tit_curve": {x: y for x, y in zip(results["tit_x"], results["tit_y"])},
        }
        del response_dict["tit_x"]
        del response_dict["tit_y"]
    except Exception as e:
        response_dict = {"idcode": idcode, "Error": str(e)}
        return jsonify(response_dict)

    response = jsonify(response_dict)
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


@direct_routes_bp.route("/queue-size")
def get_queue_size():
    # response = jsonify(len(job_queue))
    response = jsonify(0)
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


@direct_routes_severelimit_bp.route("/getFile", methods=["GET", "POST"])
def get_file():
    subID = request.json["subID"]
    ftype = request.json["file_type"]

    if ftype == "original_pdb":
        fname = f"{DIR_PATH}/pdbs/{subID}.pdb"
    elif ftype == "pdb_out":
        fname = f"{DIR_PATH}/pdbs_out/out_{subID}.pdb"

    with open(fname) as f:
        content = f.read()

    response = jsonify(content)
    response.headers.add("Access-Control-Allow-Origin", "*")

    return response
