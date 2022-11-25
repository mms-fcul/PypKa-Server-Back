import os
import asyncio
from fastapi import FastAPI, Request
from sse_starlette.sse import EventSourceResponse
from fastapi.middleware.cors import CORSMiddleware
import subprocess
from database import db_session
from models import Job, Results, Residue, Pk, Input, Protein
from pprint import pformat
import json

app = FastAPI()

origins = [
    "http://pypka.org",
    "https://pypka.org",
    "http://localhost",
    "https://localhost",
    "http://localhost:8000",
    "https://localhost:8000",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


STREAM_DELAY = 1  # second

dir_path = os.path.dirname(os.path.realpath(__file__))


@app.get("/")
async def home():
    return {"message": "ALIVE!"}


@app.get("/stream")
async def message_stream(request: Request, subID: str):
    """Stream subID log"""

    def check_logs(subID):
        # subID = str(subID)
        print(f"RECEIVED {subID}")

        with db_session() as session:
            print(f"STARTED {subID}")
            job_id = session.query(Job.job_id).filter_by(sub_id=subID).first()

            sbrun = subprocess.run(
                f"s-id pedror | grep {subID} | wc -l", shell=True, capture_output=True
            )
            in_queue = int(sbrun.stdout.decode("utf-8").strip())
            print(job_id, in_queue, job_id and not in_queue)
            if job_id and not in_queue:
                job_id = job_id[0]

                tit_x, tit_y, pI = None, None, None
                pKas = []
                pdb_out = None
                nchains, nsites = None, None

                results = (
                    session.query(
                        Results.tit_curve,
                        Results.isoelectric_point,
                        Results.pdb_out,
                        Results.pdb_out_ph,
                        Results.error,
                    )
                    .filter(Results.job_id == job_id)
                    .first()
                )

                protein = (
                    session.query(
                        Protein.nchains,
                        Protein.nsites,
                        Protein.pdb_code,
                        Protein.pdb_file,
                    )
                    .join(Input, Input.protein_id == Protein.protein_id)
                    .filter(Input.job_id == job_id)
                    .first()
                )

                params = (
                    session.query(Input.mc_set, Input.pb_set, Input.pypka_set)
                    .filter(Input.job_id == job_id)
                    .first()
                )

                print(results)

                if results:
                    tit_curve, pI, pdb_out, pdb_out_ph, error = results
                    nchains, nsites, protein_name, pdb_file = protein

                    if error:
                        return {
                            "content": {"failed": True, "log": error},
                            "subID": subID,
                            "nsites": nsites,
                            "nchains": nchains,
                            "status": "success",
                            "protein_name": protein_name,
                            "protein_pdb": pdb_file,
                        }

                    (tit_x, tit_y) = tit_curve

                    mc_set, pb_set, pypka_set = params

                    pHmin, pHmax = mc_set["pHmin"], mc_set["pHmax"]
                    ionicstr, epssol, epsin = (
                        pb_set["ionicstr"],
                        pb_set["epssol"],
                        pb_set["epsin"],
                    )

                    all_params = pformat({**pypka_set, **pb_set, **mc_set})

                pks = (
                    session.query(
                        Residue.chain, Residue.res_name, Residue.res_number, Pk.pk
                    )
                    .filter(Residue.res_id == Pk.res_id)
                    .filter(Pk.job_id == job_id)
                    .order_by(Residue.chain, Residue.res_number, Residue.res_name)
                    .all()
                )

                if pks:
                    for pk in pks:
                        chain, resname, resnumb, pka = pk
                        if pka:
                            pka = round(pka, 2)
                        else:
                            pka = "-"
                        pKas.append((chain, resname, resnumb, pka))

                if tit_x and tit_y and pKas and pI and not error:
                    response_dict = {
                        "tit_x": tit_x,
                        "tit_y": tit_y,
                        "pKas": pKas,
                        "pI": pI,
                        "nsites": nsites,
                        "nchains": nchains,
                        "pHmin": pHmin,
                        "pHmax": pHmax,
                        "ionicStrength": ionicstr,
                        "proteinDielectric": epsin,
                        "solventDielectric": epssol,
                        "params": all_params,
                        "pdb_out": True if pdb_out else False,
                        "outputFilepH": pdb_out_ph,
                        "protein_name": protein_name,
                        "original_pdb": True if pdb_file else False,
                    }

                    return {
                        "content": response_dict,
                        "subID": subID,
                        "status": "success",
                    }

        logpath = f"{dir_path}/submissions/{subID}.out"
        # logpath = f"/tmp/tmp_{subID}/LOG_{subID}"
        print("showing", logpath)

        content = []
        if os.path.isfile(logpath):
            with open(logpath) as f:
                lines = f.readlines()
                for line in lines:
                    line = line.strip()
                    if line.startswith("PB Runs"):
                        if content[-1].startswith("PB Runs"):
                            content[-1] = line
                        else:
                            content.append(line)
                    elif line.startswith("MC Run"):
                        if content[-1].startswith("MC Run"):
                            content[-1] = line
                        else:
                            content.append(line)
                    else:
                        content.append(line)
                content = "\n".join(content)

        return {
            "content": content if content else None,
            "subID": subID,
            "status": "running",
        }

    async def event_generator(subID):
        while True:
            # If client closes connection, stop sending events
            if await request.is_disconnected():
                break

            # Checks for new messages and return them to client if any
            yield {
                "event": "new_message",
                "id": "message_id",
                "data": json.dumps(check_logs(subID)),
            }

            await asyncio.sleep(STREAM_DELAY)

    return EventSourceResponse(
        event_generator(subID), headers={"Cache-Control": "public, max-age=29"}
    )
