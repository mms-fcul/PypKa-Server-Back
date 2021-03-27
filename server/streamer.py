from aiohttp import web
import socketio
import os
import time
from job_queue import Jobqueue
from db import session, Job, Results, Residue, Pk, Input, Protein
import logging
from pprint import pformat

logging.basicConfig(
    filename="/home/pedror/PypKa-Server-Back/server/socket.log", level=logging.DEBUG
)

sio = socketio.AsyncServer(cors_allowed_origins="*")
app = web.Application()
sio.attach(app)

jobqueue = Jobqueue()


async def index(request):
    """Serve the client-side application."""
    return web.Response(text="It's Alive!", content_type="text/html")


@sio.on("get pypka status", namespace="/")
async def pypka_status(sid, subID):

    logging.info(f"STARTED {subID}")
    job_id = session.query(Job.job_id).filter_by(sub_id=subID).first()
    if job_id and not jobqueue.in_queue(subID):
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
                Results.error,
            )
            .filter(Results.job_id == job_id)
            .first()
        )

        protein = (
            session.query(Protein.nchains, Protein.nsites)
            .join(Input, Input.protein_id == Protein.protein_id)
            .filter(Input.job_id == job_id)
            .first()
        )

        params = (
            session.query(Input.mc_set, Input.pb_set, Input.pypka_set)
            .filter(Input.job_id == job_id)
            .first()
        )

        if results:
            tit_curve, pI, pdb_out, error = results
            nchains, nsites = protein

            if error:
                await sio.emit(
                    "pypka finished",
                    {
                        "content": {"failed": True, "log": error},
                        "subID": subID,
                        "nsites": nsites,
                        "nchains": nchains,
                        "status": "success",
                    },
                )
                return

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
            session.query(Residue.chain, Residue.res_name, Residue.res_number, Pk.pk)
            .filter(Residue.res_id == Pk.res_id)
            .filter(Pk.job_id == job_id)
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
                "pdb_out": pdb_out,
            }

            await sio.emit(
                "pypka finished",
                {"content": response_dict, "subID": subID, "status": "success"},
            )

    logpath = f"/tmp/tmp_{subID}/LOG_{subID}"

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

    await sio.emit(
        "pypka status", {"content": content if content else None, "subID": subID}
    )


app.router.add_get("/", index)

if __name__ == "__main__":
    web.run_app(app, port="8888")
    # web.run_app(app)
