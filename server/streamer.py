from aiohttp import web
import socketio
import os
import time
from job_queue import Jobqueue
from db import session, Job, Results, Residue, Pk

sio = socketio.AsyncServer(cors_allowed_origins="*")
app = web.Application()
sio.attach(app)

jobqueue = Jobqueue()


async def index(request):
    """Serve the client-side application."""
    return web.Response(text="It's Alive!", content_type="text/html")


@sio.on("get pypka status", namespace="/")
async def pypka_status(sid, subID):

    job_id = session.query(Job.job_id).filter_by(sub_id=subID).first()
    if job_id and not jobqueue.in_queue(subID):
        job_id = job_id[0]

        tit_x, tit_y, pI, pKas = None, None, None, []

        results = (
            session.query(Results.tit_curve, Results.isoelectric_point, Results.error)
            .filter(Results.job_id == job_id)
            .first()
        )

        if results:
            tit_curve, pI, error = results

            if error:
                await sio.emit(
                    "pypka finished",
                    {
                        "content": {"failed": True, "log": error},
                        "subID": subID,
                        "status": "success",
                    },
                )
                return
            else:
                (tit_x, tit_y) = tit_curve

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
                "pI": pI
                # "params": "DEFAULT PARAMS",
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
    # web.run_app(app, port="8888")
    web.run_app(app)
