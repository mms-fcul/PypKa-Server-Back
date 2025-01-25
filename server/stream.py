import tornado.ioloop
import tornado.websocket
import subprocess
import os
import sys
from database import db_session
from models import Job, Results, Residue, Pk, Input, Protein
import logging
from pprint import pformat

logging.basicConfig(filename="socket.log", level=logging.DEBUG)

import json

dir_path = os.path.dirname(os.path.realpath(__file__))


class EchoWebSocket(tornado.websocket.WebSocketHandler):
    def open(self):
        print("WebSocket opened")

    def on_message(self, message):
        subID = message

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
                        self.write_message(
                            json.dumps(
                                {
                                    "content": {"failed": True, "log": error},
                                    "subID": subID,
                                    "nsites": nsites,
                                    "nchains": nchains,
                                    "status": "success",
                                    "protein_name": protein_name,
                                    "protein_pdb": pdb_file,
                                }
                            )
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
                        "pdb_out": pdb_out,
                        "outputFilepH": pdb_out_ph,
                        "protein_name": protein_name,
                        "original_pdb": pdb_file,
                    }

                    self.write_message(
                        json.dumps(
                            {
                                "content": response_dict,
                                "subID": subID,
                                "status": "success",
                            }
                        )
                    )
                    return

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

        self.write_message(
            json.dumps(
                {
                    "content": content if content else None,
                    "subID": subID,
                    "status": "running",
                }
            )
        )

    def on_close(self):
        print("WebSocket closed")

    def check_origin(self, origin):
        return True


def make_app():
    return tornado.web.Application(
        [
            (r"/", EchoWebSocket),
        ]
    )


def main(port):
    app = make_app()
    # app.listen(port)
    server = tornado.httpserver.HTTPServer(app)
    server.listen(port)
    # server.bind(port)
    # server.start(0)  # forks one process per cpu
    tornado.ioloop.IOLoop.current().start()


if __name__ == "__main__":
    port = sys.argv[1]
    main(port)
