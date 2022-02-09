"""deprecated: redis queue"""

import redis


class Jobqueue:
    def __init__(self):
        self.redis_conn = redis.from_url("redis://localhost:6379")
        self.queue = "queue"

    def __len__(self):
        return self.redis_conn.llen(self.queue)

    def add(self, subID):
        return self.redis_conn.rpush(self.queue, subID)

    def top(self):
        elem = self.redis_conn.lindex(self.queue, 0)
        if elem:
            elem = elem.decode()
        return elem

    def pop(self):
        self.redis_conn.lpop("queue")

    def in_queue(self, subID):
        for elem in self.redis_conn.lrange(self.queue, 0, -1):
            if elem == subID:
                return True
        return False

    def all(self):
        return [i.decode("utf-8") for i in self.redis_conn.lrange(self.queue, 0, -1)]


import os
import sys
import logging
from pypka import Titration
import time
from multiprocessing import Process, Queue

logger = logging.getLogger(__name__)
job_queue = Jobqueue()


def run_pypka(parameters, subID, get_params=False):
    """"""
    logger.info(f"{subID} added to the Queue")
    job_queue.add(subID)

    while job_queue.top() != subID:
        time.sleep(2)
        logger.info(f"{job_queue.top()} {subID}")

    logger.info(f"{subID} started")

    q = Queue()
    p = Process(
        target=launch_pypka_process,
        args=(
            subID,
            parameters,
            q,
        ),
    )
    p.start()
    p.join()
    results = q.get()

    logger.info(f"{subID} finished succesfully: {not isinstance(results, Exception)}")
    logger.info(f"{subID} removed from the Queue")
    job_queue.pop()

    logger.info(results)

    if isinstance(results, Exception):
        raise results
    else:
        tit_x, tit_y, pKs, params, params_dicts, pI = results

    response_dict = {
        "titration": [tit_x, tit_y],
        "pKas": pKs,
        "parameters": params,
        "pI": round(pI, 2),
    }

    if get_params:
        return response_dict, params_dicts
    else:
        return response_dict


def launch_pypka_process(subID, params, queue):
    new_dir = f"/tmp/tmp_{subID}"
    os.mkdir(new_dir)
    os.chdir(new_dir)
    try:
        with open(f"LOG_{subID}", "a") as f, stdout_redirected(f):
            tit = Titration(params)

        tit_x = []
        tit_y = []
        for pH, prot in tit.getTitrationCurve().items():
            tit_x.append(pH)
            tit_y.append(prot)

        pKs = []
        for site in tit:
            pK = site.pK
            if pK:
                pK = round(site.pK, 2)
            else:
                pK = "-"

            res_number = site.getResNumber()

            pKs.append([site.molecule.chain, site.res_name, res_number, pK])

        params = tit.getParameters()
        params_dicts = tit.getParametersDict()
        pI = tit.getIsoelectricPoint()

        if type(pI) == tuple:
            pI, limit, charge = pI

        queue.put((tit_x, tit_y, pKs, params, params_dicts, pI))

    except Exception as e:
        queue.put(e)


def fileno(file_or_fd):
    fd = getattr(file_or_fd, "fileno", lambda: file_or_fd)()
    if not isinstance(fd, int):
        raise ValueError("Expected a file (`.fileno()`) or a file descriptor")
    return fd


@contextmanager
def stdout_redirected(to=os.devnull, stdout=None):
    if stdout is None:
        stdout = sys.stdout

    stdout_fd = fileno(stdout)
    with os.fdopen(os.dup(stdout_fd), "wb") as copied:
        stdout.flush()
        try:
            os.dup2(fileno(to), stdout_fd)
        except ValueError:
            with open(to, "wb") as to_file:
                os.dup2(to_file.fileno(), stdout_fd)
        try:
            yield stdout
        finally:
            stdout.flush()
            os.dup2(copied.fileno(), stdout_fd)
