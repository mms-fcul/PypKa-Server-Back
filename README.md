# PypKa-Server-Back

## Dependencies

    python3 -m pip install pypka flask flask_cors
    apt install gawk gcc gfortran libgfortran4

- [Pypka](https://github.com/mms-fcul/PypKa)
- Pypka dependencies
- flask
- flask_cors

### Local Postgresql DB
Follow the instructions found [here](https://www.postgresql.org/download/linux/ubuntu/).
Then run db.sql

### Local Redis DB (no longer necessary)
On Ubuntu18/20:
```
sudo apt install redis-server
sudo sed -i s'/supervised no/supervised systemd/' /etc/redis/redis.conf
sudo systemctl restart redis.service
```

## Run
    python3 server/app.py

## Files
    - app.py: Main flask backend
      database.py  
      models.py
      pka2pI.py
      slurm.py

    - new.py: fastapi result streaming

      

    - OLD: 
      __main__.py
      pkpdb.py
      stream.py
      job_queue.py       