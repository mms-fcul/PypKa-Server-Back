# PypKa-Server-Back

[![API Status](https://img.shields.io/website-up-down-green-red/https/api.pypka.org.svg)](https://api.pypka.org)
[![License](https://img.shields.io/github/license/mms-fcul/PypKa-Server-Back)](LICENSE)
[![Issues](https://img.shields.io/github/issues/mms-fcul/PypKa-Server-Back)](https://github.com/mms-fcul/PypKa-Server-Back/issues)

## Table of Contents
- [Description](#description)
- [Installation](#installation)
- [Usage](#usage)
- [Files](#files)
- [License](#license)
- [Contact](#contact)

## Description
PypKa-Server-Back is a backend server for the [PypKa website](https://pypka.org). If you use this software in your research, please cite its [paper](https://academic.oup.com/nar/article/52/W1/W294/7645774):

@article{10.1093/nar/gkae255,
    author = {Reis, Pedro B P S and Clevert, Djork-Arné and Machuqueiro, Miguel},
    title = {PypKa server: online pKa predictions and biomolecular structure preparation with precomputed data from PDB and AlphaFold DB},
    journal = {Nucleic Acids Research},
    volume = {52},
    number = {W1},
    pages = {W294-W298},
    year = {2024},
    month = {04},
    issn = {0305-1048},
    doi = {10.1093/nar/gkae255}
}

## Installation
```sh
python3 -m pip install pipenv
pipenv install
```

### Installing PypKa
To install PypKa, follow the instructions in the [PypKa documentation](https://pypka.readthedocs.io/en/latest/getting_started/installation.html).

### Local Postgresql DB
Follow the instructions found [here](https://www.postgresql.org/download/linux/ubuntu/).
Then run `db.sql` to set up the database.


## Usage
To run the server:
```sh
python3 server/app.py  # Starts the main Flask api
python3 server/results_stream.py  # Starts the Fastapi api
```

## Files
- `server/app.py`: Main flask backend
- `server/routes_direct.py`: API endpoints for direct usage
- `server/routes_website.py`: API endpoints for the website
- `server/database.py`: Database connection and query handling
- `server/models.py`: Data models and ORM setup
- `server/const.py`: Constants and configuration settings
- `server/results_stream.py`: Streaming results to clients
- `server/slurm.py`: SLURM job scheduling and management
- `server/results_stream.py`: FastAPI result streaming

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
For questions or issues, please open an issue on GitHub or contact [pedroreis@campus.ul.pt](mailto:pedroreis@campus.ul.pt).