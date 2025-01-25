# PypKa-Server-Back

[![API Status](https://img.shields.io/website-up-down-green-red/https/api.pypka.org.svg)](https://api.pypka.org)
[![License](https://img.shields.io/github/license/mms-fcul/PypKa-Server-Back)](LICENSE)
[![Issues](https://img.shields.io/github/issues/mms-fcul/PypKa-Server-Back)](https://github.com/mms-fcul/PypKa-Server-Back/issues)

## Table of Contents
- [Description](#description)
- [Installation](#installation)
- [Usage](#usage)
- [Files](#files)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## Description
PypKa-Server-Back is a backend server for the PypKa project.

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

## Contributing
Contributions are welcome! Please read the [contributing guidelines](CONTRIBUTING.md) first.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
For questions or issues, please open an issue on GitHub or contact [your-email@example.com](mailto:your-email@example.com).