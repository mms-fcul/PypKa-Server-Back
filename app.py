from flask import Flask, jsonify, request, send_file
from flask_cors import CORS
from mailjet_rest import Client
from pprint import pprint, pformat
import datetime
import os

from pypka.pypka import Titration, getTitrableSites
import db

app = Flask(__name__)
CORS(app, resources={r'/*': {'origins': '*'}})
CONN, CUR = db.db_connect()

def save_pdb(pdbfile, subID):
    def new_name(subID):
        return 'pdbs/{0}.pdb'.format(subID)

    newfilename = new_name(subID)
    while os.path.isfile(newfilename):
        subID += '_'
        newfilename = new_name(subID)

    with open(newfilename, 'w') as f_new:
        f_new.write(pdbfile)
    return newfilename

def get_subID(request):
    subID = '{}{}'.format(datetime.datetime.today(), hash(str(request.headers)))
    return subID.replace('-', '').replace(':', '').replace(' ', '').replace('.', '')

@app.route('/getTitrableSitesNumber', methods=['POST'])
def getNumberOfTitratableSites():

    pdbfile = request.json['PDB']

    subID = get_subID(request)
    newfilename = save_pdb(pdbfile, subID)

    out_sites, chains_res = getTitrableSites(newfilename, ser_thr_titration=False)

    os.remove(newfilename)

    nchains = len(chains_res.keys())
    nres = 0
    for chain, sites in chains_res.items():
        nres += len(sites)

    response_dict = {
        'nchains': nchains,
        'nsites': nres
    }

    print(response_dict)

    response = jsonify(response_dict)
    response.headers.add('Access-Control-Allow-Origin', '*')
    return response

@app.route('/getSubID', methods=['POST'])
def getSubID():

    subID = get_subID(request)

    response = jsonify({'subID': subID})
    response.headers.add('Access-Control-Allow-Origin', '*')

    print('getSubID', response)

    return response

@app.route('/submitSim', methods=['POST'])
def submitCalculation():

    pdbfile = request.json['pdb']
    input_naming_scheme = request.json['inputNamingScheme']
    
    pHmin   = request.json['pHmin']
    pHmax   = request.json['pHmax']
    pHstep  = request.json['pHstep']
    
    epsin   = request.json['epsin']
    epsout  = request.json['epsout']
    ionic   = request.json['ionic']

    outputpKs        = request.json['outputpKs']
    outputfile       = request.json['outputfile']
    outputfilenaming = request.json['outputNamingScheme']
    outputfilepH     = request.json['outputFilepH']
    #TODO request email to front-end
    outputemail = request.json['outputemail']

    subID   = request.json['subID']

    newfilename = save_pdb(pdbfile, subID)

    if outputpKs:
        pH = '{0},{1}'.format(pHmin, pHmax)
    else:
        pH = str(outputfilepH)

    parameters = {
        'structure'     : newfilename,
        'pH'            : pH,
        'pHstep'        : pHstep,
        'epsin'         : epsin,
        'epsout'        : epsout,
        'ionicstr'      : ionic,
        'ffinput'       : input_naming_scheme,
        'scaleM'        : 2,
        'convergence'   : 0.1,
        'pbc_dimensions': 0,
        'ncpus'         : -1,
        'clean': True,
        'ser_thr_titration': False,
        'titration_output': 'titrations/titration_{0}.out'.format(subID),
        'output'        : 'pkas/pKas_{0}.out'.format(subID)
    }
    
    if outputfile:
        parameters['structure_output'] = ('pdbs_out/out_{0}.pdb'.format(subID), outputfilepH, outputfilenaming)

    pprint(parameters)

    tit = Titration(parameters)

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
            pK = '-'

        res_number = site.res_number
        if res_number > 2000:
                res_number -= 2000

        pKs.append([site.molecule.chain, site.res_name, res_number, pK])

    pdb_out = None
    if outputfile:
        with open('out_{0}.pdb'.format(subID)) as f:
            pdb_out = f.read()

    response_dict = {
        'titration': [tit_x, tit_y],
        'pKas': pKs,
        'parameters': tit.getParameters(),
        'pdb_out': pdb_out
    }

    response = jsonify(response_dict)
    response.headers.add('Access-Control-Allow-Origin', '*')


    cur_date = datetime.datetime.today()
    with open("subs.txt", 'a') as f_new:
        f_new.write('{0} {1}\n'.format(subID, cur_date))
    with open('submissions/{0}'.format(subID), 'w') as f_new:
        f_new.write(pformat(response_dict))

    db.insert_new_submission(CONN, CUR, cur_date, response_dict)

    #TODO send_email(email)
    def send_email(email)

    api_key = os.environ['MJ_APIKEY_PUBLIC']
    api_secret = os.environ['MJ_APIKEY_PRIVATE']
    mailjet = Client(auth=(api_key, api_secret), version='v3.1')
    data = {
      'Messages': [
        {
            "From": {
                "Email": "$SENDER_EMAIL",
                "Name": "PypKa"
            },
            "To": [
                {
                    "Email": email,
                    "Name": "User"
                }
                ],
                "Subject": "My first Mailjet Email!",
                "TextPart": "Greetings from Mailjet!",
                "HTMLPart": "<h3>Dear passenger 1, welcome to <a href=\"https://www.mailjet.com/\">Mailjet</a>!</h3><br />May the delivery force be with you!"
            }
        ]
    }
    result = mailjet.send.create(data=data)
    
    #temos de criar um email
    #criar um html que va diretamente para os nossos dados
    #criar uma autenticação para o api

    return response

@app.route('/getFile', methods=['POST'])
def get_file(path):
    subID = request.json['subID']
    ftype =  request.json['ftype']

    if ftype == 'titration':
        pass
    elif ftype == 'parameters':
        pass
    elif ftype == 'pKas':
        pass
    elif ftype == 'mdpdb':
        pass

    fname = 'test_file'
    return send_file(fname, cache_timeout=36000)

@app.route('/getLatestsSubmissions', methods=['POST'])
def getLatestsSubmissions():
    
    response_dict = {
    }

    response = jsonify(response_dict)
    response.headers.add('Access-Control-Allow-Origin', '*')

    return response
