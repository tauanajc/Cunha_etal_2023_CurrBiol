'''
From a list of BUSCO ids (hard-coded below as metazoa_buscos.csv, no header), query the OrthoDB v10 database to retrieve GO terms and ids, and the ancestral terms associated with each specific GO name. Output is metazoa_buscos_GOterms.csv

$ python buscos_to_GO.py

Created by Bruno de Medeiros and Tauana Cunha on February 2023
'''

import requests
import time
import pandas as pd
import numpy as np
from urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter

session = requests.Session()
retry_strategy = Retry(
    total = 10,
    backoff_factor = 1,
    status_forcelist = [ 500, 502, 503, 504 ])
adapter = HTTPAdapter(max_retries=retry_strategy)
session.mount('http://', adapter)
session.mount('https://', adapter)

with open("../data/metazoa_buscos.csv", "r") as infile:
    ids = [linha.strip() for linha in infile]
#ids = ['13at10374', '92at10374']

orthodb_results = []

for y,id in enumerate(ids): #test with first 10 lines: ids[:10]
    time.sleep(1)
    print(y)
    r = session.get(url='https://data.orthodb.org/v10/group', params= {'id':id})
    for k,i in r.json()['data'].items():
        try:
            for x in i:
                temp_dict = x.copy()
                temp_dict['ontology_type'] = k
                temp_dict['orthodb_id'] = id
                orthodb_results.append(temp_dict) 
        except:
            pass
    
tabela = pd.DataFrame(orthodb_results)

go_ancestors = []
go_obsolete = []

#print(tabela['name'])

for z,termo in enumerate(tabela['name']):
    print(z)
    if not termo is np.nan and termo.startswith('GO:'):
        time.sleep(1)
        requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/" + termo + "/ancestors?relations=is_a%2Cpart_of%2Coccurs_in%2Cregulates"
        r = session.get(url=requestURL, headers={ "Accept" : "application/json"})
        if r.json()['results'][0]['isObsolete']:
            go_obsolete.append(True)
            go_ancestors.append(np.nan)
        else:
            go_obsolete.append(False)
            go_ancestors.append(",".join(r.json()['results'][0]['ancestors']))
    else:
        go_obsolete.append(np.nan)
        go_ancestors.append(np.nan)

tabela['obsolete_term'] = go_obsolete
tabela['ancestor_terms'] = go_ancestors

tabela.to_csv('../data/metazoa_buscos_GOterms.csv')

