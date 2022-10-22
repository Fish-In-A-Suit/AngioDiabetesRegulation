#File to test the API get function to retrieve genes for
#terms in constants.py
from re import S
from xmlrpc.client import ResponseError
import requests
import urllib.parse
import json
import time
import csv

import logging

console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
file_handler = logging.FileHandler("./log_output/test_json_dump.log", 'w+')
file_handler.setLevel(logging.DEBUG)
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        file_handler,
        console_handler
    ]
)


def get_GO_genes_API(term):
    """
    This function retrieves all genes associated with each term.
    Input of GO terms must be a 1d list of GO Term Accession. e.g. ['GO:1903502','GO:1903508'].
    """
    url_term = urllib.parse.quote(term)
    parameters = {
        "use_compact_associations":True,
        "taxon":["NCBITaxon:9606"] #only for Homo Sapiens
    }
    response = requests.get(f"http://api.geneontology.org/api/bioentity/function/{url_term}/genes", params=parameters)
    logging.debug(json.dumps(response.json()['compact_associations'], indent=4))
    compact_assoc = response.json()['compact_associations']
    for item in compact_assoc:
        if item['subject'] == term: #only use directly associated genes
            genes=item['objects']
            logging.info(f"GO term: {term} -> Genes/products: {genes}")
    return genes

def get_ensembl_sequence_API(id):
    """
    This function queries ensembl for nucleotide sequence
    Input of ensembl ID's must be a 1d list. e.g. ['ENSG00000157764']
    """
    response = requests.get(f"https://rest.ensembl.org/sequence/id/{id}", headers={ "Content-Type" : "text/plain"})
    if response.ok:
        logging.debug(response.text)
        logging.info(f"Recieved sequence for id {id}.")
        return response.text
    else:
        return None

def uniprot_mapping(id_old, target='Ensembl'):
    """
    This function recieves uniprot or other ID and finds the Ensembl id.
    Input of ID's must be a 1d list. e.g. ['UniProtKB:A0A3Q1N508']
"""
    if "UniProtKB" in id_old:
        source = "UniProtKB_AC-ID"
    id = id_old.split(':')[1]
    response = requests.post(f"https://rest.uniprot.org/idmapping/run", data={'from':source, 'to':target, 'ids':id})
    logging.debug(response.json())
    jobId=response.json()['jobId']
    logging.debug(jobId)
    POLLING_INTERVAL=1
    while True:
        request = requests.get(f"https://rest.uniprot.org/idmapping/status/{jobId}")
        j = request.json()
        logging.info(j)
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(j["jobStatus"])
        else:
            if 'failedIds' in j:
                ensembl_id=None
            else:
                ensembl_id=j['results'][0]['to'].split('.')[0]
            break
    logging.info(ensembl_id)
    return ensembl_id
        
        
#DEMO TEST CODE, not to be used in production
f = open("demofile2.txt", "w+")
terms = ['GO:0001525']
for term in terms:
    genes = get_GO_genes_API(term)
    e_id = []
    seqeunces = []
    print(len(genes))
    for i in range(len(genes)):
        e_id.append(uniprot_mapping(genes[i]))
        if e_id[i] == None:
            seqeunces.append(None)
        else:
            seqeunces.append(get_ensembl_sequence_API(e_id[i]))
        out = {"term" : term, "gene" : genes[i], "ensembel_id" : e_id[i], "sequence" : seqeunces[i]}
        f.write(json.dumps(out)+"\n")
f.close()