#File to test the API get function to retrieve genes for
#terms in constants.py
import requests
import urllib.parse
import json
import time

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


def get_GO_genes_API(termlist):
    """
    This function retrieves all genes associated with each term.
    Input of GO terms must be a 1d list of GO Term Accession. e.g. ['GO:1903502','GO:1903508'].
    """
    for term in termlist:
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

def get_ensembl_sequence_API(ID_list):
    """
    This function queries ensembl for nucleotide sequence
    Input of ensembl ID's must be a 1d list. e.g. ['ENSG00000157764']
    """
    for id in ID_list:
        response = requests.get(f"https://rest.ensembl.org/sequence/id/{id}", headers={ "Content-Type" : "text/plain"})
        logging.debug(response.text)

def uniprot_mapping(id_list, target='Ensembl'):
    """
    This function recieves uniprot or other ID and finds the Ensembl id.
    Input of ID's must be a 1d list. e.g. ['UniProtKB:A0A3Q1N508']
    """

    for item in id_list:
        if "UniProtKB" in item:
            source = "UniProtKB_AC-ID"
        id = item.split(':')[1]
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
                return bool(j["results"] or j["failedIds"])
        

#get_GO_genes_API(['GO:0001525'])
#get_ensembl_sequence_API(["ENSG00000157764"])
uniprot_mapping(['UniProtKB:A0A3Q1N508'])