#File to test the API get function to retrieve genes for
#terms in constants.py
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
    Retrieves all genes associated with each term.
    Input of GO terms must be a 1d list of GO Term Accession. e.g. ['GO:1903502','GO:1903508'].
    """
    url_term = urllib.parse.quote(term)
    parameters = {
        "use_compact_associations":True,
        "taxon":["NCBITaxon:9606"] #only for Homo Sapiens
    }
    # Get JSON response for current term, read 'objects' property (array of genes) into 'genes' array
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
    Recieves uniprot or other ID and finds the Ensembl id.
    Input of ID's must be a 1d list. e.g. ['UniProtKB:A0A3Q1N508']
    https://www.uniprot.org/help/id_mapping
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
                logging.debug(f"Retrying in {POLLING_INTERVAL}s")
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

def find_genes_related_to_GO_terms(terms, destination_file=""):
    """
    Finds the genes related to the terms array and dumps the results into a json file.
    """
    for term in terms:
        term_file = str(term).replace(":","-")
        filepath=f"term_genes/{term_file}.json"
        file = open(filepath, "w+")

        genes = get_GO_genes_API(term) # get array of genes associated to a term
        e_id = [] #ensemble id
        seqeunces = []
        json_dictionaries = []
        for i in range(len(genes)):
            e_id.append(uniprot_mapping(genes[i])) # convert gene ID to Ensembl id
            if e_id[i] == None:
                seqeunces.append(None)
            else:
                seqeunces.append(get_ensembl_sequence_API(e_id[i]))
            out = {"term" : term, "gene" : genes[i], "ensembel_id" : e_id[i], "sequence" : seqeunces[i]}
           # f.write(json.dumps(out)+"\n") For file decoding purposes, json dicts need to be stored in a list and then the list written to the file as per https://stackoverflow.com/questions/21058935/python-json-loads-shows-valueerror-extra-data
            json_dictionaries.append(out)
        file.write(json.dumps(json_dictionaries)+"\n")
        file.close()

def score_genes(json_files):
    """
    Counts the number of appearances of all the genes across all specified json_files (which contain genes
    related to specific GO terms and are made by the find_genes_related_to_GO_terms function)
    """       
        
#DEMO TEST CODE, not to be used in production
"""
f = open("demofile2.json", "w+")
terms = ['GO:0001525']
for term in terms:
    genes = get_GO_genes_API(term) # get array of genes associated to a term
    e_id = []
    seqeunces = []
    json_dictionaries = []
    for i in range(len(genes)):
        e_id.append(uniprot_mapping(genes[i])) # convert gene ID to Ensembl id
        if e_id[i] == None:
            seqeunces.append(None)
        else:
            seqeunces.append(get_ensembl_sequence_API(e_id[i]))
        out = {"term" : term, "gene" : genes[i], "ensembel_id" : e_id[i], "sequence" : seqeunces[i]}
        # f.write(json.dumps(out)+"\n") For file decoding purposes, json dicts need to be stored in a list and then the list written to the file as per https://stackoverflow.com/questions/21058935/python-json-loads-shows-valueerror-extra-data
        json_dictionaries.append(out)
    f.write(json.dumps(json_dictionaries)+"\n")
f.close()
"""
terms = ['GO:0001525']
find_genes_related_to_GO_terms(terms)
# call score_genes(...) here