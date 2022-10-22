#File to test the API get function to retrieve genes for
#terms in constants.py
import requests
import urllib.parse
import json

import logging
logging.basicConfig(filename='./log_output/test_json_dump.log', filemode='w', level=logging.DEBUG)


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

get_GO_genes_API(["GO:0001525"])