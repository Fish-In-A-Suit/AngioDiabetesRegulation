#File to test the API get function to retrieve genes for
#terms in constants.py
import requests
import urllib.parse
import json
import time
import csv
import util
import os
import sys
import constants

import logging
logger = logging.getLogger(__name__)


def get_GO_genes_API(term, taxon="NCBITaxon:9606"):
    """
    Retrieves all genes associated with each term.
    Input of GO terms must be a 1d list of GO Term Accession. e.g. ['GO:1903502','GO:1903508'].
    Homo sapiens taxon is NCBITaxon:9606
    """
    logger.info("get_GO_genes_API: term = " + term)
    parameters = {
        "rows":100000
    }

    # Get JSON response for current term, read 'objects' property (array of genes) into 'genes' array
    response = requests.get(f"http://api.geneontology.org/api/bioentity/function/{term}/genes", params=parameters)
    logger.debug(json.dumps(response.json(), indent=4))
    associations = response.json()['associations']
    genes = []
    for item in associations:
        if item['subject']['taxon']['id'] == taxon and item['object']['id'] == term: #only use directly associated genes
            genes.append(item['subject']['id'])
    # IMPORTANT: Some terms (like GO:1903587) return only genes related to "subterms" (when calling http://api.geneontology.org:80 "GET /api/bioentity/function/GO%3A1903587/genes?use_compact_associations=True&taxon=NCBITaxon%3A9606 HTTP/1.1" 200 1910)
    # --> no genes associated to the term, only to subterms --> genes array can be of 0 length (and that is not an error)
    logger.info(f"Term {term} has {len(genes)} associated genes/product -> {genes}.")
    return genes

def get_ensembl_sequence_API(id):
    """
    This function queries ensembl for nucleotide sequence
    Input of ensembl ID's must be a string
    """
    response = requests.get(f"https://rest.ensembl.org/sequence/id/{id}?object_type=transcript;type=cds", headers={ "Content-Type" : "text/plain", }) # cds = c-DNA without the UTR regions; type=cdna (in this case UTR region is also kept); retrieves complementary sequence to the gene mRNA (without UTR), same as miRNA sequence (todo: preveri z genetikom)
    if response.ok:
        logger.debug(response.text)
        logger.info(f"Recieved sequence for id {id}.")
        return response.text
    else:
        return None

def get_rnacentral_sequence_API(id):
    """
    This function queries RNA Central for nucleotide sequence
    Input of RNA Central ID's must be a string
    """
    response = requests.get(f"http://rnacentral.org/api/v1/rna/{id}/?format=json")
    if response.ok:
        logger.debug(response.json())
        sequence = response.json()['sequence']
        logger.info(f"Recieved sequence for id {id} -> {sequence}.")
        return sequence
    else:
        logger.debug(f"RNACentral API error")
        return None

def uniprot_mapping_API(id_old, source='UniProtKB_AC-ID', target='Ensembl_Transcript'): # !
    """
    Recieves uniprot or other ID and finds the Ensembl id.
    Input of ID's must be a string
    https://www.uniprot.org/help/id_mapping
    """
    logger.debug(f"[uniprot_mapping]: id_old = {id_old}")
    source = "UniProtKB_AC-ID" # try to map all genes from other databases to UniProt

    # TODO: some terms (like GO-1903670) have genes that are not defined in UniProt! For example, one gene from
    # GO-1903670 has id_old ZFIN:ZDB-GENE-041014-357, throws UnboundLocalError (source referenced before assignment)
    # -> need to search through multiple databases or do we only limit uniprot?
    # possible solution: source = ""
    # but this omits any databases that are not uniprot

    if "UniProtKB" in id_old:
        source = "UniProtKB_AC-ID"
    if "ZFIN" in id_old: # gene is from Zebrafish Informatics Network -> check if a human gene ortholog exists in zfin_human-gene-orthologs.txt -> attempt to find its uniprot id
        human_gene_symbol = util.zfin_find_human_ortholog(id_old) # eg. adgrg9
        if "ZfinError" in human_gene_symbol:
            logger.debug(f"[uniprot_mapping]: ERROR! human_gene_symbol for {id_old} was not found!")
            input("Press enter to proceed.")
        else: #human ortholog exists in uniprot
            id_old = util.get_uniprotId_from_geneName(human_gene_symbol, trust_genes=False)
            logger.debug(f"id_old = {id_old}")
            if "CycleOutOfBoundsError" in id_old or id_old==0:
                return None

    id = id_old.split(':')[1]
    response = requests.post(f"https://rest.uniprot.org/idmapping/run", data={'from':source, 'to':target, 'ids':id})
    logger.debug(response.json())
    jobId=response.json()['jobId']
    logger.debug(jobId)
    POLLING_INTERVAL=1
    while True:
        request = requests.get(f"https://rest.uniprot.org/idmapping/status/{jobId}")
        j = request.json()
        logger.info(j)
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                logger.debug(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(j["jobStatus"])
        else:
            if 'failedIds' in j:
                ensembl_id=None
            else:
                ensembl_id=j['results'][0]['to'].split('.')[0]
            break
    logger.info(ensembl_id)
    return ensembl_id

def find_genes_related_to_GO_terms(terms, ask_for_overrides=True, destination_folder="term_genes"):
    """
    Finds the genes related to the terms array and dumps the results into a json file.
    """
    logger.info(f"terms array len:{len(terms)}, elements: {terms}")
    for term in terms:
        term_file = str(term).replace(":","-")
        filepath=f"{destination_folder}/{term_file}.json"
        _find_genes_related_to_GO_term(term, filepath, ask_for_overrides)


def _find_genes_related_to_GO_term(term, filepath, ask_for_overrides):
    """
    Finds the genes related to the term and dumps the results into a json file.
    """
    logger.debug(f"started gene search for GO term {term}")
    override = 0
    if os.path.isfile(filepath) and ask_for_overrides == True:
        override = input(f"File {filepath} already exists. Enter 1 to process the file again or 0 to skip:")
        if int(override) == 0: # careful! override changes type to str as input is given
            logger.info(f"Skipping file {filepath}")
            return

    genes = get_GO_genes_API(term) # get array of genes associated to a term
    e_id = [] #ensemble id
    seqeunces = []
    json_dictionaries = []
    for gene in genes:
        if 'UniProtKB' in gene:
            e_id.append(uniprot_mapping_API(gene))
            seqeunces.append(get_ensembl_sequence_API(e_id[-1]))
        elif 'ZFIN' in gene:
            human_gene_symbol = util.zfin_find_human_ortholog(gene) # eg. adgrg9
            if "ZfinError" in human_gene_symbol:
                logger.debug(f"[uniprot_mapping]: ERROR! human_gene_symbol for {gene} was not found!")
                input("Press enter to proceed.")
                e_id.append(None)
                seqeunces.append(None)
            else: #human ortholog exists in uniprot
                e_id.append(util.get_uniprotId_from_geneName(human_gene_symbol, trust_genes=False))
                logger.debug(f"id_old = {e_id[-1]}")
                if "CycleOutOfBoundsError" in e_id[-1] or e_id[-1]==0:
                    e_id[-1]=None
                    seqeunces.append(None)
                else:
                    seqeunces.append(get_ensembl_sequence_API(e_id[-1]))
        elif 'RNAcentral' in gene:
            e_id.append(gene.split(':')[1])
            seqeunces.append(get_rnacentral_sequence_API(e_id[-1]))
        else:
            e_id.append(None)
            seqeunces.append(None)
        out = {"term" : term, "product" : gene, "sequence_id" : e_id[-1], "sequence" : seqeunces[-1]}
        # f.write(json.dumps(out)+"\n") For file decoding purposes, json dicts need to be stored in a list and then the list written to the file as per https://stackoverflow.com/questions/21058935/python-json-loads-shows-valueerror-extra-data
        json_dictionaries.append(out)
    
    file = open(filepath, "w+")
    file.write(json.dumps(json_dictionaries)+"\n")
    file.close()
    logger.debug(f"finished gene search for GO term {term}")

def score_genes(json_files):
    """
    Counts the number of appearances of all the genes across all specified json_files (which contain genes
    related to specific GO terms and are made by the find_genes_related_to_GO_terms function)
    """       
"""
if __name__ == "__main__":
    # TODO: load src_data_files/trusted_genes.txt into constants.TRUSTED_GENES and implement checking to avoid asking user for input on genes already trusted
    util.load_trusted_genes("src_data_files/trusted_genes.txt")
    logger.debug(f"constants.TRUSTED_GENES length: {len(constants.TRUSTED_GENES)}")
    terms_test = ['GO:0001525']
    terms_angiogenesis_ids = util.get_array_terms("ANGIOGENESIS")
    find_genes_related_to_GO_terms(terms_angiogenesis_ids)
    # call score_genes(...) here
"""

"""
Seznam termov -> find_genes_related_to_GO_terms (za vse) -> shranjeno v term_genes
Loop cez vse jsone (prek termov) -> nova mapa za gene (angiogenesis, diabetes) - PROBLEM S PREKRIVANJI GENOV
..> ena mapa za gene, vsak gen ma number_angio pa number_diabetes glede na to pr kokih termih se za vsako pojav + mozno total score
skupn file za vse gene -> vsak gen vsebuje number_angio, number_diabetes + tocno dolocene terme k jih vsebuje

///
termi v povezavi z geni... zdej pa gene v povezavi s termi -> da je za vsak gen s kokimi termi je povezan
"""