# File to test the API get function to retrieve genes for
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
import atexit
import datetime

import logging
logger = logging.getLogger(__name__)

# global variables
FLAG_HOMOSAPIENS_ONLY = True 
FLAG_TRUST_GENES = True # trust genes in trusted_genes to be credible -> program doesnt stop at them for validation
APPROVED_DATABASES = [["UniProtKB", ["NCBITaxon:9606"]],
                      ["ZFIN", ["NCBITaxon:7955"]],
                      ["RNAcentral",["NCBITaxon:9606"]],
                      ["Xenbase",["NCBITaxon:8364"]],
                      ["MGI", ["NCBITaxon:10090"]],
                      ["RGD", ["NCBITaxon:10116"]]] #Which pairs of databases and taxons (could be multiple per database) to allow.

json_dictionaries = []
current_filepath = ""

# Remarks
#   - no need to track _current_file, since .json files in term_genes aren't saved if script terminates early / if all the json elements haven't already been computed
#   - RGD (Rat Genome Database) orthologs downloaded from: https://download.rgd.mcw.edu/data_release/
#

def get_GO_genes_API(term):
    """
    Retrieves all genes associated with each term.
    Input of GO terms must be a 1d list of GO Term Accession. e.g. ['GO:1903502','GO:1903508'].
    Homo sapiens taxon is NCBITaxon:9606
    """
    logger.info("get_GO_genes_API: term = " + term)
    parameters = {
        "rows": 10000000
    }
    
    response = requests.get(f"http://api.geneontology.org/api/bioentity/function/{term}/genes", params=parameters) # Get JSON response for current term, read 'objects' property (array of genes) into 'genes' array
    # logger.debug(json.dumps(response.json(), indent=4))
    genes = []
    associations = response.json()['associations']
    for item in associations:
        # only use directly associated genes and genes
        if item['subject']['id'] in genes: #TODO: explain code
            logger.debug(f"Gene {item['subject']['id']} already in the list. Skipping...")
        elif item['object']['id'] == term and item['subject']['taxon']['id'] == "NCBITaxon:9606":
            genes.append(item['subject']['id'])
        elif not FLAG_HOMOSAPIENS_ONLY:
            if item['object']['id'] == term and any((database[0] in item['subject']['id'] and any(taxon in item['subject']['taxon']['id'] for taxon in database[1])) for database in APPROVED_DATABASES):
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
    logger.debug(f"[get_ensembl_sequence_API] Starting Ensembl API for id {id}")
    response = requests.get(f"https://rest.ensembl.org/sequence/id/{id}?object_type=transcript;type=cds", headers={"Content-Type": "text/plain", }
                            )  # cds = c-DNA without the UTR regions; type=cdna (in this case UTR region is also kept); retrieves complementary sequence to the gene mRNA (without UTR), same as miRNA sequence (todo: preveri z genetikom)
    if response.ok:
        logger.debug(response.text)
        logger.info(f"[get_ensembl_sequence_API] Recieved sequence for id {id}.")
        return response.text
    else:
        logger.info(f"[get_ensembl_sequence_API] Failed to get sequence for id {id}")
        return None


def get_rnacentral_sequence_API(id):
    """
    This function queries RNA Central for nucleotide sequence
    Input of RNA Central ID's must be a string
    """
    logger.debug(f"[get_rnacentral_sequence_API] Starting RNACentral API for id {id}")
    response = requests.get(
        f"http://rnacentral.org/api/v1/rna/{id}/?format=json")
    if response.ok:
        logger.debug(response.json())
        sequence = response.json()['sequence']
        logger.info(f"[get_rnacentral_sequence_API] Recieved sequence for id {id} -> {sequence}.")
        return sequence
    else:
        logger.info(f"[get_rnacentral_sequence_API] RNACentral API error")
        return None


def uniprot_mapping_API(id_old, source='UniProtKB_AC-ID', target='Ensembl_Transcript'):  # !
    """
    Recieves uniprot or other ID and finds the Ensembl id.
    Input of ID's must be a string
    https://www.uniprot.org/help/id_mapping
    """
    logger.debug(f"[uniprot_mapping]: id_old = {id_old}")
    source = "UniProtKB_AC-ID"  # try to map all genes from other databases to UniProt

    # TODO: some terms (like GO-1903670) have genes that are not defined in UniProt! For example, one gene from
    # GO-1903670 has id_old ZFIN:ZDB-GENE-041014-357, throws UnboundLocalError (source referenced before assignment)
    # possible solution: source = ""
    # but this omits any databases that are not uniprot
    # TODO: some terms return multiple id, which to choose???

    id = id_old.split(':')[1]
    response = requests.post(f"https://rest.uniprot.org/idmapping/run",
                             data={'from': source, 'to': target, 'ids': id})
    logger.debug(response.json())
    jobId = response.json()['jobId']
    logger.debug(jobId)
    POLLING_INTERVAL = 1
    while True:
        request = requests.get(
            f"https://rest.uniprot.org/idmapping/status/{jobId}")
        j = request.json()
        logger.debug(j)
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                logger.debug(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(j["jobStatus"])
        else:
            if 'failedIds' in j:
                ensembl_id = None
            else:  # is index 0 really the best one?
                ensembl_id = j['results'][0]['to'].split('.')[0]#TODO:the first intex 0 automatically selects the first results (is that ok?)
            break
    logger.info(f"mapped {id_old} to {ensembl_id}")
    return ensembl_id


def find_genes_related_to_GO_terms(terms, ask_for_overrides=True, destination_folder="term_genes"):
    """
    Finds the genes related to the terms array and dumps the results into a json file.
    """
    logger.info(f"terms array len:{len(terms)}, elements: {terms}")
    for term in terms:
        term_file = str(term).replace(":", "-")
        filepath = f"{destination_folder}/{term_file}.json"
        _find_genes_related_to_GO_term(term, filepath, ask_for_overrides)


def _find_genes_related_to_GO_term(term, filepath, ask_for_overrides):
    """
    Finds the genes related to the term and dumps the results into a json file.
    """
    def _handle_load_from_crash(crash_filepath): # inner function as way of both encapsulation and code reuse
        nonlocal crash_json
        crash_json = util.read_file_as_json(crash_filepath)
        if json.dumps(crash_json) != "[]": # if crash_json isn't empty
            last_geneId = util.get_last_geneId_in_crash_json(crash_json)
            nonlocal genes
            genes = util.list_directionshrink(genes, last_geneId, forward=True) # shrink genes to start processing at the first gene after last_geneId
            logger.info(f"Crash recovery: last_geneId = {last_geneId}, genes_len = {len(genes)}")
    
    def _choose_crashfile(crash_filepaths):
        """
        Gives user the choise to choose a crash file among crash_filepaths
        """
        display_dictionary = {}
        i=0
        for path in crash_filepaths:
            display_dictionary[i] = path
            i += 1
        choice = int(input(f"Enter the number of the crashfile from {display_dictionary}"))
        return display_dictionary[choice]

    logger.debug(f"started gene search for GO term {term}")
    global current_filepath
    current_filepath = filepath

    override = 0
    if os.path.isfile(filepath) and ask_for_overrides == True:
        override = input(f"File {filepath} already exists. Enter 1 to process the file again or 0 to skip:")
        if int(override) == 0:  # careful! override changes type to str as input is given
            logger.info(f"Skipping file {filepath}")
            return
    
    genes = get_GO_genes_API(term)  # get array of genes associated to a term
    e_id = []  # ensemble id
    sequences = []
    global json_dictionaries

    # crash recovery code
    _f = filepath.split("/")
    _fn = "" # this is the filename eg. GO-0001252
    for element in _f: # finds the element with .json and assigns it to f
        if ".json" in element: _fn = element.replace(".json", "")
    
    crash_filepaths = util.get_files_in_dir("term_genes_crash", _fn)
    logger.debug(f"crash_filepaths = {crash_filepaths}")
    crash_json = ""
    crash_filepath = ""
    if (isinstance(crash_filepaths, list) and len(crash_filepaths)>=1) or crash_filepaths != "":
        crash_filepath = util.get_last_file_in_list(crash_filepaths) # get last of the crashes
    if os.path.exists(crash_filepath):
        if len(crash_filepaths)>1:
            restore_crash = int(input(f"File {crash_filepath} exists for term {term} as an option for crash recovery. Press 1 to recover data and delete the file, 2 to recover data and keep the file, 3 to display other crash files or 0 to ignore it."))
        else: 
            restore_crash = int(input(f"File {crash_filepath} exists for term {term} as an option for crash recovery. Press 1 to recover data and delete the file, 2 to recover data and keep the file or 0 to ignore it."))
        if restore_crash == 3:
            crash_filepath = _choose_crashfile(crash_filepaths)
            _handle_load_from_crash(crash_filepath)
        elif restore_crash == 1: # load crash_filepath json and delete file
            _handle_load_from_crash(crash_filepath)
            os.remove(crash_filepath)
        elif restore_crash == 2: # load crash_filepath json and keep file
            _handle_load_from_crash(crash_filepath)
        elif restore_crash == 0: # do nthn
            logger.info(f"Crash recovery not selected.")
        if crash_json != "":
            logger.info(f"Crash recovery: json appended.")
            json_dictionaries.append(crash_json)
    else:
        logger.debug(f"Crash filepath {crash_filepath} doesn't exist. Recovery not started.")

    # TODO: Code reuse using pass-by-reference, which is handy in Python 3.0 with the 'nonlocal' keyword
    # -> problem is e_id and sequences, which stay in the local scope aka cannot modify their value in the current
    # scope from another function
    # https://stackoverflow.com/questions/8447947/is-it-possible-to-modify-a-variable-in-python-that-is-in-an-outer-enclosing-b

    len_genes = len(genes)
    i = 0
    for gene in genes:  # gene is is shape prefix:id
        i += 1
        logging.info(f"[_find_genes_related_to_GO_term]: Processing gene {gene}; {i}/{len_genes}")
        if 'UniProtKB' in gene:
            e_id.append(uniprot_mapping_API(gene))
            sequences.append(get_ensembl_sequence_API(e_id[-1]))
        elif 'ZFIN' in gene:
            human_gene_symbol = util.zfin_find_human_ortholog(gene)  # eg. adgrg9
            # TODO: compute ensembl sequence!
            if "ZfinError" in human_gene_symbol:
                logger.debug(
                    f"[uniprot_mapping]: ERROR! human_gene_symbol for {gene} was not found!")
                # input("Press enter to proceed.") # no need
                e_id.append(None)
                sequences.append(None)
            else:  # human ortholog exists in uniprot
                e_id.append(uniprot_mapping_API(util.get_uniprotId_from_geneName_new(
                    human_gene_symbol, trust_genes=FLAG_TRUST_GENES)))
                logger.debug(f"id_old = {e_id[-1]}")
                if "CycleOutOfBoundsError" in e_id[-1] or e_id[-1] == 0:
                    e_id[-1] = None
                    sequences.append(None)
                else:
                    sequences.append(get_ensembl_sequence_API(e_id[-1]))
        elif 'RNAcentral' in gene:
            e_id.append(gene.split(':')[1])
            sequences.append(get_rnacentral_sequence_API(e_id[-1]))
        elif "Xenbase" in gene:
            human_gene_symbol = util.xenbase_find_human_ortholog(gene)
            if "XenbaseError" in human_gene_symbol:
                logger.info(human_gene_symbol) # error msg is human_gene_symbol
                e_id.append(None)
                sequences.append(None)
            else: # human ortholog exists in xenbase
                e_id.append(uniprot_mapping_API(util.get_uniprotId_from_geneName_new(human_gene_symbol, trust_genes=FLAG_TRUST_GENES))) # TODO: this is code repetition -> create a function and try to pass e_id by reference!!! (or there will be errors)
                logger.debug(f"id_old = {e_id[-1]}")
                if "CycleOutOfBoundsError" in e_id[-1] or e_id[-1] == 0:
                    e_id[-1] = None
                    sequences.append(None)
                else:
                    sequences.append(get_ensembl_sequence_API(e_id[-1]))
        elif "MGI" in gene:
            human_gene_symbol = util.mgi_find_human_ortholog(gene)
            if "MgiError" in human_gene_symbol:
                logger.info(human_gene_symbol) # error msg is human_gene symbol
                e_id.append(None)
                sequences.append(None)
            else: # human ortholog exists in mgi
                e_id.append(uniprot_mapping_API(util.get_uniprotId_from_geneName_new(human_gene_symbol, trust_genes=FLAG_TRUST_GENES)))
                logger.debug(f"id_old = {e_id[-1]}")
                if "CycleOutOfBoundsError" in e_id[-1] or e_id[-1] == 0:
                    e_id[-1] = None
                    sequences.append(None)
                else:
                    sequences.append(get_ensembl_sequence_API(e_id[-1]))
        elif "RGD" in gene:
            human_gene_symbol = util.rgd_find_human_ortholog(gene)
            if "RgdError" in human_gene_symbol:
                logger.info(human_gene_symbol)
                e_id.append(None)
                sequences.append(None)
            else: # human ortholog exists in rgd
                e_id.append(uniprot_mapping_API(util.get_uniprotId_from_geneName_new(human_gene_symbol, trust_genes=FLAG_TRUST_GENES)))
                logger.debug(f"id_old = {e_id[-1]}")
                if e_id[-1] == None: # request for UniProtKB:Q5YKI7 maps to None --> raises 'TypeError: argument of type NoneType is not iterable' in the elif clause (when "CycleOutOfBoundsError" is checked for in None) -> this first if check solves this
                    e_id[-1] = None
                    sequences.append(None)
                elif "CycleOutOfBoundsError" in e_id[-1] or e_id[-1] == 0:
                    e_id[-1] = None
                    sequences.append(None)
                else:
                    sequences.append(get_ensembl_sequence_API(e_id[-1]))
        else:
            input(f"No database found for {gene}. Press any key to continue.")
            e_id.append(None)
            sequences.append(None)
        out = {"term": term, "product": gene, "sequence_id": e_id[-1], "sequence": sequences[-1]}
        # f.write(json.dumps(out)+"\n") For file decoding purposes, json dicts need to be stored in a list and then the list written to the file as per https://stackoverflow.com/questions/21058935/python-json-loads-shows-valueerror-extra-data
        json_dictionaries.append(out)

    #file = open(filepath, "w+")
    #file.write(json.dumps(json_dictionaries)+"\n")
    #file.close()
    util.store_json_dictionaries(filepath, json_dictionaries)
    logger.debug(f"finished gene search for GO term {term}")

def exit_handler():
    """
    Executes last code before program exit. If any file is being processed, it's flagged.
    If there is any last IO operations etc, perform them here.
    """

    # store json dictionaries on crash
    # TODO: make snapshots (so stopping console at ctrl-c doesn't yeet all the results from previous file, offer user latest snapshot default, but allow snapshot selection)
    filename = current_filepath.split("/")[len(current_filepath.split("/"))-1].replace(".json", "") # gets the last item in path eg. GO-0001525.json
    dest = f"term_genes_crash/{filename}_{datetime.datetime.now().timestamp()}_.json"
    util.store_json_dictionaries(dest, json_dictionaries)
    logging.info("Stopping script!")


def main():
    # register listeners/handlers
    atexit.register(exit_handler)

    # set global variables -> initialise with global keyword to refer to the global variable, then set value!
    global FLAG_HOMOSAPIENS_ONLY
    FLAG_HOMOSAPIENS_ONLY = False

    global FLAG_TRUST_GENES
    FLAG_TRUST_GENES = True

    # Compare any json files for testing -> get_GO_genes_API_new with homosapiens_only=True works as expected
    # logging.info(util.json_compare("term_genes/GO-0001525.json", "term_genes/homosapiens_only,v1/GO-0001525.json"))
    # logging.info(util.json_compare("term_genes/GO-0045765.json", "term_genes/homosapiens_only,v1/GO-0045765.json"))

    # dev_test_api_download.get_GO_genes_API("GO:1903670")
    # terms_test = ['GO:1903587']
    # terms_angiogenesis_ids = util.get_array_terms("ANGIOGENESIS")

    # startup functions
    util.load_trusted_genes("src_data_files/genes_trusted.txt")
    logging.info(f"Loaded {len(constants.TRUSTED_GENES)} trusted genes.")
    util.load_human_orthologs()

    # main functions
    terms_all = util.get_array_terms("ALL")
    find_genes_related_to_GO_terms(terms_all, destination_folder="term_genes/homosapiens_only=false,v1")
    

    # showcase functions:
    # this is how to retrieve uniprotId description (function) from uniprotId:
    # logging.info(util.get_uniprotId_description("O14944"))
    # logging.info(util.get_uniprotId_description("Q9BUL8"))

    # debugging functions: RGD:1359373
    # human_gene_symbol = util.rgd_find_human_ortholog("RGD:1359373")
    # util.get_uniprotId_from_geneName_new(human_gene_symbol, trust_genes=FLAG_TRUST_GENES)

    # cleanup files, turn on only once
    # util.cleanup_crash_term_jsons()
    


if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    main()

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
