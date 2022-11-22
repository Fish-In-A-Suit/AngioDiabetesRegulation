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
# trust genes in trusted_genes to be credible -> program doesnt stop at them for validation
FLAG_TRUST_GENES = True

APPROVED_DATABASES = [["UniProtKB", ["NCBITaxon:9606"]],
                      ["ZFIN", ["NCBITaxon:7955"]],
                      ["RNAcentral", ["NCBITaxon:9606"]],
                      ["Xenbase", ["NCBITaxon:8364"]],
                      ["MGI", ["NCBITaxon:10090"]],
                      ["RGD", ["NCBITaxon:10116"]]]  # Which pairs of databases and taxons (could be multiple per database) to allow.

json_dictionaries = []
current_filepath = ""

# Remarks
#   - no need to track _current_file, since .json files in term_genes aren't saved if script terminates early / if all the json elements haven't already been computed
#   - RGD (Rat Genome Database) orthologs downloaded from: https://download.rgd.mcw.edu/data_release/
#


def get_GO_products_from_term_API(term):
    """
    Retrieves all genes associated with each term. Input of GO terms must be a 1d list of GO Term Accession. e.g. ['GO:1903502','GO:1903508'].
    Homo sapiens taxon is NCBITaxon:9606
    """
    logger.info("get_GO_genes_API: term = " + term)
    parameters = {
        "rows": 10000000
    }

    # Get JSON response for current term
    response = requests.get(
        f"http://api.geneontology.org/api/bioentity/function/{term}/genes", params=parameters)
    genes = []
    associations = response.json()['associations']
    for item in associations:
        # only use directly associated genes and genes
        # see documentation docx document (involves traversing the json)
        if item['subject']['id'] in genes:
            logger.debug(
                f"Gene {item['subject']['id']} already in the list. Skipping...")
        elif item['object']['id'] == term and item['subject']['taxon']['id'] == "NCBITaxon:9606":
            genes.append(item['subject']['id'])
        elif not FLAG_HOMOSAPIENS_ONLY:
            if item['object']['id'] == term and any((database[0] in item['subject']['id'] and any(taxon in item['subject']['taxon']['id'] for taxon in database[1])) for database in APPROVED_DATABASES):
                genes.append(item['subject']['id'])
    # IMPORTANT: Some terms (like GO:1903587) return only genes related to "subterms" (when calling http://api.geneontology.org:80 "GET /api/bioentity/function/GO%3A1903587/genes?use_compact_associations=True&taxon=NCBITaxon%3A9606 HTTP/1.1" 200 1910)
    # --> no genes associated to the term, only to subterms --> genes array can be of 0 length (and that is not an error)
    if len(genes) == 0:
        # 0 associated gene products -> save in appropriate terms_empty.txt file
        util.append_to_file(term, os.path.join(
            util.filepath_striplast(current_filepath), "terms_empty.txt"))
        logger.info(
            f"Term {term} doesn't contain any products. Saved it in {current_filepath}")
        return []

    logger.info(
        f"Term {term} has {len(genes)} associated products -> {genes}.")
    return genes


def get_ensembl_sequence_API(id):
    """
    This function queries ensembl for nucleotide sequence
    Input of ensembl ID's must be a string
    """
    logger.debug(
        f"[get_ensembl_sequence_API] Starting Ensembl API for id {id}")
    response = requests.get(f"https://rest.ensembl.org/sequence/id/{id}?object_type=transcript;type=cds", headers={"Content-Type": "text/plain", }
                            )  # cds = c-DNA without the UTR regions; type=cdna (in this case UTR region is also kept); retrieves complementary sequence to the gene mRNA (without UTR), same as miRNA sequence (todo: preveri z genetikom)
    if response.ok:
        logger.debug(response.text)
        logger.info(
            f"[get_ensembl_sequence_API] Recieved sequence for id {id}.")
        return response.text
    else:
        logger.info(
            f"[get_ensembl_sequence_API] Failed to get sequence for id {id}")
        return None


def get_rnacentral_sequence_API(id):
    """
    This function queries RNA Central for nucleotide sequence
    Input of RNA Central ID's must be a string
    """
    logger.debug(
        f"Starting RNACentral API for id {id}")
    response = requests.get(
        f"http://rnacentral.org/api/v1/rna/{id}/?format=json")
    if response.ok:
        logger.debug(response.json())
        sequence = response.json()['sequence']
        logger.info(
            f"Recieved sequence for id {id} -> {sequence}.")
        return sequence
    else:
        logger.info(f"RNACentral API error")
        return None

def find_products_related_to_GO_terms_new(terms, folder):
    """
    
    """
    def _handle_load_from_crash(crash_filepath):  # inner function as way of both encapsulation and code reuse
        nonlocal crash_json
        crash_json = util.read_file_as_json(crash_filepath)
        if json.dumps(crash_json) != "[]":  # if crash_json isn't empty
            last_GO_term = util.get_last_GO_term_in_crash_json(crash_json)
            nonlocal terms
            # shrink genes to start processing at the first gene after last_geneId
            terms = util.list_directionshrink(terms, last_GO_term, forward=True)
            logger.info(
                f"Crash recovery: last_geneId = {last_GO_term}, genes_len = {len(terms)}")

    def _choose_crashfile(crash_filepaths):
        """
        Gives user the choise to choose a crash file among crash_filepaths
        """
        display_dictionary = {}
        i = 0
        for path in crash_filepaths:
            display_dictionary[i] = path
            i += 1
        choice = int(
            input(f"Enter the number of the crashfile from {display_dictionary}"))
        return display_dictionary[choice]

    logger.info(f"terms array len:{len(terms)}, elements: {terms}")
    global json_dictionaries
    json_dictionaries=[]
    global current_filepath
    current_filepath = os.path.join(folder, "term_products.json")

    if os.path.isfile(current_filepath):
        override = input(
            f"File {current_filepath} already exists. Enter 1 to process the file again or 0 to skip:")
        if int(override) == 0:  # careful! override changes type to str as input is given
            logger.info(f"Skipping file {current_filepath}")
            return

    # crash recovery code
    _f = current_filepath.split("/")
    _fn = ""  # this is the filename eg. GO-0001252
    for element in _f:  # finds the element with .json and assigns it to f
        if ".json" in element:
            _fn = element.replace(".json", "")

    crash_filepaths = util.get_files_in_dir("term_genes_crash", _fn)
    logger.debug(f"crash_filepaths = {crash_filepaths}")
    crash_json = ""
    crash_filepath = ""
    if (isinstance(crash_filepaths, list) and len(crash_filepaths) >= 1) or crash_filepaths != "":
        crash_filepath = util.get_last_file_in_list(
            crash_filepaths)  # get last of the crashes
    if os.path.exists(crash_filepath):
        if len(crash_filepaths) > 1:
            restore_crash = int(input(
                f"File {crash_filepath} exists as an option for crash recovery. Press 1 to recover data and delete the file, 2 to recover data and keep the file, 3 to display other crash files or 0 to ignore it."))
        else:
            restore_crash = int(input(
                f"File {crash_filepath} exists as an option for crash recovery. Press 1 to recover data and delete the file, 2 to recover data and keep the file or 0 to ignore it."))
        if restore_crash == 3:
            crash_filepath = _choose_crashfile(crash_filepaths)
            _handle_load_from_crash(crash_filepath)
        elif restore_crash == 1:  # load crash_filepath json and delete file
            _handle_load_from_crash(crash_filepath)
            os.remove(crash_filepath)
        elif restore_crash == 2:  # load crash_filepath json and keep file
            _handle_load_from_crash(crash_filepath)
        elif restore_crash == 0:  # do nthn
            logger.info(f"Crash recovery not selected.")
        if crash_json != "":
            logger.info(f"Crash recovery: json appended.")
            # adding each element (instead of entire json object) to prevent the repeated nesting bug
            for i in range(len(crash_json)):
                json_dictionaries.append(crash_json[i])
    else:
        logger.info(
            f"Crash filepath {crash_filepath} doesn't exist. Recovery not started.")

    # TODO: Code reuse using pass-by-reference, which is handy in Python 3.0 with the 'nonlocal' keyword
    # -> problem is e_id and sequences, which stay in the local scope aka cannot modify their value in the current
    # scope from another function
    # https://stackoverflow.com/questions/8447947/is-it-possible-to-modify-a-variable-in-python-that-is-in-an-outer-enclosing-b


    #TODO:crash recovery
    
    for term in terms:
        term_products=_find_products_related_to_GO_term_new(term)
        json_dictionaries.append({"GO_term":term, "products":term_products})
    
    util.store_json_dictionaries(current_filepath, json_dictionaries)


def _find_products_related_to_GO_term_new(term):
    """
    
    """
    logger.info(f"started product search for GO term {term}")
    direct_products = get_GO_products_from_term_API(term)

    i = 0
    uniprot_productnames = []
    for product in direct_products:
        i+=1
        logging.debug(
            f"Processing product {product}; {i}/{len(direct_products)}")
        
        def ortholog_process(gene,product):
            if "Error" in gene:
                logger.debug(
                    f"ERROR! human_gene_symbol for {product} was not found!")
                uniprot_productnames.append(None)
            else:  # human ortholog exists in uniprot
                uniprot_productnames.append(util.get_uniprotId_from_geneName_new(
                    gene, trust_genes=FLAG_TRUST_GENES))
        
        if 'UniProtKB' in product:
            uniprot_productnames.append(product)
        elif 'ZFIN' in product:
            human_gene_symbol = util.zfin_find_human_ortholog(product)  # eg. adgrg9
            ortholog_process(human_gene_symbol,product)
        elif 'RNAcentral' in product: #TODO:what to do with miRNA sequences?, perhaps export them separately
            uniprot_productnames.append(None)
        elif "Xenbase" in product:
            human_gene_symbol = util.xenbase_find_human_ortholog(product)
            ortholog_process(human_gene_symbol,product)
        elif "MGI" in product:
            human_gene_symbol = util.mgi_find_human_ortholog(product)
            ortholog_process(human_gene_symbol,product)
        elif "RGD" in product:
            human_gene_symbol = util.rgd_find_human_ortholog(product)
            ortholog_process(human_gene_symbol,product)
        else:
            # input(f"No database found for {gene}. Press any key to continue.")
            logger.debug(f"No database found for {product}")
            unsorted_genes_filepath = os.path.join("genes_unsorted.txt") #TODO:filepath
            unsorted_genes_filepath = os.path.join(
                util.filepath_striplast(current_filepath), "genes_unsorted.txt")
            util.append_to_file(f"{term} {product}", unsorted_genes_filepath)
            uniprot_productnames.append(None)

    logger.debug(f"finished gene search for GO term {term}")
    
    return [{"direct_productID":direct_products[i], "UniprotID":uniprot_productnames[i]} 
            for i in range(len(direct_products))]




def find_genes_related_to_GO_terms(terms, ask_for_overrides=True, destination_folder="term_genes"):
    """
    Finds the genes related to the terms array and dumps the results into a json file.
    """
    logger.info(f"terms array len:{len(terms)}, elements: {terms}")
    for term in terms:
        term_file = str(term).replace(":", "-")
        filepath = f"{destination_folder}/{term_file}.json"
        global current_filepath
        current_filepath = filepath
        _find_genes_related_to_GO_term(term, filepath, ask_for_overrides)
        global json_dictionaries
        json_dictionaries = []  # reset


def _find_genes_related_to_GO_term(term, filepath, ask_for_overrides):
    """
    Finds the genes related to the term and dumps the results into a json file.
    """
    def _handle_load_from_crash(crash_filepath):  # inner function as way of both encapsulation and code reuse
        nonlocal crash_json
        crash_json = util.read_file_as_json(crash_filepath)
        if json.dumps(crash_json) != "[]":  # if crash_json isn't empty
            last_geneId = util.get_last_geneId_in_crash_json(crash_json)
            nonlocal genes
            # shrink genes to start processing at the first gene after last_geneId
            genes = util.list_directionshrink(genes, last_geneId, forward=True)
            logger.info(
                f"Crash recovery: last_geneId = {last_geneId}, genes_len = {len(genes)}")

    def _choose_crashfile(crash_filepaths):
        """
        Gives user the choise to choose a crash file among crash_filepaths
        """
        display_dictionary = {}
        i = 0
        for path in crash_filepaths:
            display_dictionary[i] = path
            i += 1
        choice = int(
            input(f"Enter the number of the crashfile from {display_dictionary}"))
        return display_dictionary[choice]

    logger.debug(f"started gene search for GO term {term}")

    override = 0
    if os.path.isfile(filepath) and ask_for_overrides == True:
        override = input(
            f"File {filepath} already exists. Enter 1 to process the file again or 0 to skip:")
        if int(override) == 0:  # careful! override changes type to str as input is given
            logger.info(f"Skipping file {filepath}")
            return
    elif os.path.isfile(filepath) and ask_for_overrides == False:
        # file exists, skip
        logger.info(f"Skipping file {filepath}")
        return

    # get array of genes associated to a term
    genes = get_GO_products_from_term_API(term)
    e_id = []  # ensemble id
    sequences = []
    global json_dictionaries

    # crash recovery code
    _f = filepath.split("/")
    _fn = ""  # this is the filename eg. GO-0001252
    for element in _f:  # finds the element with .json and assigns it to f
        if ".json" in element:
            _fn = element.replace(".json", "")

    crash_filepaths = util.get_files_in_dir("term_genes_crash", _fn)
    logger.debug(f"crash_filepaths = {crash_filepaths}")
    crash_json = ""
    crash_filepath = ""
    if (isinstance(crash_filepaths, list) and len(crash_filepaths) >= 1) or crash_filepaths != "":
        crash_filepath = util.get_last_file_in_list(
            crash_filepaths)  # get last of the crashes
    if os.path.exists(crash_filepath):
        if len(crash_filepaths) > 1:
            restore_crash = int(input(
                f"File {crash_filepath} exists for term {term} as an option for crash recovery. Press 1 to recover data and delete the file, 2 to recover data and keep the file, 3 to display other crash files or 0 to ignore it."))
        else:
            restore_crash = int(input(
                f"File {crash_filepath} exists for term {term} as an option for crash recovery. Press 1 to recover data and delete the file, 2 to recover data and keep the file or 0 to ignore it."))
        if restore_crash == 3:
            crash_filepath = _choose_crashfile(crash_filepaths)
            _handle_load_from_crash(crash_filepath)
        elif restore_crash == 1:  # load crash_filepath json and delete file
            _handle_load_from_crash(crash_filepath)
            os.remove(crash_filepath)
        elif restore_crash == 2:  # load crash_filepath json and keep file
            _handle_load_from_crash(crash_filepath)
        elif restore_crash == 0:  # do nthn
            logger.info(f"Crash recovery not selected.")
        if crash_json != "":
            logger.info(f"Crash recovery: json appended.")
            # adding each element (instead of entire json object) to prevent the repeated nesting bug
            for i in range(len(crash_json)):
                json_dictionaries.append(crash_json[i])
    else:
        logger.info(
            f"Crash filepath {crash_filepath} doesn't exist. Recovery not started.")

    # TODO: Code reuse using pass-by-reference, which is handy in Python 3.0 with the 'nonlocal' keyword
    # -> problem is e_id and sequences, which stay in the local scope aka cannot modify their value in the current
    # scope from another function
    # https://stackoverflow.com/questions/8447947/is-it-possible-to-modify-a-variable-in-python-that-is-in-an-outer-enclosing-b

    len_genes = len(genes)
    i = 0
    for gene in genes:  # gene is is shape prefix:id
        i += 1
        logging.info(
            f"Processing gene {gene}; {i}/{len_genes}")
        if 'UniProtKB' in gene:
            e_id.append(util.get_uniprotId_from_geneName_new(gene)[1])
            sequences.append(get_ensembl_sequence_API(e_id[-1]))
        elif 'ZFIN' in gene:
            human_gene_symbol = util.zfin_find_human_ortholog(
                gene)  # eg. adgrg9
            # TODO: compute ensembl sequence!
            if "ZfinError" in human_gene_symbol:
                logger.debug(
                    f"[uniprot_mapping]: ERROR! human_gene_symbol for {gene} was not found!")
                # input("Press enter to proceed.") # no need
                e_id.append(None)
                sequences.append(None)
            else:  # human ortholog exists in uniprot
                e_id.append(util.get_uniprotId_from_geneName_new(
                    human_gene_symbol, trust_genes=FLAG_TRUST_GENES)[1])
                logger.debug(f"id_old = {e_id[-1]}")
                if e_id[-1] == None:  # prevents 'TypeError: argument of type NoneType is not iterable'
                    e_id[-1] = None
                    sequences.append(None)
                elif "CycleOutOfBoundsError" in e_id[-1] or e_id[-1] == 0:
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
                # error msg is human_gene_symbol
                logger.info(human_gene_symbol)
                e_id.append(None)
                sequences.append(None)
            else:  # human ortholog exists in xenbase
                # TODO: this is code repetition -> create a function and try to pass e_id by reference!!! (or there will be errors)
                e_id.append(util.get_uniprotId_from_geneName_new(
                    human_gene_symbol, trust_genes=FLAG_TRUST_GENES)[1])
                logger.debug(f"id_old = {e_id[-1]}")
                if e_id[-1] == None:  # prevents 'TypeError: argument of type NoneType is not iterable'
                    e_id[-1] = None
                    sequences.append(None)
                elif "CycleOutOfBoundsError" in e_id[-1] or e_id[-1] == 0:
                    e_id[-1] = None
                    sequences.append(None)
                else:
                    sequences.append(get_ensembl_sequence_API(e_id[-1]))
        elif "MGI" in gene:
            human_gene_symbol = util.mgi_find_human_ortholog(gene)
            if "MgiError" in human_gene_symbol:
                # error msg is human_gene symbol
                logger.info(human_gene_symbol)
                e_id.append(None)
                sequences.append(None)
            else:  # human ortholog exists in mgi
                e_id.append(util.get_uniprotId_from_geneName_new(
                    human_gene_symbol, trust_genes=FLAG_TRUST_GENES)[1])
                logger.debug(f"id_old = {e_id[-1]}")
                if e_id[-1] == None:  # prevents 'TypeError: argument of type NoneType is not iterable'
                    e_id[-1] = None
                    sequences.append(None)
                elif "CycleOutOfBoundsError" in e_id[-1] or e_id[-1] == 0:
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
            else:  # human ortholog exists in rgd
                e_id.append(util.get_uniprotId_from_geneName_new(
                    human_gene_symbol, trust_genes=FLAG_TRUST_GENES)[1])
                logger.debug(f"id_old = {e_id[-1]}")
                # request for UniProtKB:Q5YKI7 maps to None --> raises 'TypeError: argument of type NoneType is not iterable' in the elif clause (when "CycleOutOfBoundsError" is checked for in None) -> this first if check solves this
                if e_id[-1] == None:
                    e_id[-1] = None
                    sequences.append(None)
                elif "CycleOutOfBoundsError" in e_id[-1] or e_id[-1] == 0:
                    e_id[-1] = None
                    sequences.append(None)
                else:
                    sequences.append(get_ensembl_sequence_API(e_id[-1]))
        else:
            # input(f"No database found for {gene}. Press any key to continue.")
            logger.debug(f"No database found for {gene}")
            unsorted_genes_filepath = os.path.join(
                util.filepath_striplast(filepath), "genes_unsorted.txt")
            util.append_to_file(f"{term} {gene}", unsorted_genes_filepath)
            e_id.append(None)
            sequences.append(None)
        out = {"term": term, "product": gene,
               "sequence_id": e_id[-1], "sequence": sequences[-1]}
        # f.write(json.dumps(out)+"\n") For file decoding purposes, json dicts need to be stored in a list and then the list written to the file as per https://stackoverflow.com/questions/21058935/python-json-loads-shows-valueerror-extra-data
        json_dictionaries.append(out)

    #file = open(filepath, "w+")
    # file.write(json.dumps(json_dictionaries)+"\n")
    # file.close()
    util.store_json_dictionaries(filepath, json_dictionaries)
    logger.debug(f"finished gene search for GO term {term}")


def exit_handler():
    """
    Executes last code before program exit. If any file is being processed, it's flagged.
    If there is any last IO operations etc, perform them here.
    """
    filename = current_filepath.split("/")[len(current_filepath.split("/"))-1].replace(
        ".json", "")  # gets the last item in path eg. GO-0001525.json
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
    # util.load_trusted_genes("src_data_files/genes_trusted.txt")
    util.load_list_from_file("src_data_files/genes_trusted.txt",
                             constants.TRUSTED_GENES, no_elements_in_line=2, break_character=" ")
    util.load_list_from_file(
        "term_genes/homosapiens_only=false,v1/terms_empty.txt", constants.TERMS_EMPTY)
    logging.info(
        f"Loaded {len(constants.TRUSTED_GENES)} trusted genes. trusted_genes = {constants.TRUSTED_GENES}")
    logging.info(
        f"Loaded {len(constants.TERMS_EMPTY)} empty terms. terms_empty = {constants.TERMS_EMPTY}")
    util.load_human_orthologs()

    # main functions
    # terms_all = util.get_array_terms("ALL")
    #terms = ["GO:1904204"]
    # find_genes_related_to_GO_terms(terms, ask_for_overrides = False, destination_folder="term_genes/homosapiens_only=false,v1")

    # new main functions
    #terms = ["GO:0016525"]
    terms = util.get_array_terms("ALL")
    find_products_related_to_GO_terms_new(terms, filepath="term_genes/homosapiens_only=false,v2/term_products.json")

    #response = requests.get(
    #    f"http://api.geneontology.org/api/bioentity/function/GO:0001525/genes", params={"rows": 1000000})
    #util.save_json(response.json(), "test/test.json")
    # logger.info(response.json())

    # showcase functions:
    # this is how to retrieve uniprotId description (function) from uniprotId:
    # logging.info(util.get_uniprotId_description("O14944"))
    # logging.info(util.get_uniprotId_description("Q9BUL8"))

    # this is how to sort json files into new folders by their respective array
    # util.term_sort_into_file("term_genes/homosapiens_only=false,v1", "test", constants.TERMS_DIABETES_NEGATIVE_ARRAY)

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


"""
set genov za diabetes - set genov za angiogenezo
presek
najboljših 10 al 20
poišči mRNA za te gene -> napovej miRNA
pathwayi ?
"""
