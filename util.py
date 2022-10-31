# Utility functions file
import requests
import constants
import json

import logging
import sys

_response_cycle_counter = 0
_uniprot_identifier_query_result = ""


console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
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

def get_array_terms(array_name, term_shrink=True):
    """
    Returns the specified array terms, possible options:
      - 'ALL': all arrays combined
      - 'ANGIOGENESIS': all angiogenesis arrays combined
      - 'ANGIOGENESIS-POSITIVE'
      - 'ANGIOGENESIS-NEGATIVE'
      - 'ANGIOGENESIS-GENERAL'
      - 'DIABETES': all diabetes arrays combined
      - 'DIABETES-POSITIVE'
      - 'DIABETES-NEGATIVE'
      - 'DIABETES-GENERAL'
    
    If term_shrink=True, the final list will only consist of GO:xxxxx, if false, returned list will consist of term names AND GO:xxxxx
    """
    if array_name == 'ALL':
        if term_shrink: return shrink_term_list(constants.TERMS_ANGIOGENESIS_GENERAL + constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY + constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY + constants.TERMS_DIABETES_GENERAL + constants.TERMS_DIABETES_NEGATIVE_ARRAY + constants.TERMS_DIABETES_POSITIVE_ARRAY)
        else: return constants.TERMS_ANGIOGENESIS_GENERAL + constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY + constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY + constants.TERMS_DIABETES_GENERAL + constants.TERMS_DIABETES_NEGATIVE_ARRAY + constants.TERMS_DIABETES_POSITIVE_ARRAY
    elif array_name == 'ANGIOGENESIS':
        if term_shrink: return shrink_term_list(constants.TERMS_ANGIOGENESIS_GENERAL + constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY + constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY)
        else: return constants.TERMS_ANGIOGENESIS_GENERAL + constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY + constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY
    elif array_name == 'ANGIOGENESIS-NEGATIVE':
        if term_shrink: return shrink_term_list(constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY)
        else: return constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY
    elif array_name == 'ANGIOGENESIS-POSITIVE':
        if term_shrink: return shrink_term_list(constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY)
        else: return constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY
    elif array_name == 'ANGIOGENESIS-GENERAL':
        if term_shrink: return shrink_term_list(constants.TERMS_ANGIOGENESIS_GENERAL)
        else: return constants.TERMS_ANGIOGENESIS_GENERAL
    elif array_name == 'DIABETES':
        if term_shrink: return shrink_term_list(constants.TERMS_DIABETES_GENERAL + constants.TERMS_DIABETES_NEGATIVE_ARRAY + constants.TERMS_DIABETES_POSITIVE_ARRAY)
        else: return constants.TERMS_DIABETES_GENERAL + constants.TERMS_DIABETES_NEGATIVE_ARRAY + constants.TERMS_DIABETES_POSITIVE_ARRAY
    elif array_name == 'DIABETES-NEGATIVE':
        if term_shrink: return shrink_term_list(constants.TERMS_DIABETES_NEGATIVE_ARRAY)
        else: return constants.TERMS_DIABETES_NEGATIVE_ARRAY
    elif array_name == 'DIABETES-POSITIVE':
        if term_shrink: return shrink_term_list(constants.TERMS_DIABETES_POSITIVE_ARRAY)
        else: return constants.TERMS_DIABETES_POSITIVE_ARRAY
    elif array_name == 'DIABETES-GENERAL':
        if term_shrink: return shrink_term_list(constants.TERMS_DIABETES_GENERAL)
        else: return constants.TERMS_DIABETES_GENERAL
    else:
        print(array_name + " could not be found! Returning empty array.")
        empty = []
        return empty

def read_file_as_json(filepath):
    """
    Reads the file into json
    """
    with open(filepath, "r") as read_content:
        return json.load(read_content)

def shrink_term_list(list):
    i = 0
    result_list = []
    for element in list:
        if i%2: # 0 = generic term name, 1 = GO:xxxx, 2 = generic term name, 3 = GO:xxxx -> modulo op is 1 at 1, 3, 5 etc
            result_list.append(element)
        i = i+1
    return result_list

def zfin_find_human_ortholog(gene_id, ortholog_file_path="src_data_files/zfin_human-zebrafish-gene-orthologs.txt"):
    """
    If gene_id is from the ZFIN database, searches through the zebrafish-human orthologs and returns the name of the
    symbol of the human gene ortholog.
    """
    logging.debug("[zfin_find_human_ortholog]: starting")
    #lines_firstpass = [1, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000, 21000, 22000, 23000, 24000, 25000, 26000, 27000, 28000, 29000, 30000, 31000, 32000, 33000, 34000, 35000, 36000, 37000, 38000, 39000, 40000, 41000, 42000, 43000]
    #lines_stepdown = []
    #gene_id_number = _zfin_split_row(gene_id)
    file = open(ortholog_file_path, "r")
    lines = file.readlines()
    gene_id=gene_id.split(":")[1] # eliminates ZFIN: 
    for line in lines:
        if gene_id in line:
            human_symbol = _zfin_get_human_gene_symbol_from_line(line)
            logging.debug(f"[zfin_find_human_ortholog]: Returning human symbol {human_symbol}")
            return human_symbol
    return f"[ZfinError_No-human-ortholog-found:gene_id={gene_id}"

    """ # Optimisation possibility:
    with open(ortholog_file_path, "r") as f:
        # split gene-id to get number
        for i, line in enumerate(f):
            if i in lines_firstpass or i in lines_stepdown:
                current_gene_number = _zfin_split_row(line)
                if current_gene_number < gene_id_number:
                    # still too far up, continue
                    continue
                elif current_gene_number > gene_id_number:
                    # too far down; reverse function to read backwards
                    for i in reversed(range(1000)):
                        # reverse until gene_id_number == current_gene_number
                        # PROBLEM: how to get file pointer again to read the line !?!??!
                        return 0 #delete this
    """

def _zfin_split_row(row):
    """
    Splits zfin row and retrieves a number, which is used for zfin_find_human_ortholog faster searching.
    Example: ZDB-GENE-041014-357	adgrg6	adhesion G protein-coupled receptor G6	ADGRG6	adhesion G protein-coupled receptor G6	612243	57211	13841	AA	ZDB-PUB-030905-1
    --> return: 041014
    """
    split = row.split("-")
    return int(split[2]) # the int(num) method also breaks leading zeroes (which is required)!

def _zfin_get_human_gene_symbol_from_line(line):
    """
    Splits zfin line and retrieves human gene symbol (full caps of zebrafish gene symbol)
    """
    return str(line.split("\t")[1]).upper() # split lines at tabs (spacebar is not ok!)

def get_uniprot_identifier(gene_name, prefix="UniProtKB:"):
    """
    Retrieves uniprot identifier e.g. UniProtKB:Q86SQ4, if gene_name=adgrg6
    """
    _response_cycle_counter = 0
    response = requests.get(f"https://rest.uniprot.org/uniprotkb/search?query={gene_name}+AND+organism_id:9606&format=json&fields=id,gene_names,organism_name")
    # logging.debug(f"Response = {response.text}")
    response = response.json()
    uniprot_gene_identifier = response["results"][0]["primaryAccession"]
    results_arr_len = len(response["results"])
    logging.debug(f"Gene name {gene_name} found to correspond to {uniprot_gene_identifier}. Displaying primary response {_response_cycle_counter}/{results_arr_len}: {_get_uniprot_identifier_json_nth_response(response,0)}")
    user_logic = input("Press 1 to confirm current result, 2 to cycle another result or 0 to discard.")
    logging.debug(f"[get_uniprot_identifier]: {gene_name} -> {uniprot_gene_identifier}")
    if prefix != "":
        return prefix+uniprot_gene_identifier
    else:
        return uniprot_gene_identifier

def get_uniprotId_from_geneName(gene_name, recursion=0, prefix="UniProtKB:", trust_genes=True):
    """
    Retrieves uniprot identifier from a gene symbol/name; e.g. UniProtKB:Q86SQ4, if gene_name=adgrg6

    Parameters:
      - gene_name: A gene name or symbol e.g. ADGRG6
      - recursion: An internal function parameter to iterate through an array of UniProtKB ids supplied by the response json
      - prefix: A database prefix, is prefixed before the specifi uniprot gene id
      - trust_genes: If True, all trusted genes inside trusted_genes.txt will override this function. If false, it
        notifies you that the gene was found among trusted genes and if you wish to proceed with the function.
    """
    if gene_name in constants.TRUSTED_GENES:
        if trust_genes == True:
            # return the element ahead of the gene_name, which is the previously found uniprot_id
            return "UniProtKB:" + constants.TRUSTED_GENES[constants.TRUSTED_GENES.index(gene_name)+1]
        else:
            trust_genes_user_response = int(input(f"Gene {gene_name} is found among trusted genes. Press 1 to continue with UniProt Id query, or 0 to use the trusted gene."))
            if trust_genes_user_response == 0:
                return "UniProtKB:" + constants.TRUSTED_GENES[constants.TRUSTED_GENES.index(gene_name)+1]

    global _uniprot_identifier_query_result
    uniprot_gene_identifier = ""
    if recursion == 0:
        _uniprot_identifier_query_result = requests.get(f"https://rest.uniprot.org/uniprotkb/search?query={gene_name}+AND+organism_id:9606&format=json&fields=id,gene_names,organism_name")
        _uniprot_identifier_query_result = _uniprot_identifier_query_result.json()

    uniprot_gene_identifier = _uniprot_identifier_query_result["results"][recursion]["primaryAccession"]
    results_arr_len = len(_uniprot_identifier_query_result["results"])

    logging.debug(f"Gene name {gene_name} found to correspond to {uniprot_gene_identifier}. Displaying response {recursion+1}/{results_arr_len}: {_get_uniprot_identifier_json_nth_response(_uniprot_identifier_query_result,recursion)}")
    user_logic = int(input("Press 1 to confirm current result, 2 to cycle another result or 0 to continue the program and discard all options."))
    if user_logic == 1:
        # result is confirmed, save so in TRUSTED_GENES
        logging.debug(f"[get_uniprot_identifier]: Confirmed {gene_name} -> {uniprot_gene_identifier}")
        with open("src_data_files/trusted_genes.txt", "a+") as f:
            f.seek(0) # a+ has file read pointer at bottom, f.seek(0) returns to top
            if gene_name not in f.read():
                f.seek(0,2) # move file pointer back to the end of the file
                f.write(f"{gene_name} {uniprot_gene_identifier}\n")
                f.close()
    elif user_logic == 2:
        # cycle another result
        next_recursive_step = recursion + 1
        if next_recursive_step > (results_arr_len-1): 
            # cycled through all options, return error
            logging.debug("Cycled out of options!")
            return f"[get_uniprot_identifier_new]: CycleOutOfBoundsError: Cycled through all the uniprot gene IDs for {gene_name} without confirming any"
            # TODO: major bug here, it should return error, but it returns a gene id instead
            # look at: https://stackoverflow.com/questions/11356168/return-in-recursive-function
            # and https://stackoverflow.com/questions/23543485/python-recursive-function-executing-return-statement-incorrectly 
            # TODO: handle this error in main script file
        get_uniprotId_from_geneName(gene_name, recursion=next_recursive_step)
    elif user_logic == 0:
        return f"[get_uniprot_identifier_new]: No uniprot gene IDs for {gene_name} found."
    else:
        # wrong input, recurse again
        logging.debug("[get_uniprot_identifier_new]: Wrong input! Must be 0, 1 or 2.")
        get_uniprotId_from_geneName(gene_name, recursion)

    if prefix != "":
        return prefix+uniprot_gene_identifier
    else:
        return uniprot_gene_identifier

def _get_uniprot_identifier_json_nth_response(json, nth_response):
    "Gets the entire nth element in 'results' array inside the json retrieved by get_uniprot_identifier function"
    return json["results"][nth_response]

def load_trusted_genes(trusted_genes_file_path):
    """
    Loads constants.TRUSTED_GENES list with genes from trusted_genes.txt
    """
    file = open(trusted_genes_file_path, "r")
    lines = file.readlines()
    for line in lines:
        splitlist = line.split(" ")
        for element in splitlist:
            constants.TRUSTED_GENES.append(element)
    
def sort_list_of_dictionaries(input, field, direction_reversed = True):
    """Sorts the list of dictionaries by the key "field", default direction is reversed (descending)"""
    return sorted(input, key=lambda d: d[field], reverse=direction_reversed) 

