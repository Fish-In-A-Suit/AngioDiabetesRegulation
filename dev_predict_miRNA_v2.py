# This script provides functionality to query and save miRNAs from databases and to analyse/compare the saved miRNAs
# against mRNA sequences.

# Algorithm
#   1. Find and save miRNAs
#   2. Get an array of mRNAs
#   3. Compare the mRNAs against miRNAs

import util
import constants
import os

from Bio import Entrez
from Bio import SeqIO

import logging
logger = logging.getLogger(__name__)

FLAG_EXCLUDE_PREDICTED_MRNAS = True # flag to exclude all XM_xxxxx NCBI Refseq access strings, which are experimentally computed and not yet confirmed

def find_miRNAs_with_match_strengths(miRDB_readlines, mRNA_refseq):
    """
    Finds all miRNAs and their match strengths (hsa-miR-xxx, 72.2) from miRDB_readlines for mRNA_refseq (eg. NM_xxxxx).

    Returns: list of lists (total list contains sublists, each sublist contains 2 elements - miRNA and match_strength)
    """
    result_list = []
    for line in miRDB_readlines:
        if mRNA_refseq.split(".")[0] in line: # mRNA refseqs are in the format NM_xxxxxx.a, in miRDB only NM_xxxxx is valid (strip .a)
            splitline = line.split("\t")
            miRNA = splitline[0]
            match_strength = float(splitline[2].replace("\n", ""))
            result_list.append([miRNA, match_strength]) # append a sublist [miRNA, match_strength] to result_list
    
    # order result_list by the highest match_strength
    # The sorted function takes the lists as an argument and sorts it using a lambda function as the key argument. 
    # The lambda function takes an element x and returns x[1], which is the second element in each sublist. 
    # The sorted function then sorts the list of lists by the value returned by the lambda function
    result_list = sorted(result_list, key = lambda x:x[1], reverse=True)

    return result_list


def main():
    logger.debug("Loading product mRNAs.")
    # handle the product_mRNAs_refseq_file initialisation
    product_mRNAs_refseq_file = constants.TARGET_FOLDER + "/product_mRNA_refseq.json" # default path
    if not os.path.exists(product_mRNAs_refseq_file): # if default path not available, try to get it from user input
        user_input = input(f"File {product_mRNAs_refseq_file} was not found. Check which constants.TARGET_FOLDER is specified. Press 0 to provide a custom filepath, 1 to abort or 2 to create a new product_mRNA_NCBIacc.json file.")
        if user_input == 0:
            user_input = input("Type the case-sensitive relative filepath: ")
            product_mRNAs_refseq_file = user_input
            logger.info(f"product_mRNAs_refseq_file set to {product_mRNAs_refseq_file}")
        elif user_input == 1:
            logger.info("Aborted.")
            return
        elif user_input == 2:
            util.product_mRNA_json_append_refseqIds()
            logger.info("Created new product_mRNA_NCBIacc.json file.")
        else:
            logger.info("Incorrect data specified (only 0 or 1 available). Aborting.")
            return
    else:
        logger.info(f"File {product_mRNAs_refseq_file} found.")

    product_mRNAs_refseq_scored_json = util.read_file_as_json(product_mRNAs_refseq_file)
    miRDB_miRNA_matching_results = util.readlines(constants.FILE_miRDB_miRNAs)

    # read all refseq ids (NM_xxxxx) and save them to corresponding uniprot ids
    mRNA_to_refseq_dict = {}
    for element in product_mRNAs_refseq_scored_json:
        mRNA_to_refseq_dict[element["UniprotID"]] = element["RefSeq_NT_IDs"]
    
    if FLAG_EXCLUDE_PREDICTED_MRNAS:
        mRNA_to_refseq_dict = util.remove_dict_list_elements(mRNA_to_refseq_dict, "XM")
    
    logger.debug(mRNA_to_refseq_dict)

    refseqs_to_miRNAs_dict = {} # "NM_xxxxxx": [[hsa-miR-xxx, 72.21], [...], ...], "NM_aaaaaaa": [[...], ...]
    k = 0
    for key in mRNA_to_refseq_dict:
        logger.info(f"Processing {k}/{len(mRNA_to_refseq_dict)}")
        mRNA_refseqs = mRNA_to_refseq_dict[key]
        for mRNA_refseq in mRNA_refseqs: # loop over all NM_xxxx for each mRNA
            mRNA_miRNA_matches = find_miRNAs_with_match_strengths(miRDB_miRNA_matching_results, mRNA_refseq)
            refseqs_to_miRNAs_dict[mRNA_refseq] = mRNA_miRNA_matches
        k += 1
    

    i = 0
    for element in product_mRNAs_refseq_scored_json: # each element is a single UniProtKB (from product_mRNA_refseq.json)
        logger.info(f"Processing {i}/{len(product_mRNAs_refseq_scored_json)}")

        current_element_refseqs = util.remove_list_elements(element["RefSeq_NT_IDs"],"XM") # retrieves all NM_xxxx refseqs of the current element
        current_refseqs_to_miRNAs_dict = {}
        # fill up current_refseqs_to_miRNAs_dict (since the "full" refseqs_to_miRNAs_dict contains ALL refseqs of all mRNAs and you only need matches for the current refseqs for the current uniprot element)
        for key in refseqs_to_miRNAs_dict:
            if key in current_element_refseqs:
                current_refseqs_to_miRNAs_dict[key] = refseqs_to_miRNAs_dict[key] # copy into current_refseqs_to_miRNAs_dict
        
        product_mRNAs_refseq_scored_json[i]["miRNA_matches"] = current_refseqs_to_miRNAs_dict
        i += 1 
        current_refseqs_to_miRNAs_dict = {} # reset current dict
    
    # save the modified product_mRNA_refseq_scored_json
    util.save_json(product_mRNAs_refseq_scored_json, constants.TARGET_FOLDER+"/product_mRNA_miRDB-predicted-miRNA-matching")

    return

if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    main()