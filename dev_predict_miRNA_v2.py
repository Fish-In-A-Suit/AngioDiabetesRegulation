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
    product_miRNA_matching_results = util.readlines(constants.FILE_miRDB_miRNAs)

    # TODO: optimisation suggestion: convert constants.FILE_miRDB_miRNAs file into a json file, where each element would be
    # accessed by the mRNA target product id (eg. NM_000842) and it would hold the references to all hsa-miR-xxxx miRNAs and the respective
    # match strenghts.
    # 
    # [
    #   {
    #   "product_mRNA_id": "NM_000842"
    #   "miRNAs": {
    #        "miRNAid": "hsa-miR-xxx"
    #        "matchStrength": 75.9
    #       },
    #       {
    #        ...
    #       },
    #       ...
    #   }
    # ]

    return

if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    main()