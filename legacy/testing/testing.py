import legacy.util as util
import json
import xml.etree.ElementTree as ET
import sys
import legacy.constants as constants
import os

import logging
logger = logging.getLogger(__name__)

def find_products_related_to_GO_terms_new(terms, GO_term):
    def _handle_load_from_crash():
        nonlocal terms
        terms = util.list_directionshrink(terms, GO_term, forward=True)
    
    logger.debug(f"terms len before crash handle: {len(terms)}")
    _handle_load_from_crash()
    logger.debug(f"terms len after crash handle: {len(terms)}")
            

def main():
    # test_json = util.read_file_as_json("term_genes_crash/product-search-crash_1669038434.568223_.json")
    # logger.info(test_json[-1]["GO_term"])

    # terms_all = util.get_array_terms("ALL")
    # find_products_related_to_GO_terms_new(terms_all, "GO:0043534")

    # myTree = ET.parse("src_data_files/uniprot_sprot_human.xml")
    # i = 0
    # for element in myTree.iter("entry"):
    #    logger.debug(element.attrib)
    #    if i == 10:
    #        sys.exit()
    #    i+=1

    
    
    #util.replace_in_file(os.path.join(constants.TARGET_FOLDER, "input.txt"), os.path.join(constants.TARGET_FOLDER, "input1.txt"), " ", "\t")
    #util.generate_go_names_in_input_file(os.path.join(constants.TARGET_FOLDER, "input1.txt"))
    #return 0 

    start_time = util.get_time()
    # util.get_uniprotids_from_json("gene_scores/test_score_homosapinesonly=false,v2-term_enums,cross_section,top10.json")[0]
    # util.get_identifier_values_from_json("gene_scores/test_score_homosapinesonly=false,v2-term_enums,cross_section,top10.json","gene")[0]
    end_time = util.compute_time(start_time)
    logger.debug(end_time)

    # CALL THESE FUNCTIONS TO SAVE CUDA MRNA-MIRNA MATCH STRENGTH RESULTS AS JSON
    # sequence_comparison_results_json = util.load_sequence_comparison_results("src_data_files/sequence_comparison_results.txt")
    # util.save_json(sequence_comparison_results_json, os.path.join(constants.TARGET_FOLDER, "sequence_comparison_results.json"))

    # util.save_mirbase_hsap_miRNAs_for_cpp("src_data_files/miRNAdbs/mirbase_miRNA_06-02-2023.dat")

    #
    # match_strength_test = util.compare_miRNA_mRNA_match_strength_single("ATTTC", "TATAG")
    # logger.debug(f"match strength is: {match_strength_test}")

    util.init_CUDA_sequence_comparison_results_json()
    util.init_MI_HSA_miRNA_mapping()
    util.init_mRNA_and_miRNA_sequence_dicts()

    miRNA_id = util.map_mi_hsa("", "hsa-let-7f-1", 1)
    logger.debug(f"miRNA_id = {miRNA_id}")

    # this code showcases how to access elements in constants.CUDA_sequence_comparison_results_json.
    # i = 0
    # logger.debug(f"-------")
    # for dict_element in constants.CUDA_sequence_comparison_results_json:
    #    if i == 0:
    #        logger.debug(f"{dict_element['MI0000060'][1]['UniProtKB:Q0VGL1']}")
    #    i += 1

    # logger.debug(util.find_CUDA_miRNA_mRNA_match_strength('hsa-let-7f-1', 'UniProtKB:Q0VGL1'))


    logger.debug(util.compare_python_CUDA_miRNA_mRNA_match_strength("MI0000060", 'UniProtKB:Q0VGL1'))
    logger.debug(util.compare_python_CUDA_miRNA_mRNA_match_strength("MI0000060", 'UniProtKB:Q96E22'))

    logger.debug("---")

    logger.debug(util.compare_miRNA_mRNA_match_strength_single_v2("MI0000060", 'UniProtKB:Q0VGL1', debugLog=True, fileWrite=True))

    total_ops = len(constants.miRNA_id_sequence_dict) * len(constants.mRNA_id_sequence_dict)
    i = 0

    writelines = []
    for miRNA_id, miRNA_sequence in constants.miRNA_id_sequence_dict.items():
        logger.debug(f"{i}/{len(constants.miRNA_id_sequence_dict)}: {miRNA_id}")
        writelines.append(f"{miRNA_id}:\n")
        i+=1
        for mRNA_id, mRNA_sequence in constants.mRNA_id_sequence_dict.items():
            match_strength = util.compare_miRNA_mRNA_match_strength_single_v2(miRNA_sequence, mRNA_sequence, debugLog=False)
            writelines.append(f"\t{mRNA_id}: {format(match_strength, '.4f')}\n")
        # write all mRNAs for this miRNA
        with open("python_sequence_comparison_score.txt", "a") as f:
            for line in writelines:
                f.write(line)
        writelines = [] 

    """ # slower way
    with open("python_sequence_comparison_score.txt", "a") as f:
        for miRNA_id, miRNA_sequence in constants.miRNA_id_sequence_dict.items():
            logger.debug(f"{miRNA_id}:")
            f.write(f"{miRNA_id}:\n")
            for mRNA_id, mRNA_sequence in constants.mRNA_id_sequence_dict.items():
                match_strength = util.compare_miRNA_mRNA_match_strength_single_v2(miRNA_sequence, mRNA_sequence, debugLog=False)
                logger.debug(f"    - {mRNA_id}: {match_strength}")
                f.write(f"\t{mRNA_id}: {format(match_strength, '.4f')}\n")
                i+=1
    """

    #string = "Word"
    #string_list = [*string]
    #logger.info(string_list)


if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    main()
