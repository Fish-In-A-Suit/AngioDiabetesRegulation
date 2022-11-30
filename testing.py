import util
import json
import xml.etree.ElementTree as ET
import sys

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

    util.extract_n_elements_from_json("gene_scores/test_score_homosapinesonly=false,v2-term_enums,cross_section.json", 10, "gene_scores/test_score_homosapinesonly=false,v2-term_enums,cross_section,top10.json")


if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    main()
