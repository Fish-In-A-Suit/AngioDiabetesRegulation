import util
import constants
import json
import sys
import os
import math

import logging
logger = logging.getLogger(__name__)

def _map_enum(enum_longname):
    """
    Used in _score_term_enums to map longer term enums to shorter enums to not clutter the json
    """
    if enum_longname == "terms_angiogenesis_positive": return "nterms_angio+"
    elif enum_longname == "terms_angiogenesis_negative": return "nterms_angio-"
    elif enum_longname == "terms_angiogenesis_general": return "nterms_angio0"
    elif enum_longname == "terms_diabetes_positive": return "nterms_dia+"
    elif enum_longname == "terms_diabetes_negative": return "nterms_dia-"
    elif enum_longname == "terms_diabetes_general": return "nterms_dia0"
    elif enum_longname == "terms_diabetes_all": return "nterms_dia"
    elif enum_longname == "terms_angiogenesis_all": return "nterms_angio"
    elif enum_longname == "notfound": return "notfound"

def _score_term_enums(terms, process_list, term_list):
    """
    Finds which arrays the terms belong to in constants.py
    Returns the following set:
    {
        nterms_angio: x, nterms_dia: x, 
        nterms_angio+: x, nterms_angio-: x, nterms_angio0: x,
        nterms_dia+: x, nterms_dia-: x, nterms_dia0: x
    }
    """
    result_set ={}
    for process in process_list: #define the result_set
        result_set[f"nterms_{process['name']}"] = 0
        result_set[f"nterms_{process['name']}+"] = 0
        result_set[f"nterms_{process['name']}-"] = 0
        result_set[f"nterms_{process['name']}0"] = 0

    
    #result_set = {
    #    "nterms_angio": 0,
    #    "nterms_dia": 0,
    #    "nterms_angio+": 0,
    #    "nterms_angio-": 0,
    #    "nterms_angio0": 0,
    #    "nterms_dia+": 0,
    #    "nterms_dia-": 0,
    #    "nterms_dia0": 0
    #}
        
    for term in terms:
        term_entry = next((item for item in term_list if item["GO_id"] == term), None) #search the list of dictionaries and return the index
        if term_entry is None:
            continue
        else:
            result_set[f"nterms_{term_entry['process']}"] += 1
            result_set[f"nterms_{term_entry['process']}{term_entry['direction']}"] += 1
    logger.debug(f"term_enum_scores: {result_set}")
    return result_set

def score_genes_v2(product_list, process_list, term_list, destination_filename="product_scores.json" , source_folder=constants.TARGET_FOLDER, use_cross_section=False):
    """
    Counts the number of appearances of all the genes across all specified json_files (which contain genes
    related to specific GO terms and are made by the find_genes_related_to_GO_terms function)

    Parameters:
      - TODO: explain others
      - use_cross_section: If True, only genes that have both nterms_angio >= 1 and nterms_dia >= 1 are appended to final json.
                           If False, all genes regardless of having bother angiogenetic or diabetic influence are appended to final json.
    """
    
    #source_filepath = os.path.join(source_folder, "terms_direct_products.json")
    #logger.info(f"Finding all genes from file: {source_filepath}")
    gene_set = set() # Set data type used instead of List, because a Set cannot have multiple occurences of the same element
    term_genes = [] # array of all genes across all terms; structure = {[term1, genes1], [term2, genes2], ... [term_n, genes_n]}
    #product_list = util.read_file_as_json(source_filepath)
    for term in product_list:
        if term["products"] != []:
            genes = [product["UniprotID"] for product in term["products"] if product["UniprotID"] is not None]
        else:
            genes = []
        term_genes.append([term["GO_term"], genes])
        gene_set.update(genes)
    logger.info(f"Found {len(gene_set)} different genes.")
    #logger.debug(f"[term_genes]: {term_genes}")
    
    json_gene_scores=[]
    for gene in gene_set:
        #destination_file = f"{destination_folder}/{gene}.json"
        gene_count_result = _score_gene_basic(gene,term_genes)
        term_enum_score_set = _score_term_enums(gene_count_result[1], process_list, term_list) # add nterms_angio, nterms_dia, nterms_angio+, nterms_angio-, nterms_angio0, nterms_dia+, nterms_dia-, nterms_dia0
        al_score = _score_gene_al(term_enum_score_set, process_list)
        al_e_score = _score_gene_al_e(term_enum_score_set, process_list)
        al_corr_score = _score_gene_al_corr(term_enum_score_set, process_list)
        out = {"product": gene, "count": gene_count_result[0], "terms:": gene_count_result[1], "term_enum_scores": term_enum_score_set, "al_score": al_score,  "al_e_score": al_e_score, "al_corr_score": al_corr_score}
        if use_cross_section == True:
            if all([(term_enum_score_set[f"nterms_{line['name']}"] != 0) for line in process_list]):
                json_gene_scores.append(out)
        else:
            json_gene_scores.append(out)
    logger.info(f"[gene_scores]: {json_gene_scores}")

    json_gene_scores = util.sort_list_of_dictionaries(json_gene_scores, "al_corr_score")
    util.save_json(json_gene_scores, os.path.join(source_folder, destination_filename))
    logger.info (f"Done with scoring!")


def score_genes_v1(allowed_term_ids, destination_file, source_folder="term_genes", use_cross_section=False):
    logger.info(f"Finding all genes from terms: {allowed_term_ids}")
    gene_set = set() # Set data type used instead of List, because a Set cannot have multiple occurences of the same element
    term_genes = [] # array of all genes across all terms; structure = {[term1, genes1], [term2, genes2], ... [term_n, genes_n]}
    for term in allowed_term_ids:
        if term.replace("-",":") in constants.TERMS_EMPTY: # some terms have 0 genes, don't process these
            continue
        genes = _import_genes_from_term_json(term, source_folder)
        term_genes.append([term, genes])
        gene_set.update(genes) # updates values in gene_set with genes, without duplicating existing elements
    
    logger.info(f"Found {len(gene_set)} different genes.")
    logger.debug(f"[term_genes]: {term_genes}")

    json_gene_scores=[]
    for gene in gene_set:
        #destination_file = f"{destination_folder}/{gene}.json"
        gene_count_result = _score_gene_basic(gene,term_genes)
        term_enum_score_set = _score_term_enums(gene_count_result[1]) # add nterms_angio, nterms_dia, nterms_angio+, nterms_angio-, nterms_angio0, nterms_dia+, nterms_dia-, nterms_dia0
        out = {"gene": gene, "count": gene_count_result[0], "terms": gene_count_result[1], "term_enum_scores": term_enum_score_set}
        if use_cross_section == True:
            if term_enum_score_set["nterms_angio"] != 0 and term_enum_score_set["nterms_dia"] != 0:
                json_gene_scores.append(out)
        else:
            json_gene_scores.append(out)
    logger.info(f"[gene_scores]: {json_gene_scores}")

    logger.info(f"[gene_scores]: {json_gene_scores}")
    json_gene_scores = util.sort_list_of_dictionaries(json_gene_scores, "count")
    util.save_json(json_gene_scores, destination_file)
    logger.info (f"Done with scoring!")

def _score_gene_basic(gene, term_genes_list):
    """
    Scores the 'gene' in question by finding how many times they occur in term_genes_list, and also finds which
    terms they correspond to.

    Parametes:
      - gene: a UniprotId of a specific gene
      - term_genes_list: a list of all terms and genes in the form of {[term1, genes1], [term2, genes2], ... , [term_n, genes_n]}
    
    Returns:
      - [0] = count: an int specifying how many times a gene was involved in term_genes_list
      - [1] = terms_involved_in: a list with all the terms a gene is involved in
    """
    # TODO
    # perhaps we will have multiple scores/indicators
    # scores=[]
    count=0
    terms_involved_in=[]
    for element in term_genes_list: # term_genes_list = {[term1,genes1], [term2,genes2],...}
        if gene in element[1]: # if gene in genes_n of term_n
            count += 1
            terms_involved_in.append(element[0]) # also append which term the gene was involved in
    logger.debug(f"Gene {gene} has {count} occurences.")
    return count, terms_involved_in

def _score_gene_al(term_scores, process_list):
    """
    With the current setting, top genes will be those which have positive angiogenetic and diabetic terms.
    Score is lower for genes with negative angiogenetic and diabetic terms.
    """
    wanted_sign = {}
    opposite_sign = {}
    for process in process_list:
        wanted_sign[process["name"]] = process["direction"]
        if wanted_sign == "+":
            opposite_sign[process["name"]] = "-"
        else:
            opposite_sign[process["name"]] = "+"

    score = (
        (sum([term_scores[f"nterms_{process['name']}{wanted_sign[process['name']]}"] for process in process_list]) * 10 +
        sum([term_scores[f"nterms_{process['name']}{opposite_sign[process['name']]}"] for process in process_list]) * -10) *
        sum([term_scores[f"nterms_{process['name']}0"] for process in process_list]) * 1 
    )
    return score

def _score_gene_al_e(term_scores, process_list):
    wanted_sign = {}
    opposite_sign = {}
    for process in process_list:
        wanted_sign[process["name"]] = process["direction"]
        if wanted_sign == "+":
            opposite_sign[process["name"]] = "-"
        else:
            opposite_sign[process["name"]] = "+"

    score = (
        (2**(sum([term_scores[f"nterms_{process['name']}{wanted_sign[process['name']]}"] for process in process_list])) * min([term_scores[f"nterms_{process['name']}{wanted_sign[process['name']]}"] for process in process_list]) -
        8**(sum([term_scores[f"nterms_{process['name']}{opposite_sign[process['name']]}"] for process in process_list])) * max([term_scores[f"nterms_{process['name']}{opposite_sign[process['name']]}"] for process in process_list])) *
        (sum([term_scores[f"nterms_{process['name']}0"] for process in process_list])) * 1
    )
    return score
    
def _score_gene_al_corr(term_scores, process_list):
    wanted_sign = {}
    opposite_sign = {}
    for process in process_list:
        wanted_sign[process["name"]] = process["direction"]
        if wanted_sign == "+":
            opposite_sign[process["name"]] = "-"
        else:
            opposite_sign[process["name"]] = "+"

    score = (
        min(math.prod([term_scores[f"nterms_{process['name']}"] for process in process_list]), 1) * (
        (min(math.prod([term_scores[f"nterms_{process['name']}{wanted_sign[process['name']]}"] for process in process_list]), 1) * max(1-(sum([term_scores[f"nterms_{process['name']}{opposite_sign[process['name']]}"] for process in process_list])), 0) + 
        max(-(math.prod([term_scores[f"nterms_{process['name']}{opposite_sign[process['name']]}"] for process in process_list])), -1) * max(1-(sum([term_scores[f"nterms_{process['name']}{opposite_sign[process['name']]}"] for process in process_list])), 0)) * 10 +
        (sum([term_scores[f"nterms_{process['name']}{wanted_sign[process['name']]}"]**(1/2) for process in process_list]) - sum([term_scores[f"nterms_{process['name']}{wanted_sign[process['name']]}"]**(1/2) for process in process_list]) )
        ) * max(sum([term_scores[f"nterms_{process['name']}0"]**(1/3) for process in process_list]), 1)
        )
    return score

def _import_genes_from_term_json(term, source_folder):
    """
    TODO
    """
    term_file = str(term).replace(":", "-")
    input_json = util.read_file_as_json(f"{source_folder}/{term_file}.json")
    logger.debug(input_json)
    genes = []
    for block in input_json:
        genes.append(block["product"])
    logger.debug(f"From file {source_folder}/{term_file}.json we found {len(genes)} genes: {genes}")
    return genes

def main():
    util.load_list_from_file(os.path.join(constants.TARGET_FOLDER, "terms_empty.txt"), constants.TERMS_EMPTY)
    #util.load_human_orthologs()

    product_list = util.read_file_as_json(os.path.join(constants.TARGET_FOLDER, "terms_direct_products.json"))
    inputs = util.read_input_file()
    process_list = inputs["processes"]
    terms_list = inputs["GO_terms"]
    
    score_genes_v2(product_list, process_list, terms_list, source_folder=constants.TARGET_FOLDER, use_cross_section=True)

    # dest_filename_v1 = "XXXXXX"
    # score_genes_v1(util.get_array_terms("ALL"), dest_filename_v1, source_folder="term_genes/homosapiens_only=false,v1", use_cross_section=True)

    #util.scoring_results_postprocess("gene_scores/test_score_homosapinesonly=false,v2-term_enums,cross_section.json")

    # util.load_json_by_terms("term_genes/homosapiens_only=false,v1", terms_all)
    # termfiles_angiogenesis = util.load_json_by_terms("term_genes/homosapiens_only=false,v1", util.get_array_terms("ANGIOGENESIS"))
    # termfiles_diabetes = util.load_json_by_terms("term_genes/homosapiens_only=false,v1", util.get_array_terms("DIABETES"))

    # util.find_term_corresponding_array("GO:0001525")
        

if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    main()

