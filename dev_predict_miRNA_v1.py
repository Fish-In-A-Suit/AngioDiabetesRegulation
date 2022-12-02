import util
import constants
import requests
import jellyfish
import os
import itertools

import logging
logger = logging.getLogger(__name__)

def _find_all_unique_substrings_from_strings(mRNAs, length):
    logger.debug(f"Finding all unique substrings with length {length} from {len(mRNAs)} mRNAs.")
    substrings = set()
    for mRNA in mRNAs:
        substrings.update(_find_all_unique_substrings_from_string(mRNA, length))
    logger.debug(f"Found {len(substrings)} unique substrings: {substrings}")
    return substrings

def _find_all_unique_substrings_from_string(mRNA, length):
    substrings = set()
    for i in range(len(mRNA)-length):
        substrings.update([mRNA[i:i+length]])
    logger.debug(f"Recieved string {mRNA}  ->  {substrings}")
    return list(substrings)

def find_indices_of_substring(full_string, sub_string):
    return [index for index in range(len(full_string)) if full_string.startswith(sub_string, index)]

def _overlap_substring_on_mRNAs(substring, mRNAssubstrings, treshold_to_accept):
    #TODO:exception if length of substring and all mRNA substrings is not same
    logger.debug(f"Starting overlap analysis of substring {substring}")
    miRNA_results=[]
    for i in range(len(mRNAssubstrings)):
        logger.debug(f"Trying {i}. mRNAsubsequences: {mRNAssubstrings[i]}")
        miRNA_results.append(_overlap_substring_on_mRNAsubstrings(substring, mRNAssubstrings[i], treshold_to_accept))
    return miRNA_results

def _overlap_substring_on_mRNAsubstrings(substring, mRNAsubstrings):
    #logger.debug(f"Overlaping {substring} on {mRNAsubstrings}")
    length = len(substring)
    temp_miRNA_results = []
    for j in range(len(mRNAsubstrings)):
            overlap = 1-(jellyfish.levenshtein_distance(substring, mRNAsubstrings[j]) / length)
            #logger.debug(f"Overlap of substring {substring} on mRNAsubstring {mRNAsubstrings[j]} is {overlap}")
            temp_miRNA_results.append([mRNAsubstrings[j], overlap])
    return temp_miRNA_results

def _find_best_overlap_miRNA_on_mRNA(miRNA, mRNA):
    if mRNA == None:
        return 0
    mRNAsubstrings = _find_all_unique_substrings_from_string(mRNA, len(miRNA))
    miRNA_overlaps = _overlap_substring_on_mRNAsubstrings(miRNA, mRNAsubstrings)
    only_overlaps = [element[1] for element in miRNA_overlaps]
    best_overlap = max(only_overlaps)

    logger.debug(f"Best overlap of miRNA ({miRNA}) on mRNA ({mRNA}) is {best_overlap})")
    return best_overlap

def _score_miRNA(overlap_result, productnames, product_score):
    score = 0
    for overlap_element in overlap_result["miRNA_overlaps"]:
        current_productname = overlap_element["productID"]
        overlap = overlap_element["overlap"]
        score_index = productnames.index(current_productname)

        score += overlap * product_score[score_index]
    return score


def _generate_all_miRNA_permutations_with_repetitions(lenght):
    nucleotides = ["A", "T", "C", "G"]
    permutations = list(map("".join, itertools.product(nucleotides, repeat=lenght)))
    logger.debug(f"There are {len(permutations)} permutations with length {lenght}")
    return permutations

def predict_miRNAs(productnames, mRNAs, product_scores, length, treshold_to_accept, target_folder):
    """
    Aljosa&Ladi algorithm for finding miRNA of certain length from a list of mRNA.
    Because this is a scientific program, return all possible predictions with some sort of score.

    How it works: 
    1   find all possible substrings with length=(10) in all mRNA sequences and combine then in a set (to remove duplicates)
    1.5 add variability to substrings by replacing up to treshold_to_accept % of elements in each substring with random elements. 
    (1 and 1.5 ensure we have all possible substrings that have at least treshold_to_accept % fit with at least one mRNA sequence)
    2   in every mRNA sequence try to find each possible substring or different 
        substring that is at least treshold_to_accept similar to the wanted.

    treshold_to_accept - what is the minimal percentage of the miRNA length that must fit to a mRNA in order to be considered as 'fitting'
    """

    logger.info(f"Starting brute-force prediction of miRNA (length {length}) overlaps for (productID, score): {list(map(lambda i, s: (i, s), productnames, product_scores))}")

    miRNA_candidates = _generate_all_miRNA_permutations_with_repetitions(length)

    json_overlap_results = []
    for i in range(len(miRNA_candidates)):
        miRNA = miRNA_candidates[i]
        logger.info(f"Overlaping miRNA -> {miRNA} ({i}/{len(miRNA_candidates)})")

        mirna_overlap_results=[]

        for j in range(len(mRNAs)):
            mRNA = mRNAs[j]
            
            overlap = _find_best_overlap_miRNA_on_mRNA(miRNA, mRNA)
            if overlap >= treshold_to_accept:
                mirna_overlap_results.append({"productID" : productnames[j], "overlap" : overlap})
        
        temp_json = { "miRNA" : miRNA, "miRNA_overlaps" : mirna_overlap_results} #TODO: what if mirna_overlap_results is []??
        json_overlap_results.append(temp_json)

    for i in range(len(json_overlap_results)):
        mirna_score = _score_miRNA(json_overlap_results[i], productnames, product_scores)
        
        json_overlap_results[i]["miRNA_score"] = mirna_score

    json_overlap_results = util.sort_list_of_dictionaries(json_overlap_results, "miRNA_score")
    util.save_json(json_overlap_results, os.path.join(target_folder, "miRNA_overlap.json"))
    logger.info (f"Done with predicting miRNA!")

def main():
    mrna_filepath = os.path.join(constants.TARGET_FOLDER, "product_mRNA.json")
    #product_list = ["UniProtKB:Q16613", "UniProtKB:O15530", "UniProtKB:Q9Y243"]
    product_list = util.get_identifier_values_from_json(mrna_filepath, "UniprotID")[0]

    mRNAs = util.get_identifier_values_from_json(mrna_filepath, "mRNA")[0] #mRNAs = ["ABABABABABABABABABAB","ABCABCABABABABABABABABABABABCABC","ABCABCABCABCABCABC"]
    #predicted_miRNAs = predict_miRNAs(mRNAs, 10, 0.5)

    scores_json = util.read_file_as_json(os.path.join(constants.TARGET_FOLDER, "product_scores.json"))
    scores = []
    for product in product_list:
        product_scores_index = next((index for (index, d) in enumerate(scores_json) if d["product"] == product), None)
        scores.append(scores_json[product_scores_index]["al_corr_score"])


    #product_scores_index = next((index for (index, d) in enumerate(product_scores) if d["product"] == productID), None)
    #corr_score = product_scores[product_scores_index]["al_corr_score"]   
    # 
    #product_list = ["UniProtKB:Q16613", "UniProtKB:O15530", "UniProtKB:Q9Y243"]
    #mRNAs = ["GAACATCTGCTACTACAGCCTTGCAGCCCGGAGTCCCGGATTTTACTGGTTCCCGTGCCTGCGGACAGGC","ATTGCTGGGGCTCCGCTTCGGGGAGGAGGACGCTGAGGAGGCGCCGAGCCGCGC","CCAAACCCTAAAGCTGATATCACAAAGTACCATTTCTCCAAGTTGGGGGCTCAGAGGGGAGTCATCATGAGCGA"]
    #scores = [0.5, -1, 1]
    
    predict_miRNAs(product_list, mRNAs, scores, length=4, treshold_to_accept=0.5, target_folder=constants.TARGET_FOLDER) 
    #WARNING: this is a brute force method pathfinder, with extensive debug output! do not use with length more than 5.
    #It will create large log files and take up a significant amount of ram!

if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    main()
