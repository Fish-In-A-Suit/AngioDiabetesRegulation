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

def _overlap_substring_on_mRNAs(substring, mRNAs):
    #TODO:exception if length of substring and all mRNA substrings is not same
    logger.debug(f"Starting overlap analysis of substring {substring}")
    miRNA_results=[]
    for i in range(len(mRNAs)):
        #logger.debug(f"Trying {i}. mRNAsubsequences: {mRNAssubstrings[i]}")
        miRNA_results.append(_overlap_substring_on_mRNA(substring, mRNAs[i]))
    return miRNA_results

def _overlap_substring_on_mRNA(substring, mRNA):
    """
    Overlaps substring over mRNA (also a string), returns overlap score as levenstein distance
    """
    
    mRNA_length = len(mRNA)
    substring_length = len(substring)
    # logger.debug(f"Overlapping substring {substring} over type {type(mRNA)}")
    overlap = 1-(jellyfish.levenshtein_distance(substring, mRNA) - max(mRNA_length-substring_length, 0)) / substring_length
    return overlap

def _score_miRNA(overlap_result, productnames, product_score):
    """
    Computes a final score for a miRNA sequence from all possible overlap scores.

    Parameters:
      - overlap_result: a json of all overlap results (obtained after string overlapping onto mRNAs) for a given miRNA
                        in the form of 
        {
            "miRNA": ATGCCTA, 
            "miRNA_overlaps": {
               ["productID":aaa, "overlap":a], 
               ["productID":bbb, "overlap":b]
            }
        }
      
      - productnames: a list of uniprotid productnames obtained by util.get_identifier_values_from_json(mrna_filepath, "UniprotID")[0]
      - product_score: a list of uniprotid product scores (same length as productnames just contains their respective al_corr scores)
    """
    score = 0
    for overlap_element in overlap_result["miRNA_overlaps"]:
        current_productname = overlap_element["productID"]
        overlap = overlap_element["overlap"]
        score_index = productnames.index(current_productname)

        score += overlap * product_score[score_index]
    return score


def generate_all_miRNA_permutations_with_repetitions(length):
    """
    Generates all the possible miRNA permutations (with repetitions (eg AATT) included) of the specified length.
    """
    nucleotides = ["A", "T", "C", "G"]
    # https://stackoverflow.com/questions/3099987/generating-permutations-with-repetitions
    permutations = list(map("".join, itertools.product(nucleotides, repeat=length)))
    logger.debug(f"There are {len(permutations)} permutations with length {length}")
    return permutations

def predict_miRNAs(productnames, mRNAs, product_scores, length, treshold_to_accept, target_folder, debug_loop_break=1):
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
    
    Parameters:
      - productnames: a list of uniprotIds obtained by util.get_identifier_values_from_json(mrna_filepath, "UniprotID")[0]
      - mRNAs: a list of all mRNA sequences obtained by util.get_identifier_values_from_json(mrna_filepath, "mRNA")[0]
      - product_scores: a list of all respective uniprotId product scores
      - length: the length of miRNA permutations
      - threshold_to_accept: the minimal percentage of miRNA length that must fit to a mRNA to be considered as a possible fit
      - target_folder: where to save miRNA_overlap.json
      - debug_loop_break: used for breaking the top level for loop after debug_loop_break iterations
    """

    logger.info(f"Starting brute-force prediction of miRNA (length {length}) overlaps for (productID, score): {list(map(lambda i, s: (i, s), productnames, product_scores))}")

    t = util.get_time()
    miRNA_candidates = generate_all_miRNA_permutations_with_repetitions(length)
    logger.debug(f"Permutation calculation: {util.compute_time(t)} seconds. Size of object = {util.get_size(miRNA_candidates)}")

    json_overlap_results = [] # this is the final result that is saved
    logger.info(f"Total miRNA candidates: {len(miRNA_candidates)}")
    start_time = util.get_time()
    for i in range(len(miRNA_candidates)): # for each miRNA
        if(i >= debug_loop_break): 
            logger.debug(f"{i} iterations took {util.compute_time(start_time)}")
            break # breaks out of the top-level for loop after debug_loop_break iterations
        miRNA = miRNA_candidates[i]
        logger.info(f"Overlaping miRNA -> {miRNA} ({i}/{len(miRNA_candidates)})")
        
        # fills up after the mRNA for loop
        # each element represents overlap scores of all mRNAs for a single miRNA
        # mirna_overlap_results[1] = {"productID": xxx, "overlap": x} ...

        mirna_overlap_results=[] # fills up after the mRNA for loop

        for j in range(len(mRNAs)): # for each mRNA
            mRNA = mRNAs[j]
            if mRNA == None: continue
            overlap = _overlap_substring_on_mRNA(miRNA, mRNA) # levenstein score
            if overlap >= treshold_to_accept:
                # TODO: in the future, maybe accept, score and save all predicted miRNAs, then check which of the predicted miRNA sequences already exist as discovered and what their predicted score is
                mirna_overlap_results.append({"productID" : productnames[j], "overlap" : overlap})
        
        # each element in json_overlap_results represents a single miRNA and it's scoring results (obtained by overlapping it over mRNAs)
        # {
        #   "miRNA": ATGCCTA, 
        #   "miRNA_overlaps": {
        #       ["productID":aaa, "overlap":a], 
        #       ["productID":bbb, "overlap":b]
        #   }
        # }
        temp_json = { "miRNA" : miRNA, "miRNA_overlaps" : mirna_overlap_results} #TODO: what if mirna_overlap_results is []??
        # TODO: maybe compute _score_miRNA here to avoid doubling for loops ?
        json_overlap_results.append(temp_json)

    for i in range(len(json_overlap_results)):
        mirna_score = _score_miRNA(json_overlap_results[i], productnames, product_scores)
        
        json_overlap_results[i]["miRNA_score"] = mirna_score

    json_overlap_results = util.sort_list_of_dictionaries(json_overlap_results, "miRNA_score")
    util.save_json(json_overlap_results, os.path.join(target_folder, "miRNA_overlap.json"))
    logger.info (f"Done with predicting miRNA!")

def main():
    constants.PERMUTATIONS_LEN_FOUR = util.get_permutations_with_repetitions(4, ["A","T","C","G"])
    logger.debug(constants.PERMUTATIONS_LEN_FOUR)
    
    mrna_filepath = os.path.join(constants.TARGET_FOLDER, "product_mRNA.json")
    #product_list = ["UniProtKB:Q16613", "UniProtKB:O15530", "UniProtKB:Q9Y243"]
    product_list = util.get_identifier_values_from_json(mrna_filepath, "UniprotID")[0]

    mRNAs = util.get_identifier_values_from_json(mrna_filepath, "mRNA")[0] #mRNAs = ["ABABABABABABABABABAB","ABCABCABABABABABABABABABABABCABC","ABCABCABCABCABCABC"]
    
    # in product_mRNA.json that includes all mRNAs, found 2 mRNAs with null
    #_d_mRNAs_types = []
    #for i,mRNA in enumerate(mRNAs):
    #    _d_mRNAs_types.append(f"i={i} {type(mRNA)}")
    #logger.debug(f"mRNAs types (any Nones -> check algorithm for possible errors) = {_d_mRNAs_types}")

    scores_json = util.read_file_as_json(os.path.join(constants.TARGET_FOLDER, "product_scores.json"))
    scores = []
    for product in product_list:
        product_scores_index = next((index for (index, d) in enumerate(scores_json) if d["product"] == product), None)
        scores.append(scores_json[product_scores_index]["al_corr_score"])
    logger.debug(f"Found {len(mRNAs)} mRNAs with the following {len(scores)} product scores: {scores}")

    #product_scores_index = next((index for (index, d) in enumerate(product_scores) if d["product"] == productID), None)
    #corr_score = product_scores[product_scores_index]["al_corr_score"]   
    # 
    #product_list = ["UniProtKB:Q16613", "UniProtKB:O15530", "UniProtKB:Q9Y243"]
    #mRNAs = ["GAACATCTGCTACTACAGCCTTGCAGCCCGGAGTCCCGGATTTTACTGGTTCCCGTGCCTGCGGACAGGC","ATTGCTGGGGCTCCGCTTCGGGGAGGAGGACGCTGAGGAGGCGCCGAGCCGCGC","CCAAACCCTAAAGCTGATATCACAAAGTACCATTTCTCCAAGTTGGGGGCTCAGAGGGGAGTCATCATGAGCGA"]
    #scores = [0.5, -1, 1]
    
    predict_miRNAs(product_list, mRNAs, scores, length=13, treshold_to_accept=0.5, target_folder=constants.TARGET_FOLDER) 
    #WARNING: this is a brute force method pathfinder, with extensive debug output! do not use with length more than 5.
    #It will create large log files and take up a significant amount of ram!

if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    main()
