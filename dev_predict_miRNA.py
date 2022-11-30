import util
import constants
import requests
import jellyfish
import os

import logging
logger = logging.getLogger(__name__)

def get_mrna(gene_list, to_be_inhibited, target_folder):
    
    ensembl_ids = util.get_ensembl_ids_from_uniprot_ids(gene_list)
    mRNAs = _get_ensembl_mRNA_sequences(ensembl_ids)

    json_mRNAs = []
    for i in range(len(gene_list)):
        out = {"UniprotID" : gene_list[i], "EnsemblID" : ensembl_ids[i], "ToBeInhibited":to_be_inhibited[i], "mRNA" : mRNAs[i]}
        json_mRNAs.append(out)
    util.save_json(json_mRNAs, os.path.join(target_folder, "product_mRNA.json"))
    return json_mRNAs

def _get_ensembl_mRNA_sequences(ensembl_ids):
    mRNAs=[]
    for id in ensembl_ids:
        mRNAs.append(_get_ensembl_sequence_API(id))
    return mRNAs

def _get_ensembl_sequence_API(id, type="cdna"):
    """
    This function queries ensembl for nucleotide sequence
    Input of ensembl ID's must be a string

    type: either cds (without UTR) or cdna (with UTR)
    """
    logger.debug(
        f"Starting Ensembl API for id {id}")
    response = requests.get(f"https://rest.ensembl.org/sequence/id/{id}?object_type=transcript;type={type}", headers={"Content-Type": "text/plain", }
                            )  # cds = c-DNA without the UTR regions; type=cdna (in this case UTR region is also kept); retrieves complementary sequence to the gene mRNA (without UTR), same as miRNA sequence (todo: preveri z genetikom)
    if response.ok:
        logger.debug(response.text)
        logger.info(
            f"Recieved sequence for id {id}.")
        return response.text
    else:
        logger.info(
            f"Failed to get sequence for id {id}")
        return None
    

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
    return substrings

def find_indices_of_substring(full_string, sub_string):
    return [index for index in range(len(full_string)) if full_string.startswith(sub_string, index)]

def _overlap_substring_on_mRNAs(substring, mRNAssubstrings, treshold_to_accept):
    #TODO:exception if length of substring and all mRNA substrings is not same
    logger.debug(f"Starting overlap analysis of substring {substring}")
    miRNA_results=[]
    for i in range(len(mRNAssubstrings)):
        logger.debug(f"Trying {i}. mRNAsubsequences: {mRNAssubstrings[i]}")
        miRNA_results.append(_overlap_substring_on_mRNA(substring, mRNAssubstrings[i], treshold_to_accept))
    return miRNA_results

def _overlap_substring_on_mRNA(substring, mRNAsubstrings, treshold_to_accept):
    logger.debug(f"Overlaping {substring} on {mRNAsubstrings}")
    length = len(substring)
    temp_miRNA_results = []
    for j in range(len(mRNAsubstrings)):
            overlap = 1-(jellyfish.levenshtein_distance(substring, mRNAsubstrings[j]) / length)
            logger.debug(f"Overlap of substring {substring} on mRNAsubstring {mRNAsubstrings[j]} is {overlap}")
            if overlap >= treshold_to_accept:
                temp_miRNA_results.append([mRNAsubstrings[j], overlap])
    return temp_miRNA_results

def predict_miRNAs(mRNAs, to_be_inhibited, length, treshold_to_accept):
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
    mRNAsubsequences = [] #with this the function _find_all_unique_substrings_from_string is only called once instead of each time in the loop of _overlap_substring_on_mRNAs
    for mRNA in mRNAs:
        mRNAsubsequences.append(list(_find_all_unique_substrings_from_string(mRNA, length)))

    unique_mRNA_subsequences = set()
    unique_mRNA_subsequences.update([item for sublist in mRNAsubsequences for item in sublist])

    #substrings = _add_variability_to_substrings(substrings,treshold_to_accept)


    substring_results = []
    for substring in unique_mRNA_subsequences:
        miRNA_res = _overlap_substring_on_mRNAs(substring, mRNAsubsequences, treshold_to_accept)
        substring_results.append([substring, miRNA_res])
    logger.debug(substring_results)

def main():
    mrna_filepath = "term_genes/homosapiens_only=false,v2/product_mRNA.json"
    #gene_list = ["UniProtKB:Q16613", "UniProtKB:O15530", "UniProtKB:Q9Y243"]
    gene_list = util.get_identifier_values_from_json(mrna_filepath, "gene")[0]
    to_be_inhibited = [1, 1, 1]

    mRNAs = util.get_identifier_values_from_json(mrna_filepath, "mRNA")[0] #mRNAs = ["ABABABABABABABABABAB","ABCABCABABABABABABABABABABABCABC","ABCABCABCABCABCABC"]
    #predicted_miRNAs = predict_miRNAs(mRNAs, to_be_inhibited, 10, 0.5)


if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    main()
