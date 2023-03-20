import util
import constants
import os

import logging
logger = logging.getLogger(__name__)


def score_miRNAs(overlaps, scores_dict, treshold):
    logger.info(f"Started scoring miRNAs")
    mirna_scores = []
    for mirna in overlaps:
        basic_score = 0
        inhibited = []
        for product in overlaps[mirna]:
            product_score = scores_dict[product]
            overlap = overlaps[mirna][product]
            if overlap >= treshold: inhibited.append(product)
            basic_score += _score_mirna_basic(product_score, overlap, treshold)

        out = {"mirna": mirna, "basic_score": basic_score, "products_inhibited": inhibited}
        mirna_scores.append(out)

    return mirna_scores

def _score_mirna_basic(product_score, overlap, treshold):
    if overlap < treshold:
        a = -1
    elif overlap >= treshold:
        a = 1

    modifier = a * product_score
    return modifier

def main():
    """
    Idealy the prediction format should minimally be:
    [
        {
            "mirna_id": hsa-let-7a-1,
            "overlaps": {
                "UniProtKB:Q0VGL1": 0.2125,
                "UniProtKB:Q96GR2": 0.2125
            }

        }
    """
    #for product_mRNA_miRDB-predicted-miRNA-matching.json: TODO some post processing?
    input = util.read_file_as_json(os.path.join(constants.TARGET_FOLDER, "product_mRNA_miRDB-predicted-miRNA-matching.json"))
    predict_dict = {}
    for mrna in input:
        if len(mrna["RefSeq_NT_IDs"]) is not 0:
            for mirna in mrna["miRNA_matches"][mrna["RefSeq_NT_IDs"][0]]:
                if mirna[0] not in predict_dict:
                    predict_dict[mirna[0]] = {}
                predict_dict[mirna[0]][mrna["UniprotID"]] = mirna[1]
    #logger.debug(predict_dict)
    #for sequence_comparison_results.json:
    """input = util.read_file_as_json(os.path.join(constants.TARGET_FOLDER, "sequence_comparison_results.json"))
    #make prediction dictionary:
    predict_dict = {}
    for mirna in input:
        for mirnaname in mirna:
            predict_dict[mirna[mirnaname][0]] =  mirna[mirnaname][1] #mirna[mirnaname][1] #TODO change the format of json!
    logger.debug(predict_dict["hsa-let-7a-1"])"""

    #get scores
    scores = util.read_file_as_json(os.path.join(constants.TARGET_FOLDER, "product_scores.json"))
    #make a scores dictionary:
    scores_dict = {}
    for product in scores:
        scores_dict[product["product"]] = product["al_corr_score"]
    #logger.debug(scores_dict)

    mirna_scores = score_miRNAs(predict_dict, scores_dict, treshold=0.9)
    
    logger.info(f"[mirna_scores]: {mirna_scores}")

    mirna_scores = util.sort_list_of_dictionaries(mirna_scores, "basic_score")
    util.save_json(mirna_scores, os.path.join(constants.TARGET_FOLDER, "mirna_scores.json"))
    logger.info (f"Done with scoring!")
        

if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    main()