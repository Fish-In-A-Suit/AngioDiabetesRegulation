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
    #input = util.read_file_as_json(os.path.join(constants.TARGET_FOLDER, "product_mRNA_miRDB-predicted-miRNA-matching.json"))
    #for sequence_comparison_results.json:
    input = util.read_file_as_json(os.path.join(constants.TARGET_FOLDER, "sequence_comparison_results.json"))
    #make prediction dictionary:
    predict_dict = {}
    for mirna in input:
        for mirnaname in mirna:
            predict_dict[mirna[mirnaname][0]] =  mirna[mirnaname][1] #mirna[mirnaname][1] #TODO change the format of json!
    logger.debug(predict_dict["hsa-let-7a-1"])

    #get scores
    scores = util.read_file_as_json(os.path.join(constants.TARGET_FOLDER, "product_scores.json"))
    #make a scores dictionary:
    scores_dict = {}
    for product in scores:
        scores_dict[product["product"]] = product["al_corr_score"]
    logger.debug(scores_dict)

    mirna_scores = score_miRNAs(predict_dict, scores_dict, treshold=0.3)
    
    logger.info(f"[mirna_scores]: {mirna_scores}")

    mirna_scores = util.sort_list_of_dictionaries(mirna_scores, "basic_score")
    util.save_json(mirna_scores, os.path.join(constants.TARGET_FOLDER, "mirna_scores.json"))
    logger.info (f"Done with scoring!")

    #score_genes_v2(source_folder=constants.TARGET_FOLDER, use_cross_section=True)

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