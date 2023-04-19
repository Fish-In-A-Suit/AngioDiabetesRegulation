import util as util
import constants as constants
import requests
import jellyfish
import os
from util import Timer

import logging
logger = logging.getLogger(__name__)

def get_mrna(gene_list, target_folder):
    """
    Loops through gene_list and constructs a .json file, where each element has structure:
    {
        "UniprotID":
        "EnsemblID":
        "mRNA":
    }
    File is saved as target_folder/product_mRNA.json
    Return: the whole json object
    """
    ensembl_ids = get_ensembl_ids_from_uniprot_ids(gene_list)
    mRNAs = _get_ensembl_mRNA_sequences(ensembl_ids)

    json_mRNAs = []
    for i in range(len(gene_list)):
        out = {"UniprotID" : gene_list[i], "EnsemblID" : ensembl_ids[i], "mRNA" : mRNAs[i]}
        json_mRNAs.append(out)
    util.save_json(json_mRNAs, os.path.join(target_folder, "product_mRNA.json"))
    return json_mRNAs

def get_ensembl_ids_from_uniprot_ids(gene_list):
    ensembl_ids=[]
    for gene in gene_list:
        ensembl_ids.append(_get_ensembl_id_from_uniprot_id(gene))
    return ensembl_ids

def _get_ensembl_id_from_uniprot_id(gene):
    gene_id = gene.split(":")[1]
    up_query = util._uniprot_query_API(gene_id, type="prot")
    ensembl_id = util._return_ensembl_from_id_and_uniprot_query(gene_id, up_query)
    return ensembl_id

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

def main():
    # manual settings
    # gene_list = ["UniProtKB:Q16613", "UniProtKB:O15530", "UniProtKB:Q9Y243"]
    # gene_list = util.get_uniprotids_from_json("gene_scores/test_score_homosapinesonly=false,v2-term_enums,cross_section_postprocess.json")[0]
    timer = Timer()
    source_file = os.path.join(constants.TARGET_FOLDER, "product_scores.json")
    gene_list = util.get_identifier_values_from_json(source_file, "product")[0]

    # Maybe get all genes from term_products.json into gene_list.
    # Might not be needed since we are only interested in top 10. Commented function call below executes the desired algorithm.
    # gene_list = util.get_uniprotids_from_json("gene_scores/test_score_homosapinesonly=false,v2-term_enums,cross_section.json")[0]

    get_mrna(gene_list, constants.TARGET_FOLDER)

    # appends mRNA RefSeq__NT accession IDs (from uniprotkb_human_idmapping.dat) to the result array, makes a new file called product_mRNA_NCBIacc
    util.product_mRNA_json_append_refseqIds()

    # create the mRNAs file to be processed by the c++ code (for miRNA annealing)
    util.save_mRNAs_for_cpp(constants.TARGET_FOLDER + "/product_mRNA_refseq.json")

    timer.print_elapsed_time()


if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    main()
