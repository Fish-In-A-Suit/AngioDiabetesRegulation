import util
import constants
import requests

import logging
logger = logging.getLogger(__name__)

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

def _get_ensembl_sequence_API(id):
    """
    This function queries ensembl for nucleotide sequence
    Input of ensembl ID's must be a string
    """
    logger.debug(
        f"Starting Ensembl API for id {id}")
    response = requests.get(f"https://rest.ensembl.org/sequence/id/{id}?object_type=transcript;type=cdna", headers={"Content-Type": "text/plain", }
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
    
def predict_miRNAs(mRNAs, to_be_inhibited, length, treshold_to_accept):
    """
    Aljosa&Ladi algorithm for finding miRNA of certain length from a list of mRNA

    treshold_to_accept - what is the minimal percentage of the miRNA length that must fit to a mRNA in order to be considered as 'fitting'
    """



def main():
    gene_list = ["UniProtKB:Q16613", "UniProtKB:O15530", "UniProtKB:Q9Y243"]
    to_be_inhibited = [1, 1, 1]
    ensembl_ids = get_ensembl_ids_from_uniprot_ids(gene_list)
    mRNAs = _get_ensembl_mRNA_sequences(ensembl_ids)
    #predicted_miRNAs = predict_miRNAs(mRNAs, to_be_inhibited, 10, 0.5)


if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    main()
