from asyncio.log import logger
import util
import constants
import json

import logging
logger = logging.getLogger(__name__)



def score_genes(allowed_term_ids, destination_file, source_folder="term_genes"):
    """
    Counts the number of appearances of all the genes across all specified json_files (which contain genes
    related to specific GO terms and are made by the find_genes_related_to_GO_terms function)
    """
    logger.info(f"Finding all genes from terms: {allowed_term_ids}")
    gene_set = set()
    term_genes = []
    for term in allowed_term_ids:
        genes = _import_genes_from_term_json(term, source_folder)
        term_genes.append([term, genes])
        gene_set.update(genes)
    logger.info(f"Found {len(gene_set)} different genes.")
    logger.debug(f"[term_genes]: {term_genes}")
    

    json_gene_scores=[]
    for gene in gene_set:
        #destination_file = f"{destination_folder}/{gene}.json"
        gene_count=_score_gene_basic(gene,term_genes)
        out = {"gene": gene, "count": gene_count}
        json_gene_scores.append(out)

    logger.info(f"[gene_scores]: {json_gene_scores}")

    file = open(destination_file, "w+")
    file.write(json.dumps(json_gene_scores)+"\n")
    file.close()

    logger.info (f"Done with scoring!")


def _score_gene_basic(gene, term_genes_list):
    """
    TODO
    """
    #perhaps we will have multiple scores/indicators
    #scores=[]
    count=0
    for term in term_genes_list:
        if gene in term[1]: count += 1
    logger.debug(f"Gene {gene} has {count} occurences.")
    return count
    
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
    filename = "gene_scores/test_score.json"
    terms = ["GO:1903587", "GO:1903670"] #change to gene list from constants
    score_genes(terms, filename)

if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    main()

