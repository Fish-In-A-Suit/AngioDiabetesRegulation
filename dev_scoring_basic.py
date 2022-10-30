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
    gene_set = set() # Set data type used instead of List, because a Set cannot have multiple occurences of the same element
    term_genes = [] # array of all genes across all terms; structure = {[term1, genes1], [term2, genes2], ... [term_n, genes_n]}
    for term in allowed_term_ids:
        genes = _import_genes_from_term_json(term, source_folder)
        term_genes.append([term, genes])
        gene_set.update(genes) # updates values in gene_set with genes, without duplicating existing elements
    logger.info(f"Found {len(gene_set)} different genes.")
    logger.debug(f"[term_genes]: {term_genes}")
    
    json_gene_scores=[]
    for gene in gene_set:
        #destination_file = f"{destination_folder}/{gene}.json"
        gene_count_result = _score_gene_basic(gene,term_genes)
        out = {"gene": gene, "count": gene_count_result[0], "terms:": gene_count_result[1]}
        json_gene_scores.append(out)
    logger.info(f"[gene_scores]: {json_gene_scores}")

    file = open(destination_file, "w+")
    file.write(json.dumps(json_gene_scores)+"\n")
    file.close()

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
    terms = ["GO:1903587", "GO:1903670"] # TODO: change to gene list from constants
    score_genes(terms, filename)

if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    main()
