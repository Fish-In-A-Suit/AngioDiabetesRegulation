# Main script file

"""
Introduction to Gene Ontology:
    - GO Annotation = a statement about the function of a particular gene (created by associating
      a gene / gene product with a GO Term). It includes statements about
         (1) how the gene functions at the molecular level
         (2) where in the cell it functions and
         (3) what biological pathways it helps carry out.
      At the minimum, an annotation consists of:
         (1) Gene product (may be a protein, RNA, etc)
         (2) GO term
         (3) Reference
         (4) Evidence (go evidence codes: http://geneontology.org/docs/guide-go-evidence-codes/)
      Additionally, annotations can be equipped with 'annotation qualifiers': NOT, contributes_to, colocalizes_with
"""

     
import logging
logger = logging.getLogger(__name__)

import util


def main():
   run_get_products = True
   run_score_products = False
   run_get_mRNA = False
   run_predict_miRNA_overlap = False
   run_score_miRNA = False
   
   inputs = util.read_input_file()

   if run_get_products:
      terms = util.get_array_terms_from_input_list(inputs["GO_terms"])
      import dev_go_download
      dev_go_download.main(terms)

   if run_score_products:
      import dev_scoring_products_basic
      dev_scoring_products_basic.main()

   

if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    main()





