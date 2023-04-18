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
   run_get_products = False
   run_score_products = False
   run_get_mRNA = False
   run_predict_miRNA_overlap = False # predict miRNAs based on brute-force, TODO
   run_predict_miRNA_miRDB = True # predict miRNAs from the miRDB database
   run_score_miRNA = False
   
   # inputs being read from file implementation
   # inputs = util.read_input_file()

   if run_get_products:
      terms = util.get_array_terms("ALL") # old implementation of the angio-diabetes problem
      # terms = util.get_array_terms_from_input_list(inputs["GO_terms"]) 
      import dev_go_download
      dev_go_download.main(terms)
      # dev_go_download.main(terms, debug_shorten_terms_count=20, debug=True)

   if run_score_products:
      import dev_scoring_products_basic
      # dev_scoring_products_basic.main()
      dev_scoring_products_basic.main(use_legacy=True)
   
   if run_get_mRNA:
      import dev_mrna_download
      dev_mrna_download.main()  

   if run_predict_miRNA_miRDB:
      import dev_predict_miRNA_v2
      dev_predict_miRNA_v2.main()

if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    main()





