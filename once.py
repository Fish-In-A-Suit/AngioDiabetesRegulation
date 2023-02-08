import util

import logging
logger = logging.getLogger(__name__)

# This file contains the code for any functions that you would like to run only once.

def main():
    util.save_mRNAs_for_cpp("test_run_2/product_mRNA_refseq.json")

if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    main()