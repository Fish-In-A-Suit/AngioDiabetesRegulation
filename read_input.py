import legacy.util as util
import legacy.constants as constants
import os

import logging
logger = logging.getLogger(__name__)


def main():
    out=util.read_input_file()
    logger.debug(out["processes"])
    terms = util.get_array_terms_from_input_list(out["GO_terms"])
    print(terms)
        

if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    main()