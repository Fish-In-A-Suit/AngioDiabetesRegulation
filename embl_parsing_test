from Bio import SeqIO
import util

import logging
logger = logging.getLogger(__name__)

# This file showcases how to parse a file that stores multiple records in an embl file format.
def main():
    records = list(SeqIO.parse("src_data_files/miRNAdbs/mirbase_miRNA_06-02-2023.dat", "embl"))
    i = 0
    """
    for record in records:
        record_attribs = dir(record)
        logger.debug(f"record_attribs = {record_attribs}")
        logger.debug(f"Record id: {record.id}")
        logger.debug(f"     - name: {record.name}")
        logger.debug(f"     - description: {record.description}")
        logger.debug(f"     - seq: {record.seq}")
        logger.debug(f"     - dbxrefs: {record.dbxrefs}")
        logger.debug(f"     - features: {record.features}")
        logger.debug(f"     - annotations: {record.annotations}")
        #logger.debug(f"     - letter_annotations: {record.letter_annotations}")
        #logger.debug(f"     - lineage: {record.lineage}")
        i += 1
    logger.info(f"Read {i} records.")
    """

    util.save_mirbase_hsap_miRNAs("src_data_files/miRNAdbs/mirbase_miRNA_06-02-2023.dat")

    return

if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    main()