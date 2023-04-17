from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .Model import Product
import requests
import os
import time
from typing import Dict
import gzip
import logging

logger = logging.getLogger(__name__)

class miRDB60predictor:
    def __init__(self):
        # set the filepath to the miRDB prediction result file
        self._filepath = "src_data_files/miRNAdbs/miRDB_v6.0_prediction_result.txt.gz"
        # check if the file exists and download it if necessary
        self._check_file()
        # read the file into memory and decode the bytes to utf-8
        with gzip.open(self._filepath, "rb") as read_content:
            self._readlines = [line.decode("utf-8") for line in read_content.readlines()]
        # log the first 10 lines of the file
        logger.info(self._readlines[:10])
    
    def _check_file(self):
        _max_retries = 5
        _retry_delay = 2
        # create the directory where the file will be saved if it doesn't exist
        os.makedirs(os.path.dirname(self._filepath), exist_ok=True)
        if not os.path.exists(self._filepath):
            # download the file from the miRDB website with retries
            url = "https://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz"
            for i in range(_max_retries):
                try:
                    response = requests.get(url)
                    response.raise_for_status()
                except requests.exceptions.HTTPError as e:
                    logger.warning(f"Error downloading file from {url}: {e}")
                    if i < _max_retries - 1:
                        logger.warning(f"Retrying in {_retry_delay} seconds...")
                        time.sleep(_retry_delay)
                    else:
                        raise Exception(f"Failed to download file from {url} after {self._max_retries} attempts") from e

    def predict_from_product(self, product: Product, threshold: float = 0.0) -> Dict[str, float]:
        """
        Finds all miRNAs and their match strengths (hsa-miR-xxx, 72.2) from miRDB_readlines for mRNA_refseq (e.g. NM_xxxxx).

        :param product: Product object with refseq_nt_id attribute
        :param threshold: Minimum match strength to include in result_list
        :return: Dictionary containing miRNAs as keys and match strengths as values
        """
        if not product.refseq_nt_id:
            return None

        result_dict = {}

        # Iterate over each line in the list of read lines
        for line in self._readlines:
            # Check if product.refseq_nt_id is present in the line
            if product.refseq_nt_id.split(".")[0] in line:
                # Split the line by tabs to extract miRNA and match_strength
                miRNA, _, match_strength = line.strip().split("\t")
                # Convert match_strength to float
                match_strength = float(match_strength)
                # Add miRNA and match_strength to result_dict if match_strength >= threshold
                if match_strength >= threshold:
                    result_dict[miRNA] = match_strength

        return result_dict
