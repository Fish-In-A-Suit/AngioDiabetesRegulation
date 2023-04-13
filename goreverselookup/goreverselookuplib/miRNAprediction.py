from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .Model import Product
import requests
import os
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
        # create the directory where the file will be saved if it doesn't exist
        os.makedirs(os.path.dirname(self._filepath), exist_ok=True)
        if not os.path.exists(self._filepath):
            # download the file from the miRDB website
            url = "https://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz"
            response = requests.get(url)
            # save the file to the specified filepath
            with open(self._filepath, "wb") as f:
                f.write(response.content)
            # log a message indicating the file has been downloaded
            logger.info(f"Downloaded miRDB_v6.0_prediction_result.txt.gz to {self._filepath}")

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
