import requests
import json
import os
from tqdm import tqdm
import logging
logger = logging.getLogger(__name__)

class go_analysis_model:
    def __init__(self, target_process=[], go_terms_definitions=[], target_folder="./", homosapiensonly=True):
        self.target_folder = target_folder
        self.homosapiensonly = homosapiensonly
        self.target_process = target_process #list of dicts
        self.go_terms_definitions = go_terms_definitions #list of dicts
        #flags
        self.products_download_updated = False
        self.product_scores_updated = False
        self.mrna_download_updated = False
        self.mrna_overlaps_updated = False
        self.mrna_scores_updated = False
        self.report_updated = False
        logger.info(f"Initialized model")

    def check_go_dict_format(self, go_dict):
        if not all(k in go_dict for k in ("GO_id","process","direction","weight")):
            raise KeyError(f"GO_dict: {go_dict} does not include all the required elements")
        if not ("GO:" in go_dict["GO_id"]):
            raise ValueError(f"missing 'GO:' in GO term id {go_dict['GO_id']}")
        #TODO maybe check weight to be 0-1 or data types?

    def add_goterm(self, entries:list|dict|tuple):
        if any(isinstance(el, list | dict | tuple) for el in entries):
            logger.info(f"Detected multiple entries input")
            for sub_entry in entries:
                self.__add_single_goterm(sub_entry)
        else:
            logger.info(f"Detected single entry input")
            self.__add_single_goterm(entries)

    def __add_single_goterm(self, new_entry:list|dict|tuple):
        """
        Can be 1d-list, tuple or dict:
        List of format:
        [GO_id, process, direction, weight]
        Tuple of format:
        (GO_id, process, direction, weight)
        Dict of format:
        {"GO_id":"GO:1903589", "process":"angio", "direction":"+", "weight":1}
        """

        #Extract data from input enqiry. Try as dict and as list or tuple. Also check if all requered fields are present
        try:
            go_id = new_entry["GO_id"]
            process = new_entry["process"]
            direction = new_entry["direction"]
            weight = new_entry["weight"]            
        except TypeError:
            try:
                go_id = new_entry[0]
                process = new_entry[1]
                direction = new_entry[2]
                weight = new_entry[3]
            except IndexError:
                raise IndexError(f"input entry: {new_entry} does not include all the required elements")
        except KeyError:
            raise KeyError(f"input entry: {new_entry} does not include all the required elements")
        #Check if GO_id is already in go_terms_definitions
        if not any(d['GO_id'] == go_id for d in self.go_terms_definitions):
            raise ValueError(f"{go_id} already in go_terms_definitions list.")

        dict_to_be_appended = {"GO_id":go_id, "process":process, "direction":direction, "weight":weight}
        #if not isinstance(self.check_go_dict_format(dict_to_be_appended), Exception):
        self.check_go_dict_format(dict_to_be_appended)

        self.go_terms_definitions.append(dict_to_be_appended)

        logger.info(f"Appended {dict_to_be_appended} to GO term list")
        #change flags
        self.products_download_updated = False
        self.product_scores_updated = False
        self.mrna_download_updated = False
        self.mrna_overlaps_updated = False
        self.mrna_scores_updated = False
        self.report_updated = False

    def load_data(self, which="all"):
        """
        Which data to load?
        config
        go_terms
        productlist
        productscores
        mRNA
        miRNAoverlap
        miRNAscores
        report
        """


    def save_data(self, which="all"):
        """
        Which data to save?
        config
        go_terms
        productlist
        productscores
        mRNA
        miRNAoverlap
        miRNAscores
        report
        """


        
        