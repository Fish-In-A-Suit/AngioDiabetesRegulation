import requests
import json
import os
from tqdm import tqdm
import logging
logger = logging.getLogger(__name__)

class go_analysis_model:
    def __init__(self, config_filepath = None, target_processes=[], go_terms_definitions=[], target_folder="./", homosapiensonly=True):
        self.target_folder = target_folder
        self.homosapiensonly = homosapiensonly
        self.target_processes = target_processes #list of dicts
        self.go_terms_definitions = go_terms_definitions #list of dicts
        self.config_filepath = os.path.join(target_folder, "config.txt") if config_filepath is None else config_filepath
        #flags - maybe better use hash?
        #self.products_download_updated = False
        #self.product_scores_updated = False
        #self.mrna_download_updated = False
        #self.mrna_overlaps_updated = False
        #self.mrna_scores_updated = False
        #self.report_updated = False
        logger.info(f"Initialized model")

    def check_go_dict_format(self, go_dict):
        if not all(k in go_dict for k in ("GO_id","process","direction","weight","description")):
            raise KeyError(f"GO_dict: {go_dict} does not include all the required elements")
        if not ("GO:" in go_dict["GO_id"]):
            raise ValueError(f"missing 'GO:' in GO term id {go_dict['GO_id']}")
        #TODO maybe check weight to be 0-1 or data types?

    def add_goterm(self, entries:list|dict|tuple):
        def __add_single_goterm(new_entry:list|dict|tuple):
            """
            Can be 1d-list, tuple or dict:
            List of format:
            [GO_id, process, direction, weight]
            Tuple of format:
            (GO_id, process, direction, weight)
            Dict of format:
            {"GO_id":"GO:1903589", "process":"angio", "direction":"+", "weight":1}
            Optional 5th element is description.
            """

            #Extract data from input enqiry. Try as dict and as list or tuple. Also check if all requered fields are present
            try:
                go_id = new_entry["GO_id"]
                process = new_entry["process"]
                direction = new_entry["direction"]
                weight = new_entry["weight"]
                desctiption = new_entry["description"] if "description" in new_entry.keys() else ""         
            except TypeError:
                try:
                    go_id = new_entry[0]
                    process = new_entry[1]
                    direction = new_entry[2]
                    weight = new_entry[3]
                    desctiption = new_entry[4] if len(new_entry) == 5 else ""
                except IndexError:
                    raise IndexError(f"input entry: {new_entry} does not include all the required elements")
            except KeyError:
                raise KeyError(f"input entry: {new_entry} does not include all the required elements")
            
            #Check if GO_id is already in go_terms_definitions
            if any(d['GO_id'] == go_id for d in self.go_terms_definitions):
                raise ValueError(f"{go_id} already in go_terms_definitions list.")

            dict_to_be_appended = {"GO_id":go_id, "process":process, "direction":direction, "weight":weight, "description":desctiption}
            #if not isinstance(self.check_go_dict_format(dict_to_be_appended), Exception):
            self.check_go_dict_format(dict_to_be_appended)

            self.go_terms_definitions.append(dict_to_be_appended)
            logger.info(f"Appended {dict_to_be_appended} to GO term list")
            self.save_data("config")
        
        if any(isinstance(el, list | dict | tuple) for el in entries):
            logger.info(f"Detected multiple entries input")
            for sub_entry in entries:
                __add_single_goterm(sub_entry)
        else:
            logger.info(f"Detected single entry input")
            __add_single_goterm(entries)

    def add_target_process(self, entries:list|dict|tuple):
        def __add_single_target_process(new_entry):
            """
            Can be 1d-list, tuple or dict:
            List of format:
            [process, direction]
            Tuple of format:
            (process, direction)
            Dict of format:
            {"process":"angio", "direction":"+"}
            """
            try:
                process = new_entry["process"]
                direction = new_entry["direction"]         
            except TypeError:
                try:
                    process = new_entry[0]
                    direction = new_entry[1]
                except IndexError:
                    raise IndexError(f"input entry: {new_entry} does not include all the required elements")
            except KeyError:
                raise KeyError(f"input entry: {new_entry} does not include all the required elements")
            
            #Check if the process is already in target_process list
            if any(d['process'] == process for d in self.target_processes):
                raise ValueError(f"{process} already in target_processes list.")

            dict_to_be_appended = {"process": process, "direction":direction}
            self.target_processes.append(dict_to_be_appended)
            self.save_data("config")


        if any(isinstance(el, list | dict | tuple) for el in entries):
            logger.info(f"Detected multiple entries input")
            for sub_entry in entries:
                __add_single_target_process(sub_entry)
        else:
            logger.info(f"Detected single entry input")
            __add_single_target_process(entries)

    def load_data(self, which="all"):
        """
        Which data to load?
        config - includes settings, desired processes and go_terms
        productlist
        productscores
        mRNA
        miRNAoverlap
        miRNAscores
        report
        """
        if which in ("config","all"):
            #This data is tab separated.
            #It is always 
            line_element_delimiter='\t'
            """
            returns dictionary of the contents of the input file
            {
            setting_name:value,
            processes:[{"name":"process_name1", "direction":"+"}],
            GO_terms: [{"GO_id":"GO:1903589", "process":"angio", "direction":"+", "weight":1}]
            }
            """
            def process_comment(line, comment_delimiter="#", line_keep_delimiter="###"):
                """
                Processes a comment in the line: returns the part of the line before the comment.
                Example line: 'attg 4 binding_strength # this is an example' -> returns: 'attg 4 binding_strength '
                
                Parameters:
                - line: the line whose comment to process
                - comment_delimiter: the character used to denote a comment
                - line_keep_delimiter: a special set of characters to denote a "logic line", the result of which you should keep
                                        logic lines should be marked with "###" at the start. For a logic line, the returned result is
                                        line without the line_keep_delimiter
                """
                if line_keep_delimiter in line:
                    return line.replace(line_keep_delimiter, "")
                
                if comment_delimiter in line:
                    return line.split(comment_delimiter)[0]
                else:
                    return line
                
            with open(self.config_filepath, "r") as read_content:
                read_lines = read_content.read().splitlines()[2:] #skip first 2 lines
                section = "" #what is the current section i am reading
                self.go_terms_definitions = [] #initialise the variables
                self.target_processes = [] #initialise the variables
                for line in read_lines:
                    line = process_comment(line)
                    if "settings" in line:
                        section = "settings"
                        continue
                    elif "processes" in line:
                        section = "process"
                        continue
                    elif "GO_terms" in line:
                        section = "GO"
                        continue
                    if section == "settings":
                        chunks = line.split(line_element_delimiter)
                        if chunks[0] == "target_folder":
                            self.target_folder = chunks[1]
                        elif chunks[1] == "homosapiensonly":
                            self.homosapiensonly = chunks[1]
                    elif section == "process":
                        chunks = line.split(line_element_delimiter)
                        self.target_processes.append({"process":chunks[0], "direction":chunks[1]})
                    elif section == "GO":
                        chunks = line.split(line_element_delimiter)
                        if len(chunks) == 5:
                            self.go_terms_definitions.append({"GO_id":chunks[0], "process":chunks[1], "direction":chunks[2], "weight":chunks[3], "description": chunks[4]})
                        else:
                            self.go_terms_definitions.append({"GO_id":chunks[0], "process":chunks[1], "direction":chunks[2], "weight":chunks[3], "description": ""})
            logger.info(f"Succesfully loaded config from {self.config_filepath}")

    def generate_hash(self, object):
        """
        do we need it?
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
        if which in ("config", "all"):
            #append setings that need to be saved - declaration is manual for now.
            settings = {}
            settings["target_folder"] = self.target_folder
            settings["homosapiensonly"] = self.homosapiensonly
            #create directory path, ignore if already exists
            os.makedirs(os.path.dirname(self.config_filepath), exist_ok=True)          
            # Writing to file
            with open(self.config_filepath, "w") as file:
                #First lets write the header:
                file.write("# --- Config file for Gene Ontoloy Reverse Lookup --- #" + "\n")
                file.write("###[name] is the start of new section and includes description of data, data is tab delimited" + "\n")
                #Then goes the ###seting header and then the settings data
                file.write("###settings"+"\n")
                for key in settings.keys():
                    file.write(f"{key}\t{settings[key]}\n")
                #Then goes the target processes declaration
                file.write("###processes [proces name] [to be expressed + or suppressed -]"+"\n")
                for process in self.target_processes:
                    file.write(f"{process['process']}\t{process['direction']}\n")
                #finally go the go terms declarations
                file.write("###GO_terms [GO id] [process] [upregulated + or downregulated - or general 0] [weight 0-1]"+"\n")
                for go_term in self.go_terms_definitions:
                    file.write(f"{go_term['GO_id']}\t{go_term['process']}\t{go_term['direction']}\t{go_term['weight']}\t{go_term['description']}\n")
            logger.info(f"Succesfuly saved the config to {self.config_filepath}")


        
        