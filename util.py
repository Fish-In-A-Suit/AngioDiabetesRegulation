# Utility functions file
import requests
import constants
import json
import os
import collections
import shutil
import time
from Bio import SeqIO
import itertools

import logging
logger = logging.getLogger(__name__)

import sys

_response_cycle_counter = 0
_uniprot_query_result = ""

_uniprot_genes_readlines = ""
_zfin_ortholog_readlines = ""
_xenbase_ortholog_readlines = ""
_mgi_ortholog_readlines = ""
_rgd_ortholog_readlines = ""

def read_input_file(filepath=os.path.join(constants.TARGET_FOLDER, "input.txt"), line_element_delimiter='\t'):
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
        
    output = {}
    temp_processes = []
    temp_GO = []
    with open(filepath, "r") as read_content:
        read_lines = read_content.read().splitlines()[2:] #skip first 2 lines
        section = "" #what is the current section i am reading
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
                output[chunks[0]] = chunks[1]
            elif section == "process":
                chunks = line.split(line_element_delimiter)
                if "processes" in output.keys():
                    temp_processes = output["processes"]
                temp_processes.append({"name":chunks[0], "direction":chunks[1]})
                output["processes"] = temp_processes
            elif section == "GO":
                chunks = line.split(line_element_delimiter)
                if "GO_terms" in output.keys():
                    temp_GO = output["GO_terms"]
                temp_GO.append({"GO_id":chunks[0], "process":chunks[1], "direction":chunks[2], "weight":chunks[3]})
                output["GO_terms"] = temp_GO
    
    return output


def get_array_terms_from_input_list(input_list, filter=[]):
    output = [element["GO_id"] for element in input_list if all(x in element.values() for x in filter)]
    return output


def get_array_terms(array_name, term_shrink=True):
    """
    Returns the specified array terms, possible options:
      - 'ALL': all arrays combined
      - 'ANGIOGENESIS': all angiogenesis arrays combined
      - 'ANGIOGENESIS-POSITIVE'
      - 'ANGIOGENESIS-NEGATIVE'
      - 'ANGIOGENESIS-GENERAL'
      - 'DIABETES': all diabetes arrays combined
      - 'DIABETES-POSITIVE'
      - 'DIABETES-NEGATIVE'
      - 'DIABETES-GENERAL'
    
    If term_shrink=True, the final list will only consist of GO:xxxxx, if false, returned list will consist of term names AND GO:xxxxx
    """
    def _shrink_term_list(list):
        i = 0
        result_list = []
        for element in list:
            if i%2: # 0 = generic term name, 1 = GO:xxxx, 2 = generic term name, 3 = GO:xxxx -> modulo op is 1 at 1, 3, 5 etc
                result_list.append(element)
            i = i+1
        return result_list

    if array_name == 'ALL':
        if term_shrink: return _shrink_term_list(constants.TERMS_ANGIOGENESIS_GENERAL + constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY + constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY + constants.TERMS_DIABETES_GENERAL + constants.TERMS_DIABETES_NEGATIVE_ARRAY + constants.TERMS_DIABETES_POSITIVE_ARRAY)
        else: return constants.TERMS_ANGIOGENESIS_GENERAL + constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY + constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY + constants.TERMS_DIABETES_GENERAL + constants.TERMS_DIABETES_NEGATIVE_ARRAY + constants.TERMS_DIABETES_POSITIVE_ARRAY
    elif array_name == 'ANGIOGENESIS':
        if term_shrink: return _shrink_term_list(constants.TERMS_ANGIOGENESIS_GENERAL + constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY + constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY)
        else: return constants.TERMS_ANGIOGENESIS_GENERAL + constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY + constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY
    elif array_name == 'ANGIOGENESIS-NEGATIVE':
        if term_shrink: return _shrink_term_list(constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY)
        else: return constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY
    elif array_name == 'ANGIOGENESIS-POSITIVE':
        if term_shrink: return _shrink_term_list(constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY)
        else: return constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY
    elif array_name == 'ANGIOGENESIS-GENERAL':
        if term_shrink: return _shrink_term_list(constants.TERMS_ANGIOGENESIS_GENERAL)
        else: return constants.TERMS_ANGIOGENESIS_GENERAL
    elif array_name == 'DIABETES':
        if term_shrink: return _shrink_term_list(constants.TERMS_DIABETES_GENERAL + constants.TERMS_DIABETES_NEGATIVE_ARRAY + constants.TERMS_DIABETES_POSITIVE_ARRAY)
        else: return constants.TERMS_DIABETES_GENERAL + constants.TERMS_DIABETES_NEGATIVE_ARRAY + constants.TERMS_DIABETES_POSITIVE_ARRAY
    elif array_name == 'DIABETES-NEGATIVE':
        if term_shrink: return _shrink_term_list(constants.TERMS_DIABETES_NEGATIVE_ARRAY)
        else: return constants.TERMS_DIABETES_NEGATIVE_ARRAY
    elif array_name == 'DIABETES-POSITIVE':
        if term_shrink: return _shrink_term_list(constants.TERMS_DIABETES_POSITIVE_ARRAY)
        else: return constants.TERMS_DIABETES_POSITIVE_ARRAY
    elif array_name == 'DIABETES-GENERAL':
        if term_shrink: return _shrink_term_list(constants.TERMS_DIABETES_GENERAL)
        else: return constants.TERMS_DIABETES_GENERAL
    else:
        print(array_name + " could not be found! Returning empty array.")
        empty = []
        return empty

def get_files_in_dir(dir_filepath, searchstring = ""):
    """
    Searches for files in directory. Use searchstring parameter to append all files that have a specific string in their name.
    Returns a single or a list of filepaths in the form dir_filepath/file
    """
    if not check_file_exists(dir_filepath):
        return []
    if searchstring == "":
        return os.listdir(dir_filepath)
    else: # search for searchstring, if file includes searchstring, append filepath
        result_filepaths = []
        for f in os.listdir(dir_filepath):
            if searchstring in f:
                result_filepaths.append(f"{dir_filepath}/{f}")
        return result_filepaths

def get_last_file_in_list(list):
    """
    Organizes the file by timestamps and returns the most recent file
    """
    timestamps = {}
    index = 0
    for f in list:
        timestamp = f.split("/")[-1].split("_")[1]
        # timestamp = timestamp.split(".")[0] no longer needed
        timestamps[index] = timestamp
        index += 1
    ordered_dict = collections.OrderedDict(timestamps)
    logger.debug(f"ordered_dict = {ordered_dict}")
    if len(ordered_dict) > 0: return list[len(ordered_dict)-1]
    else: return "empty" # it's okay, because file "empty" doesn't exist

def read_file_as_json(filepath):
    """
    Reads the file into json
    """
    with open(filepath, "r") as read_content:
        return json.load(read_content)

def readlines(filepath):
    """
    Reads the lines of the specified filepath
    """
    with open(filepath, "r") as read_content:
        return read_content.readlines()

def writelines(filepath, lines):
    """
    Writes 'lines' to filepath file
    """
    file = open(filepath, "w")
    file.writelines(lines)
    file.close()

def append_to_file(to_append, filepath, append_if_exists=False, add_linebreaks=True):
    """
    Appends to_append to filepath.
      - append_if_exists: If True, appends to_append regardless if the same entry already exists. Otherwise appends only if to_append is a unique entry
      - add_linebreaks: If True, adds \n add the end of to_append (if \n isn't already inside to_append)
    """
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "a+") as f:
        logger.debug(f"Opened file {filepath}")
        if "\n" not in to_append and add_linebreaks == True:
            to_append = f"{to_append}\n"
        if append_if_exists == True:
            # appends to_append regardless if the same entry already exists
            f.write(to_append)
            f.close()
        else:
            # appends only if there is not a similar entry
            f.seek(0) # a+ has file read pointer at bottom, f.seek(0) returns to top
            if to_append not in f.read():
                f.seek(0,2) # move file pointer back to the end of the file
                f.write(to_append)
                f.close()
                logger.debug(f"Written {to_append} to file {filepath}.")

def filepath_striplast(filepath):
    """
    Steps back one step in the file tree. Example: dirA/dirB/fileC --> dirA/dirB/
    """
    result_filepath = ""
    if "/" in filepath:
        split = filepath.split("/")
        for i in range(len(split)-1):
            result_filepath = os.path.join(result_filepath, split[i])
    logger.debug(f"result_filepath = {result_filepath}")
    return result_filepath

def load_list_from_file(filepath, array, no_elements_in_line=1, break_character=" "):
    """
    Reads a file line by line. If a line contains multiple elements, adjust appropriate no_elements_in_line and break character
    """
    try:
        file = open(filepath, "r")
        lines = file.readlines()
        for line in lines:
            if "\n" in line:
                line = line.replace("\n","")
            if no_elements_in_line > 1:
                splitlist = line.split(break_character)
                for e in splitlist:
                    array.append(e)
            else:
                array.append(line)
    except OSError:
        logger.debug(f"{filepath} does not exist!")

def append_to_filepath(srcfilepath, append):
    """
    Appends the 'append' to the file before the filetype.
    Example usage: util.append_to_file("term_genes/GO-0001525.json", ";params=homosapiens_only,v1")
    --> result: term_genes/GO-0001525;params=homosapiens_only,v1.json
    """
    split=srcfilepath.split(".")
    dstfilepath = split[0] + append + "." + split[1]
    os.rename(srcfilepath, dstfilepath)

def term_sort_into_file(src_folder, dest_folder, term_list):
    """
    Reads src_file, compare dir, copy existing files to dest_file
    """
    files_copied = []
    files_not_found = []
    for term in term_list:
        if ":" in term: term = term.replace(":","-")
        src_term_filepath = f"{src_folder}/{term}.json"
        if os.path.exists(src_term_filepath):
            files_copied.append(src_term_filepath)
            shutil.copy(src_term_filepath, dest_folder)
        else:
            files_not_found.append(src_folder)
    logger.info(f"Finsihed file copy. Files copies = {files_copied}, files not found = {files_not_found}")
    
def save_json(json_object, filepath):
    """
    Stores json_object to file at filepath
    """
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as f:
        json.dump(json_object, f, indent=4)
    logger.info(f"Stored json to file {filepath}")

def save_txt(string, filepath):
    """
    Stores string to file at filepath
    """
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as f:
        f.write(string)
    logger.info(f"Stored string to file {filepath}")

def load_json_by_terms(src_folder, terms):
    """
    Returns the json files from src_folder that correspond to term names
    """
    result_file_list = []
    missing_files = []
    for term in terms:
        _fi = ""
        if ":" in term: _fi = term.replace(":","-") # switch from GO:xxxx to GO-xxxx
        else: _fi = term 
        if ".json" not in _fi: _fi = f"{_fi}.json" # append .json
        if os.path.exists(f"{src_folder}/{_fi}"): result_file_list.append(f"{src_folder}/{_fi}")
        else: missing_files.append(f"{src_folder}/{_fi}")
    if len(missing_files) > 0:
        logger.info(f"Missing files in {src_folder}: {missing_files}")
        input("Missing files may be caused by a term only including subterm-genes and not any term specific genes. Press any key to continue")
    logger.info(f"Returning {len(result_file_list)} result files;  -> {result_file_list}")
    return result_file_list

def load_sequence_comparison_results(filepath = "src_data_files/sequence_comparison_results.txt"):
    if os.path.isfile(filepath):
        logger.info(f"Loading sequence comparison results from {filepath}")
        file_lines = []
        load_list_from_file(filepath, file_lines) # read all file lines from filepath into file_lines
        
        sequence_comparison_results_json = [] # master list, containing all current_miRNA_element_dict elements
        current_miRNA_element_dict = {} # a dictionary for a single miRNA element
        current_miRNA_mirbase_id = "" # eg. MI0000060
        current_miRNA_othername = "" # eg. hsa-let-7a-1
        current_miRNA_mRNA_matches_dict = {} 
        i = 0
        for line in file_lines:
            if "#" in line or "//" in line:
                continue
            if "\t" not in line:
                if i != 0:
                    # not the first iteration -> create & save current_miRNA_element_dict to sequence_comparison_results_json
                    current_miRNA_element_list = [current_miRNA_othername, current_miRNA_mRNA_matches_dict]
                    current_miRNA_element_dict[current_miRNA_mirbase_id] = current_miRNA_element_list
                    sequence_comparison_results_json.append(current_miRNA_element_dict)
                    # reset
                    current_miRNA_element_dict = {}
                    current_miRNA_mRNA_matches_dict = {}

                current_miRNA_identifiers_split = line.strip().split(",")
                current_miRNA_othername = current_miRNA_identifiers_split[1]
                current_miRNA_mirbase_id = current_miRNA_identifiers_split[0]
            if "\t" in line:
                line = line.replace("\t", "")
                mRNAid = line.split(": ")[0]
                match_strength = float(line.split(": ")[1])
                current_miRNA_mRNA_matches_dict[mRNAid] = match_strength
            i += 1
        return sequence_comparison_results_json
    else:
        logger.info(f"ERROR! Filepath {filepath} does not exist!")




def store_json_dictionaries(filepath, dictionaries):
    """
    Writes the json dictionaries to file at filepath
    """
    if json.dumps(dictionaries) != "[]":
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        file = open(filepath, "w+")
        json.dump(dictionaries, file, indent=4)
        file.close()
    else: logger.info("JSON for analysis progress not stored, as it is empty.")

def json_compare(file1, file2):
    """
    Compares the contents of file1 and file2 json files. Used during development to check for json
    file similarity/differences.

    Returns: True if no differences, False if differences
    """
    file1_json = read_file_as_json(file1)
    file2_json = read_file_as_json(file2)
    if file1_json == file2_json:
        logger.debug(f"Compare {file1} and {file2}: True")
        return True
    else:
        logger.debug(f"Compare {file1} and {file2}: False")
        return False

def cleanup_crash_term_jsons():
    """
    Deletes all except the last jsons for each term in term_genes_crash
    """
    logging.info("Cleaning up term_genes_crash.")
    jsons_dict = {} # dictionary in form of <string> term : <list> filepath snapshots
    all_files = os.listdir("term_genes_crash")
    logger.debug(f"all_files = {all_files}")
    previous_term = ""
    first_iteration = True
    term_snapshots = [] # list of all crash files belonging to specific term at different timestamps
    iteration = 0
    for f in all_files:
        current_term = f.split("_")[0]
        logger.debug(f"iteration: {iteration}, curterm = {current_term}, prevterm = {previous_term}")
        if current_term != previous_term and first_iteration!=True: # new term
            previous_term = current_term
            # if len(term_snapshots)>=1: jsons_dict[current_term] = term_snapshots
            jsons_dict[current_term] = term_snapshots
            term_snapshots = [] # reset  
        elif first_iteration == True: 
            term_snapshots.append(f)
            previous_term = current_term
        else: # previous_term same as current_term
            term_snapshots.append(f)
        
        if iteration == (len(all_files)-1): # final iteration, also append; MUST BE KEPT OUT OF THE IF CASCADE
            jsons_dict[current_term] = term_snapshots
        
        logger.debug(f"iteration: {iteration}, term_snapshots = {term_snapshots}, jsons_dict = {jsons_dict}")

        first_iteration = False
        iteration += 1
    logging.debug(f"jsons_dict = {jsons_dict}")

    for term,snapshots in jsons_dict.items():
        # extract timestamps
        snapshots_timestamps = []
        for filename in snapshots:
            snapshots_timestamps.append(get_timestamp_from_crashterm(filename))
        # order timestamps
        snapshots_timestamps.sort()
        # delete all but the last timestamp files
        snapshots_timestamps.pop() # removes last element
        for filename in snapshots:
            file_timestamp = get_timestamp_from_crashterm(filename)
            if file_timestamp in snapshots_timestamps and term in filename:
                os.remove(f"term_genes_crash\\{filename}")
                logging.info(f"Removed {filename}")

def merge_similar_term_products_in_json(crash_json):
    """
    Used in crash recovery. The new workflow uses temp_dictionary, which may append the same term and its products (eg GO:xxxx) to json_dictionaries
    but the very same term can already exist is json_dictionaries from crash recovery. Therefore, it prevents term duplications in the saved json.
    """
    previous_term = ""
    previous_term_products = []
    
    for index,element in enumerate(crash_json):
        current_term = element["GO_term"]
        current_term_products = element["products"]

        if current_term == previous_term: # merge term products
            merged_products = current_term_products + previous_term_products
            # pop TODO: POP BOTH!!! previous too
            crash_json.pop(index)
            # insert one index back
            logger.debug(f"index = {index}, crash_json = {crash_json}")
            crash_json[index-1]["GO_term"] = current_term
            crash_json[index-1]["products"] = merged_products
        
        previous_term = current_term
        previous_term_products = current_term_products
    return crash_json

def get_timestamp_from_crashterm(crashfile_string):
    """
    Extracts the first part of the timestamp. Eg GO-0001525_1668012531.381136_.json -> 1668012531
    """
    timestamp = crashfile_string.split("_")[1]
    timestamp_final = int(timestamp.split(".")[0])
    return timestamp_final

def get_last_geneId_in_crash_json(json):
    """
    Gets the Id of the last gene in supplied json. Used in the crash recovery algorithm
    """
    geneId = json[len(json)-1]["product"] # get last item in array (len(json[0])-1) and then "product", which is geneId
    logger.debug(f"geneId = {geneId}")
    return geneId
    
def get_last_GO_term_in_crash_json(json):
    """
    Gets the Id of the last gene in supplied json. Used in the crash recovery algorithm
    """
    GO_term = json[len(json)-1]["GO_term"] # get last item in array (len(json[0])-1) and then "GO_term", which is geneId
    logger.debug(f"GO_term = {GO_term}")
    return GO_term

def list_directionshrink(list, reference_element, forward=True):
    """
    Returns list elements either up to or from the reference_element.
        If forward = True, keeps all elements in list starting after the reference_element
        If forward = False, keeps all elements in list starting before the reference_element

    Example: 
    ls = [0,1,2,3,4,5,6]
    list_directionshrink(ls, 2, forward=True) --> [3,4,5,6]
    """
    result_list = []
    start_append_forward = False
    stop_append_backward = False
    
    if reference_element not in list:
        return list

    for element in list:
        if element == reference_element:
            start_append_forward = True
            stop_append_backward = True
        if start_append_forward == True and forward == True:
            result_list.append(element)
        if stop_append_backward == False and forward == False:
            result_list.append(element)
    result_list.remove(reference_element)
    logger.debug(f"[list_directionshrink]: Shrank from {len(list)} to {len(result_list)} elements.")
    return result_list

def sort_list_of_dictionaries(input, field, direction_reversed = True):
    """Sorts the list of dictionaries by the key "field", default direction is reversed (descending)"""
    return sorted(input, key=lambda d: d[field], reverse=direction_reversed)

def get_dict_key_at_index(dictionary, index):
    """
    Returns the value of the dictionary key at specified index.
    """
    keys_list = list(dictionary.keys()) # call list(dict) on a dictionary to return a list of its keys
    key_at_index = keys_list[index]
    return key_at_index

def get_dict_key_at_value(dictionary, in_value):
    """
    Returns the key of the specified value in the dictionary.
    """
    for key,value in dictionary.items():
        if value == in_value:
            return key
    logger.debug(f"ERROR! Couldnt find value {in_value} during dictionary search. Returning None.")
    return None

def remove_dict_list_elements(src_dict, string_to_remove):
    """
    Removes the elements which contain 'string_to_remove' from the lists inside a dictionary. The dictionary has to have the form:
    key1: list1, key2: list2, key3: list3; example 'UniProtKB:Q0VGL1': ['NM_001008395.3', 'XM_017012198.1'], 'UniProtKB:Q96GR2': ['NM_001199377.1', 'NM_015162.4'];
    if you call this function using "XM" as 'string_to_remove', it produces: 'UniProtKB:Q0VGL1': ['NM_001008395.3'], 'UniProtKB:Q96GR2': ['NM_001199377.1', 'NM_015162.4'];

    Returns the dictionary with the lists without the elements that contain string_to_remove
    """
    result_dict = {}
    for key in src_dict:
        current_list = src_dict[key]
        current_result_list = []
        for element in current_list:
            if string_to_remove not in element:
                current_result_list.append(element)
        result_dict[key] = current_result_list
        current_result_list = []
    return result_dict


def remove_list_elements(src_list, string_to_remove):
    """
    Removes all the elements which contain 'string_to_remove' from the src_list
    Returns the list without the removed elements.
    """
    result = []
    for element in src_list:
        if string_to_remove not in element:
            result.append(element)
    return result
    

def _return_ensembl_from_id_and_uniprot_query(uniprotId, query):
    logger.debug(f"Starting retrival of ensemblId for uniprotId {uniprotId}")
    index = next((index for (index, d) in enumerate(query["results"]) if d["primaryAccession"] == uniprotId), None)
    ensembl_index_list=[]
    xref_arr_length = len(query["results"][index]["uniProtKBCrossReferences"]) # the count of cross-referenced databases
    for i in range(xref_arr_length):
        if query["results"][index]["uniProtKBCrossReferences"][i]["database"] == "Ensembl":
            ensembl_index_list.append(i)

    if len(ensembl_index_list) == 0:
        enId = None
    elif len(ensembl_index_list) == 1:
        enId=query["results"][index]["uniProtKBCrossReferences"][ensembl_index_list[0]]["id"].split(".")[0]
    elif len(ensembl_index_list) > 1:
        if any("isoformId" in query["results"][index]["uniProtKBCrossReferences"][i] for i in ensembl_index_list):
            for i in ensembl_index_list:
                if "-1" in query["results"][index]["uniProtKBCrossReferences"][i]["isoformId"]:
                    enId=query["results"][index]["uniProtKBCrossReferences"][i]["id"].split(".")[0]
        else:
            enId=query["results"][index]["uniProtKBCrossReferences"][1]["id"].split(".")[0] 
                
    logger.info(f"uniprotId {uniprotId} -> ensemblId {enId}")
    return enId

def _uniprot_query_API(gene_name, type="gene"):
        """
        Finds Uniprot Ids belonging to gene_name, returns result as json
        """
        logger.debug(f"Starting uniprot info query API for gene name {gene_name}")
        if type == "gene":
            _uniprot_identifier_query_result = requests.get(f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene_name}+AND+organism_id:9606&format=json&fields=accession,gene_names,organism_name,reviewed,xref_ensembl")
        elif type == "prot":
            _uniprot_identifier_query_result = requests.get(f"https://rest.uniprot.org/uniprotkb/search?query={gene_name}+AND+organism_id:9606&format=json&fields=accession,gene_names,organism_name,reviewed,xref_ensembl,protein_name")
        if _uniprot_identifier_query_result.text == "{'results': []}":
            # empty response, may be because ortholog was found, but human ortholog has different name than the gene from the file - denotes an error in file parsing
            raise Exception(f"No uniprot identifier query result for {gene_name} found.")
        logger.debug(f"type = {type}, query result json: {_uniprot_identifier_query_result.json()}")
        return _uniprot_identifier_query_result.json()
    
def get_uniprotId_description(uniprotId):
    """
    Returns a description of the specified uniprotId. 
    Example: https://rest.uniprot.org/uniprotkb/O14944.txt, this function parses lines after '-!- FUNCTION'
    
    Parameters:
      - uniprotId: either just the Id (e.g. O14944) or a prefixed Id (e.g. UniProtKB:Q9BUL8)
    """
    # TODO: You can also get interaction data, for example: uniprotId O14944 also Interacts with EGFR and ERBB4
    # -> parse this from "-!- SUBUNIT" part of the response text
    response=""
    if ":" in uniprotId:
        id = uniprotId.split(":")[1]
        response = requests.get(f"https://rest.uniprot.org/uniprotkb/{id}.txt")
    else:
        response = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprotId}.txt")
    split = response.text.split("\n")
    result_lines = []
    _readline_flag = False
    for line in split:
        if _readline_flag == False and "-!-" in line and "FUNCTION" in line:
            _readline_flag = True
        if _readline_flag == True and "-!-" in line and "FUNCTION" not in line:
            _readline_flag = False
        if _readline_flag == True:
            result_lines.append(line.replace("CC", "").strip())
    uniprotId_description = " ".join(result_lines).replace("-!-","").strip()
    return uniprotId_description

def get_uniprotId_from_geneName_new(gene_name, trust_genes=True):
    """
    Retrieves UniProt Identifier from a gene symbol/name; e.g. UniProtKB:Q86SQ4, if gene_name=adgrg6. 
    This function no longer uses recursion and also checks for verification status of the returned Ids.

    Parameters:
      - gene_name: A gene name or symbol e.g. ADGRG6
      - prefix: A database prefix, is prefixed before the specifi uniprot gene id
      - trust_genes: If True, all trusted genes inside genes_trusted.txt will override this function. If false, it
        notifies you that the gene was found among trusted genes and if you wish to proceed with the function.
    """
    
    def _get_uniprot_identifier_json_nth_response(json, nth_response):
        "Gets the entire nth element in 'results' array inside the json retrieved by get_uniprot_identifier function"
        return json["results"][nth_response]

    prefix="UniProtKB:" # do we need it?

    if gene_name in constants.TRUSTED_GENES:
        if trust_genes == True:
            # return the element ahead of the gene_name, which is the previously found uniprot_id
            logger.info(f"Gene {gene_name} is found among trusted genes.")
            uniprotId = constants.TRUSTED_GENES[constants.TRUSTED_GENES.index(gene_name)+1]
            return prefix+uniprotId
        else:
            trust_genes_user_response = int(input(f"Gene {gene_name} is found among trusted genes. Press 1 to continue with UniProt Id query, or 0 to use the trusted gene."))
            if trust_genes_user_response == 0:
                uniprotId = constants.TRUSTED_GENES[constants.TRUSTED_GENES.index(gene_name)+1]
                return prefix+uniprotId

    _uniprot_query_result = _uniprot_query_API(gene_name)

    uniprot_geneIds_dictionary = {} # stores uniprotId : reviewed_status pairs
    results_arr_len = len(_uniprot_query_result["results"])

    # if only one result, auto accept it
    if results_arr_len == 1:
        uniprotId = _uniprot_query_result["results"][0]["primaryAccession"]
        logger.info(f"Auto accepted {gene_name} -> {uniprotId}. Reason: Only 1 result.")
        if prefix != "": return prefix + uniprotId
        else: return uniprotId

    for i in range(results_arr_len):
        # reviewed is either TrEMBL or SwissProt. TrEMBL should be set to False and SwissProt reviews should be set to True
        is_reviewed = False if "TrEMBL" in _uniprot_query_result["results"][i]["entryType"] else True
        # has transcript is either True or False. We want to have transcripts.
        has_transcript = True if _uniprot_query_result["results"][i]["uniProtKBCrossReferences"] else False

        geneId = _uniprot_query_result["results"][i]["primaryAccession"]
        uniprot_geneIds_dictionary[geneId] = {"is_reviewed":is_reviewed, "has_transcript":has_transcript}# add to dict
    logger.info(f"uniprot_geneIds_dictionary = {uniprot_geneIds_dictionary}")

    # Used to check the gene_names field of the request against the gene_name parameter supplied by the function. If gene_name (parameter) is not
    # among gene_names, then such element should not be analysed further.
    # "As per 'uniprot-mapping-multiple-results-issue': implemented functionality to dismiss results that have different genes as the specified gene (idk why uniprot returns such genes nevertheless)"
    geneIds_impactGenes_dictionary = {} # stores <str> uniprotId : <list> impact_genes; used later to check that uniprotId contains relevant gene
    for j in range(results_arr_len):
        impactGenes_all = []
        impactGenes_synonyms = "" # init to prevent errors
        impactGenes_primary = _uniprot_query_result["results"][j]["genes"][0]["geneName"]["value"]
        if "synonyms" in _uniprot_query_result["results"][j]["genes"][0]:
            impactGenes_synonyms = _uniprot_query_result["results"][j]["genes"][0]["synonyms"][0]["value"]
        if isinstance(impactGenes_primary, list):
            for g in impactGenes_primary: impactGenes_all.append(g)
        else: impactGenes_all.append(impactGenes_primary)
        if isinstance(impactGenes_synonyms, list):
            for g in impactGenes_synonyms: impactGenes_all.append(g)
        else: impactGenes_all.append(impactGenes_synonyms)
        geneIds_impactGenes_dictionary[get_dict_key_at_index(uniprot_geneIds_dictionary,j)] = impactGenes_all
    logger.debug(f"geneIds_impactGenes_dictionary = {geneIds_impactGenes_dictionary}")

    # Eliminate entries where gene_name is not among "impact genes" - final check to mitigate 'uniprot-mapping-multiple-results-issue'
    # Also checks for autoaccept (if onyl one of the UniprotIds is reviewed and has transcript) to prevent multiple for loop reiterations
    pop_keys = [] # indices of items to eliminate
    pop_keys_indices = []
    reviewedId_single = ""
    NO_reviewed_Ids = 0
    for key, value in geneIds_impactGenes_dictionary.items():
        _d_review_status = uniprot_geneIds_dictionary[key]["is_reviewed"]
        _d_has_transcript = uniprot_geneIds_dictionary[key]["has_transcript"]
        logger.debug(f"Autoaccept process: key = {key}, value = {value}, review status = {_d_review_status}, has_transcript = {_d_has_transcript}")
        if gene_name not in value: 
            pop_keys.append(key)
            pop_keys_indices.append(i)
        else: # gene_name is found among geneNames request return field -
            # if uniprot_geneIds_dictionary[key]["is_reviewed"] == True and uniprot_geneIds_dictionary[key]["has_transcript"] == True: # check if UniprotId is reviewed #TODO: check AND as condition here!
            if uniprot_geneIds_dictionary[key]["is_reviewed"] == True:    
                if uniprot_geneIds_dictionary[key]["has_transcript"] == True:
                    logger.debug(f"{key} has_transcript = {_d_has_transcript}")
                    # TODO: check the meaning of has_transcript field on the code execution
                    # Some elements have reviewed = True and has_transcript = False -> meaning ?
                if reviewedId_single == "": # check if no other reviewed id was added before
                    reviewedId_single = key
                    NO_reviewed_Ids += 1

    # eliminate faulty elements in uniprot_geneIds_dictionary 
    for pop_key in pop_keys:
        uniprot_geneIds_dictionary.pop(pop_key)
        geneIds_impactGenes_dictionary.pop(pop_key)
    logger.debug(f"removed {len(pop_keys)} faulty ids (gene_name did not match with geneNames field). New uniprot_geneIds_dictionary len = {len(uniprot_geneIds_dictionary)}")

    # Autoaccept if only one of the UniprotIds is reviewed
    if NO_reviewed_Ids == 1:
        logger.info(f"Found a single reviewed UniProt Id with a transcript for gene_name {gene_name}: {reviewedId_single}")
        if prefix != "": return prefix+reviewedId_single
        else: return reviewedId_single
    
    # Return None if there are no reviewed ids (TODO: implement a flag for this)
    # TODO: implement a file where you store such instances
    if NO_reviewed_Ids == 0:
        logger.info(f"No reviewed UniProt Ids for gene_name {gene_name} found, returning None")
        return None
    
    # Either 2+ or 0 verified UniprotIds -> begin user cycling options
    results_arr_len = len(uniprot_geneIds_dictionary) # update after eliminations
    i = 0
    next_step = 1
    for uprId, info in uniprot_geneIds_dictionary.items(): # used camelcase for differentiating inter-for-loop variables (from function variables)
        is_reviewed=info["is_reviewed"]
        has_transcript=info["has_transcript"]
        logger.info(f"Gene name {gene_name} found to correspond to {uprId} (reviewed = {is_reviewed}, has transcript = {has_transcript}). Displaying response {i}/{results_arr_len}: {_get_uniprot_identifier_json_nth_response(_uniprot_query_result,i)}")
        logger.info(get_uniprotId_description(uprId))
        
        user_logic = int(input("Press 1 to confirm current result, 2 to cycle another result, or 0 to continue the program and discard all options."))
        if user_logic == 1: # result is confirmed, save so in TRUSTED_GENES
            logger.info(f"Confirmed {gene_name} -> {uprId}")
            append_to_file(f"{gene_name} {uprId}\n", "src_data_files/genes_trusted.txt")
            # return
            # TODO: handle a case if a user wants to confirm multiple reviewed proteins -> store them in a list and return when next_step > (results_arr_len-1)
            if prefix != "": return prefix + uprId
            else: return uprId
        elif user_logic == 2: #cycle next result
            if next_step > (results_arr_len - 1):
                logger.info("Cycled out of options!")
                return f"CycleOutOfBoundsError: Cycled through all the uniprot gene IDs for {gene_name} without confirming any"
            i += 1 
            next_step += 1
            continue # skip to the next for loop element
        elif user_logic == 0: #TODO: debate the use of this, add return / redo / start loop over functionality
            logger.info("No uniprot geneIds selected.")
            return f"No uniprot gene Ids for {gene_name} selected."
        else:
            # wrong input
            raise Exception(f" Wrong input!")
    logger.info("No uniprot geneIds selected")
    return f"No uniprot gene Ids for {gene_name} selected."

def zfin_find_human_ortholog(gene_id):
    """
    If gene_id is from the ZFIN database, searches through the zebrafish-human orthologs and returns the name of the
    symbol of the human gene ortholog.

    Returns:
      - [0]: gene symbol
      - [1]: long name of the gene
    """
    def _zfin_get_human_gene_symbol_from_line(line, improved_algorithm=True):
        """
        Splits zfin line and retrieves human gene symbol (full caps of zebrafish gene symbol)
        """
        if improved_algorithm == True:
            # better, as zfin human ortholog sometimes has different name than the zebrafish gene
            human_symbol = line.split("\t")[3]
            human_gene_name = line.split("\t")[4]
            return human_symbol, human_gene_name # look at zfin orthologs txt file (in src_data_files) -> when you higlight a row, you see a TAB as '->' and a SPACEBAR as '.' -> splitting at \t means every [3] linesplit element is the human gene name
        else:
            human_symbol = str(line.split("\t")[1]).upper()
            human_gene_name = str(line.split("\t")[4])
            return human_symbol, human_gene_name # split lines at tabs (spacebar is not ok!)

    gene_id=gene_id.split(":")[1] # eliminates 'ZFIN:' 
    for line in _zfin_ortholog_readlines:
        if gene_id in line:
            e = _zfin_get_human_gene_symbol_from_line(line, improved_algorithm=True)
            human_symbol = e[0]
            human_gene_name = e[1]
            logger.info(f"[ Returning human symbol {human_symbol} and {human_gene_name}")
            return human_symbol, human_gene_name
    return [f"ZfinError_No-human-ortholog-found:gene_id={gene_id}"]

def xenbase_find_human_ortholog(gene_id):
    """
    Attempts to find a human ortholog from the xenbase database.
    Parameters:
      - gene_id: eg. Xenbase:XB-GENE-495335 or XB-GENE-495335
    Returns: 
      - [0]: symbol of the human ortholog gene (eg. rsu1) or 'XenbaseError_no-human-ortholog-found'
      - [1]: long name of the gene
    """
    def _xenbase_get_human_symbol_from_line(line):
        """Splits xenbase line at tabs and gets human gene symbol (in full caps)"""
        symbol = str(line.split("\t")[2]).upper()
        name = str(line.split("\t")[3])
        return symbol, name

    gene_id_short = ""
    if ":" in gene_id: gene_id_short = gene_id.split(":")[1]
    else: gene_id_short = gene_id
    
    for line in _xenbase_ortholog_readlines:
        if gene_id_short in line:
            e = _xenbase_get_human_symbol_from_line(line)
            human_symbol = e[0]
            human_gene_name = e[1]
            logger.info(f"Found human ortholog {human_symbol}, name = {human_gene_name} for xenbase gene {gene_id}")
            return human_symbol, human_gene_name
    return [f"[XenbaseError_No-human-ortholog-found:gene_id={gene_id}"]

def mgi_find_human_ortholog(gene_id):
    """
    Attempts to find a human ortholog from the mgi database.
    Parameters: gene-id eg. MGI:MGI:98480
    Returns: symbol of the human ortholog gene or "MgiError_no-human-ortholog-found".
    
    Note: Cannot return longer gene name from the MGI .txt file, since it doesn't contain the longer name
    """
    def _mgi_get_human_symbol_from_line(line, line_index):
        """
        Splits mgi line at tabs and gets human gene symbol
        """
        split = line.split("\t")
        if split[1] != "human":
            # try i+2 to check one line further down
            line = _mgi_ortholog_readlines[line_index+2]
            second_split = line.split("\t")
            if second_split[1] == "human":
                logger.debug(f"Found keyword 'human' on secondpass line querying.")
                return second_split[3]
            else:
                # this still means no human ortholog!
                # example: MGI:2660935 (Prl3d2) contains no "human" (neither i+1 nor i+2), also checked uniprot and no human gene for prl3d2 exists
                return f"[MgiError_No-human-ortholog-found:gene_id={gene_id}"
        else: return split[3]

    logger.debug(f"Starting MGI search for {gene_id}")
    gene_id_short = ""
    if ":" in gene_id:
        split = gene_id.split(":")
        if len(split) == 3: gene_id_short = split[2] # in case of MGI:xxx:xxxxx
        elif len(split) == 2: gene_id_short = split[1] # in case of MGI:xxxxx
    else: gene_id_short = gene_id

    i = 0
    for line in _mgi_ortholog_readlines:
        if gene_id_short in line:
            # if "mouse" gene smybol is found at line i, then human gene symbol will be found at line i+1
            logger.debug(f"i = {i}, gene_id_short = {gene_id_short}, line = {line}")
            human_symbol = _mgi_get_human_symbol_from_line(_mgi_ortholog_readlines[i+1], i)
            logger.info(f"Found human ortholog {human_symbol} for mgi gene {gene_id}")
            return human_symbol # return here doesnt affect line counter 'i', since if gene is found i is no longer needed
        i += 1
    return f"[MgiError_No-human-ortholog-found:gene_id={gene_id}"

def rgd_find_human_ortholog(gene_id):
    """ 
    Attempts to find a human ortholog from the RGD (rat genome database) 
    Returns: human gene symbol

    Note: longer name of the gene cannot be returned, since it is not specified in the rgd txt file
    """
    def _rgd_get_human_symbol_from_line(line):
        """ Splits rgd line at tabs and gets human gene smybol """
        # also clears whitespace from linesplit (which is split at tab). Some lines in RGD db text file had whitespace instead of \t -> clear whitespace from array to resolve
        # example: linesplit = ['Ang2', '1359373', '497229', '', '', '', '', 'Ang2', '1624110', '11731', 'MGI:104984', 'RGD', '\n']
        linesplit = line.split("\t")
        result_list = [] 
        for element in linesplit: 
            if element != "":
                result_list.append(element)
        return result_list[3]

    gene_id_short = ""
    if ":" in gene_id: gene_id_short = gene_id.split(":")[1]
    else: gene_id_short = gene_id

    i = 0
    for line in _rgd_ortholog_readlines:
        if gene_id_short in line:
            splitline_debug = line.split("\t")
            human_symbol = _rgd_get_human_symbol_from_line(line)
            logger.info(f"Found human ortholog {human_symbol} for RGD gene {gene_id}")
            return human_symbol
    return f"[RgdError_No-human-ortholog-found:gene_id={gene_id}"

def load_human_orthologs():
    """
    This function should be called at runtime once to load the ortholog mapping txt files into proper variables.
    """
    global _uniprot_genes_readlines
    _uniprot_genes_readlines = readlines("src_data_files/uniprotkb_human_idmapping.dat") # from https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/, file HUMAN_9606_idmapping.dat.gz  
    global _zfin_ortholog_readlines
    _zfin_ortholog_readlines = readlines("src_data_files/zfin_human_ortholog_mapping.txt")
    global _xenbase_ortholog_readlines
    _xenbase_ortholog_readlines = readlines("src_data_files/xenbase_human_ortholog_mapping.txt")
    global _mgi_ortholog_readlines
    _mgi_ortholog_readlines = readlines("src_data_files/mgi_human_ortholog_mapping.txt")
    global _rgd_ortholog_readlines
    _rgd_ortholog_readlines = readlines("src_data_files/rgd_human_ortholog_mapping.txt")

def find_term_corresponding_array(term):
    """
    Finds the corresponding array from constants.py for the given GO term. 
    Parameters:
      - term: A GO term; GO:xxxxx or GO-xxxxx
    
    Returns: a list with 2 elements:
      - [0]: array enum from constants.TERM_ARRAY_ENUMS_SINGLE, specifying the exact array enum (either angio+, angio- or angio_general, dia+, dia- or dia_general)
      - [1]: array enum from constants.TERM_ARRAY_ENUMS_MULTIPLE, specifying the general process (either angiogenesis or diabetes)
    """
    result = []
    logger.debug(f"term type = {type(term)}, term = {term}")
    if "-" in term: 
        term = term.replace("-",":") # replace - for : which is used in constants.py
    logger.debug(f"len terms single = {len(constants.TERM_ARRAYS_ENUMS_SINGLE)}")
    for array_enum, terms in constants.TERM_ARRAYS_ENUMS_SINGLE.items():
        if term in terms:
            result.append(array_enum)
    if len(result) == 0: result.append("notfound") # if no items were added in previous for loop, append "/"
    
    # len(result) is now definitely 1
    for array_enum, terms in constants.TERM_ARRAY_ENUMS_MULTIPLE.items():
        if term in terms:
            result.append(array_enum)
    if len(result) == 1: result.append("notfound") # if no items were added in previous for loop, append "/"
    logger.debug(f"term {term} array enum query: [0] = {result[0]}, [1] = {result[1]}")
    return result
        
def handle_load_from_crash(full_list, crash_filepath, crash_flag, additional_lists=[]):
    crash_json = read_file_as_json(crash_filepath)
    if json.dumps(crash_json) != "[]": # if crash_json isn't empty
        if crash_flag == 1:
            last_geneId = get_last_geneId_in_crash_json(crash_json)
            logger.info(f"Crash recovery: last_geneId = {last_geneId}, genes_len = {len(full_list)}")
            return list_directionshrink(full_list, last_geneId, forward=True)
        if crash_flag == 2:
            # TODO: PROCESS LADI'S NEW ALGORITHM TO SHORTEN ARRAY OF PRODUCTS
            # RETURN BOTH DIRECTIONSHRINK FOR TERMS AND FOR LAST PRODUCT !!! (add products list in additional_lists)
            last_element = crash_json[-1]["products"][-1]
            direct_productId = last_element["direct_productID"]
            uniprotId = last_element["UniprotID"]

def get_last_product_in_crash_json(crash_json):
    """
    Used in crash recovery. Returns:
      - [0]: GO_term: the last GO term before the crash happened
      - [1]: direct_productId: the direct product Id (can be ZFIN, RGD etc)
      - [2]: uniprotId: the uniprot id of the last product in crash_json
    """
    GO_term = crash_json[-1]["GO_term"]
    last_element = crash_json[-1]["products"][-1]
    direct_productId = last_element["direct_productID"]
    uniprotId = last_element["UniprotID"]
    return GO_term, direct_productId, uniprotId

def get_pre_last_product_in_crash_json(crash_json):
    """
    Used in crash recovery. Returns the element before the last element in crash json similarly to get_last_product_in_crash_json
    """
    # TODO: CHECK IF CRASH_JSON CAN RETURN -2 EG IF LEN(CRASH_JSON) > 1
    if len(crash_json) > 1:
        GO_term = crash_json[-2]["GO_term"]
        last_element = crash_json[-2]["products"][-1]
        direct_productId = last_element["direct_productID"]
        uniprotId = last_element["UniprotID"]
        return GO_term, direct_productId, uniprotId
    else: return [-1]

def choose_crashfile(crash_filepaths):
    """
    Gives user the choise to choose a crash file among crash_filepaths
    """
    display_dictionary = {}
    i = 0
    for path in crash_filepaths:
        display_dictionary[i] = path
        i += 1
    choice = int(input(f"Enter the number of the crashfile from {display_dictionary}"))
    return display_dictionary[choice]

def uniprot_find_human_symbol(gene_id):
    """
    Finds a gene symbol (eg. adgrg) from gene_id (UniProtKB:XXXXX)
    """
    if ":" in gene_id: gene_id = gene_id.split(":")[1]
    for line in _uniprot_genes_readlines:
        if gene_id in line:
            line = line.replace("\n", "")
            split = line.split("\t")
            gene_symbol = split[2]
            return gene_symbol
    return "UniprotError - Gene symbol not found."

def product_mRNA_json_append_refseqIds(product_mRNA_json_filepath=constants.TARGET_FOLDER+"/product_mRNA.json", result_filepath=constants.TARGET_FOLDER+"/product_mRNA_refseq.json"):
    """
    Appends the NCBI accession ids (eg. NM_001008395.3) to each element inside product_mRNA. NCBI accession ids are necessary for intercompatibility with
    the miRDB_v6.0_prediction_result.txt file (which contains miRNA match strengths against known mRNAs (mRNAs presented in NCBI accession format)
    """
    product_mRNA_json = read_file_as_json(product_mRNA_json_filepath)
    uniprot_human_idmapping = readlines("src_data_files/uniprotkb_human_idmapping.dat")

    # uniprotkb_human_idmapping: if stripped uniprot id is in line and if RefSeq_NT is the second line element, append third line element
    # (refseq id) to array
    uniprotIds_to_refseqIds = []
    i = 0
    for element in product_mRNA_json:
        logger.info(f"Processing {i}")
        uniprotId_stripped = element["UniprotID"].split(":")[1] # UniProtKB:Q0VGL1 -> Q0VGL1
        for line in uniprot_human_idmapping:
            if uniprotId_stripped in line and "RefSeq_NT" in line:
                uniprotIds_to_refseqIds.append(line.split("\t")[2].replace("\n", "")) # append the third line element (refseq id)
        
        # after all lines are read: append uniprotIds_to_refseqIds_dict to the current UniProtKB Id, reset the list
        product_mRNA_json[i]["RefSeq_NT_IDs"] = uniprotIds_to_refseqIds # create RefSeq_NT_IDs field for i'th element in product_mRNA and store refseqids
        uniprotIds_to_refseqIds = []
        i += 1
    
    # save new json
    save_json(product_mRNA_json, result_filepath)

def find_gene_symbol_and_name_from_id(gene_id):
    if "UniProtKB" in gene_id:
        gene_symbol = uniprot_find_human_symbol(gene_id)
        gene_longname = "" # todo: need to access uniprot_sprot.xml, which takes ages to load
        return gene_symbol, gene_longname
    elif "ZFIN" in gene_id:
        gene_symbol = zfin_find_human_ortholog(gene_id)[0]
        gene_longname = zfin_find_human_ortholog(gene_id)[1]
        return gene_symbol, gene_longname
    elif "Xenbase" in gene_id:
        gene_symbol = xenbase_find_human_ortholog(gene_id)[0]
        gene_longname = xenbase_find_human_ortholog(gene_id)[0]
        return gene_symbol, gene_longname
    elif "RGD" in gene_id:
        gene_symbol = rgd_find_human_ortholog(gene_id) # not [0] because MGI and RGD cannot return gene longname
        return gene_symbol, ""
    elif "MGI" in gene_id:
        gene_symbol = mgi_find_human_ortholog(gene_id)
        return gene_symbol, ""

def scoring_results_postprocess(score_results_filepath):
    """
    Adds product names and descriptions.
    """
    score_results_json = read_file_as_json(score_results_filepath)
    final_json_elements = []
    for i,element in enumerate(score_results_json):
        id = element["gene"]
        count = element["count"]
        symbol = ""
        long_name = ""
        description = ""
        terms = element["terms"] 
        term_enum_scores = element["term_enum_scores"]

        symbol = find_gene_symbol_and_name_from_id(id)[0]
        long_name = find_gene_symbol_and_name_from_id(id)[1]
        if "UniProtKB" in id:
            description = get_uniprotId_description(id)
        else:
            description = get_uniprotId_description(get_uniprotId_from_geneName_new(symbol))
        out = {"gene": id, "count": count, "symbol": symbol, "long_name": long_name, "description": description, "terms": terms, "term_enum_scores": term_enum_scores}
        final_json_elements.append(out)
    
    _fn = score_results_filepath.replace(".json", "")
    save_json(final_json_elements, f"{_fn}_postprocess.json")

def get_ensembl_ids_from_uniprot_ids(gene_list):
    """
    Gets a gene_list of UniprotIds, returns a list of appropriate EnsemblIds
    """
    def _get_ensembl_id_from_uniprot_id(gene):
        gene_id = gene.split(":")[1]
        up_query = _uniprot_query_API(gene_id, type="prot")
        ensembl_id = _return_ensembl_from_id_and_uniprot_query(gene_id, up_query)
        return ensembl_id

    ensembl_ids=[]
    for gene in gene_list:
        ensembl_ids.append(_get_ensembl_id_from_uniprot_id(gene))
    return ensembl_ids

def extract_n_elements_from_json(json_filepath, n, destination_filepath):
    """
    Takes n elements from the start of the json and saves them into destination_filepath. Also returns computed json object.
    Note: n should be from 1-infinity. The 0 convention is handled in code. So if n=10, 10 elements are returned (from 0 to 9)
    """
    json_object = read_file_as_json(json_filepath)
    results = []
    for index, element in enumerate(json_object):
        if index < n:
            results.append(element)
        else:
            break
    store_json_dictionaries(destination_filepath, results)
    return results

def get_uniprotids_from_json(json_filepath, uniprot_id_identifier="gene"):
    """
    Note: this function is deprecated. Rather use get_identifier_values_from_json with the "gene" identifier
    Parameters:
      - json_filepath
      - uniprot_id_identifier: the json identifier which holds the value of the element's uniprotId

    Returns:
      - [0]: a list of uniprotIds of elements inside json_filepath
      - [1]: the entire json object
    """
    #json_obj = read_file_as_json(json_filepath)
    #uniprotIds = []
    #for el in json_obj:
    #    uniprotIds.append(el[uniprot_id_identifier])
    #logger.debug(f"UniprotIds: {uniprotIds}")
    #return uniprotIds, json_obj
    return get_identifier_values_from_json(json_filepath, uniprot_id_identifier)

def get_identifier_values_from_json(json_filepath, identifier):
    """
    Loops through all top-level json elements, queries identifier of each element and appends the value to list

    Returns:
      - [0]: list of identifier values
      - [1]: entire json obj

    Example: you want to get all UniprotIds (of all elements) inside a json -> use get_identifier_values_from_json(filepath, "gene")
    """
    json_obj = read_file_as_json(json_filepath)
    identifier_values = []
    for element in json_obj:
        identifier_values.append(element[identifier])
    # logger.debug(f"Values for identifier {identifier}: {identifier_values}") #can clog the console
    logger.debug(f"Returned {len(json_obj)} elements")
    return identifier_values, json_obj

def read_embl_file_records(embl_filepath):
    """
    Reads the records inside an embl filepath (either a .embl file or a .dat file with records stored in the embl format) into a list
    of records, each record can have it's fields queried by record.id, record.seq (for sequence), record.annotation etc.

    Returns: a list of embl records
    """
    records = list(SeqIO.parse(embl_filepath, "embl"))
    return records

def save_mirbase_hsap_miRNAs_for_cpp(embl_mirbase_download_filepath, destination_filepath = "src_data_files/miRNAdbs/mirbase_miRNA_hsa-only.txt", convert_uracil_thymine = True):
    """
    Loops through the records in embl_mirbase_download_filepath, saves the miRNAs with the 'Homo Sapiens' in their DE (record.desc) field
    (note: a Homo-Sapiens miRNA will also contain the hsa prefix in record.name) in the destination_filepath.txt file in the format:

        AC (record.id), ID (record.name), SEQ (record.seq), base pair count

        MI0000060 \t hsa-let-7a-1 \t ugggaugagg uaguagguug uauaguuuua gggucacacc caccacuggg agauaacuau acaaucuacu gucuuuccua \t 80
    
    Parameters:
      - @param embl_mirbase_download_filepath: the filepath to a .dat file downloaded from https://mirbase.org/ftp.shtml (eg. src_data_files/miRNAdbs/mirbase_miRNA_06-02-2023.dat)
      - @param destination_filepath: the filepath to the .txt file that will contain the parsed homo sapiens miRNAs only
      - @param convert_uracil_thymine: if true, it converts uracil (U) in a sequence (if Us are detected) to thymine (T). If false, uracils are left unchanged.
    """
    if os.path.exists(destination_filepath):
        user_in = input(f"{destination_filepath} already exists. Press 0 to abort or 1 to process the mirbase embl raw file again.")
        if user_in == 0:
            logger.info("Aborted.")
            return
        logger.info(f"Processing {destination_filepath}")

    with open(destination_filepath, "a") as f:
        records = read_embl_file_records(embl_mirbase_download_filepath)
        count = len(records)
        i=0
        for record in records:
            logger.info(f"{i}/{count}")
            if "hsa" in record.name or "Homo" in record.description:
                sequence = record.seq
                if convert_uracil_thymine == True:
                    if "U" in sequence:
                        sequence = convert_thymine_uracil(sequence, 0)
                f.write(f"{record.id}\t{record.name}\t{sequence}\t{len(record.seq)}\n")
            i+=1

def save_mRNAs_for_cpp(product_mRNA_refseq_json_filepath, destination_filepath=constants.TARGET_FOLDER+"/product_mRNAs_cpp.txt"):
    """
    Loops through the product_mRNA_refseq_json_filepath file (usually stored in test_run_X/product_mRNA_refseq_json). From each json element,
    UniProtID, mRNA sequence and refseq_NT_ids are extracted and stored on the same line. 

    File structure:
    UniProtId \t mRNA_sequence \t refseq_NT_id 1 \t refseq_NT_ID 2 \t ...

    As you can see, the first and second line elements will always be UniProtId and mRNA sequence, respectively. Then, one ore more
    refseq_NT_ids may follow (some mRNAs have 8+ refseq IDs which denote different transcription variants)

    The mRNA_refseq file is created in the current constants.TAGRET_FOLDER after running dev_mrna_download.py
    """
    # handle product_mRNA_refseq_json_filepath
    if not check_file_exists(product_mRNA_refseq_json_filepath):
        user_input = input(f"File {product_mRNA_refseq_json_filepath} does not exist. Press 0 to abort (and make sure to run dev_mrna_download.py) next time. Press 1 to write a custom relative path to a mrna_refseq json file.")
        if user_input == 0:
            return -1
        elif user_input == 1:
            user_input = input(f"Write the relative path to mrna refseq file and press enter:")
            product_mRNA_refseq_json_filepath = user_input
            if not check_file_exists(product_mRNA_refseq_json_filepath): return -1
        else:
            logger.info("Invalid input. Can only be 0 or 1.")
            return -1

    # handle destination_filepath
    if os.path.exists(destination_filepath):
        user_input = input(f"File {destination_filepath} already exists. Press 0 to continue and overwrite the file, 1 to abort or 2 to change the name of destination_filepath")
        # if user_input == 0, do nothing
        if user_input == 1:
            logger.info("Aborting")
            return -1
        elif user_input == 2:
            dest_filepath_input = input(f"Write the relative filepath to the destination file and press enter:")
            destination_filepath = dest_filepath_input
            if os.path.exists(destination_filepath): 
                logger.info("That destination path exists too. Aborting")
                return -1

    # read each element into a single line representation
    logger.info(f"Determined mRNA refseq json filepath: {product_mRNA_refseq_json_filepath}")
    jsonObj = read_file_as_json(product_mRNA_refseq_json_filepath)
    result_lines = []
    for element in jsonObj: # read each element into a single line, separator = tab (\t)
        # debug
        element_uniprot_id = element["UniprotID"]
        element_mRNA_seq = element["mRNA"]
        element_refseq_NT_ids = element["RefSeq_NT_IDs"]
        logger.debug(f"unipr = {element_uniprot_id}, mRNA = {element_mRNA_seq}, refseq = {element_refseq_NT_ids}")

        # handle None (bug fix)
        if element_uniprot_id == None: element_uniprot_id = "None"
        if element_mRNA_seq == None: element_mRNA_seq = "None"

        # main line creation logic
        line_string = element_uniprot_id + "\t" + element_mRNA_seq # uniprotid and mRNA_sequence are the first two line elements
        for refseq_NT_id in element["RefSeq_NT_IDs"]: # append all refseqs from the NCBI Nucleotide database (that are not XM, since XM are not confirmed)
            if "XM" not in refseq_NT_id:
                line_string = line_string + "\t" + refseq_NT_id
        line_string = line_string + "\n"
        result_lines.append(line_string)
    
    # write the lines to destination_filepath
    f = open(destination_filepath, "w")
    f.writelines(result_lines)
    f.close()


def get_time():
    return time.time()

def compute_time(first_time, print_message=False):
    """
    Computes and displays time between first_time and now

    Parameters:
      - first_time: use get_time() at a start of a function, then call compute_time with that time from the end of the function
      - print_message: if True, prints message to logger.debug
    """
    now = time.time()
    if print_message: logger.debug("%s seconds" % (now - first_time))
    return "%s seconds" % (now - first_time)

def compare_miRNA_mRNA_match_strength_single(miRNA,mRNA, debugLog = True):
    """
    Compares the match strength of miRNA and mRNA. This function contains the same algorithm which the CUDA kernels use
    to brute-force compute the best match miRNA-mRNA match strength. Algorithm anneals entire miRNA nucleotide sequence on
    all possible mRNA starting nucleotides.

    @param miRNA: a miRNA sequence or id (a MI or HSA id)
    @param mRNA: a mRNA sequence or id (uniprot id)
    @param debugLog: if True, will print the forward and backward comparison scores.

    @return: the best match strength between miRNA and mRNA
    """
    
    def _run_miRNA_mRNA_match_strength_compare_logic(miRNA,mRNA):
        miRNA_size = len(miRNA)
        mRNA_size = len(mRNA)

        _max_match_strength = 0.0
        i = 0
        total_opcounts = 0
        for character in mRNA:
            if i > (mRNA_size - miRNA_size):
                break
            mRNA_substring = mRNA[i:int(i+miRNA_size)] # CUDA block does this
            num_matches = 0
            for j in range(0, miRNA_size): # CUDA thread does this
                total_opcounts += 1
                if (miRNA[j] == 'A' and mRNA_substring[j] == 'T') or (miRNA[j] == 'T' and mRNA_substring[j] == 'A') or (miRNA[j] == 'C' and mRNA_substring[j] == 'G') or (miRNA[j] == 'G' and mRNA_substring[j] == 'C'):
                    num_matches += 1
            current_match_score = num_matches/miRNA_size
            if current_match_score > _max_match_strength:
                _max_match_strength = current_match_score
            i += 1
        
        # processing is valid, if total_opcounts / ((mRNA_size-miRNA_size)*miRNA_size) == 1
        if total_opcounts / ((mRNA_size - miRNA_size)*miRNA_size) == 1:
            logger.debug("opcount score OK")
        else:
            logger.debug(f"total_opcounts = {total_opcounts}, mRNA_size = {mRNA_size}, miRNA_size = {miRNA_size} (diff = {mRNA_size - miRNA_size}), opscore = {total_opcounts / ((mRNA_size - miRNA_size)*miRNA_size)}")

        return _max_match_strength

    # check if miRNA is id
    _d_miRNA_id = ""
    _d_mRNA_id = ""
    if "hsa" in miRNA or "MI" in miRNA:
        _d_miRNA_id = miRNA
        miRNA = find_miRNA_sequence(miRNA)
    # check if mRNA is id
    if "Prot" in mRNA:
        _d_mRNA_id = mRNA
        mRNA = find_mRNA_sequence(mRNA)

    forward_match_strength = _run_miRNA_mRNA_match_strength_compare_logic(miRNA,mRNA) # forward comparison
    backward_match_strength = _run_miRNA_mRNA_match_strength_compare_logic(miRNA,mRNA[::-1]) # backward comparison - flip one sequence
    if debugLog == True:
        logger.debug(f"{_d_miRNA_id} :: {_d_mRNA_id} match scores:")
        logger.debug(f"    - forward: {forward_match_strength}")
        logger.debug(f"    - backward: {backward_match_strength}")

    result_match_strength = 0
    if forward_match_strength > backward_match_strength:
        result_match_strength = forward_match_strength
    else:
        result_match_strength = backward_match_strength
    return result_match_strength

    # miRNA_size = len(miRNA)
    # mRNA_size = len(mRNA)   
    # i = 0
    # maximum_match_strength = 0.0
    # for character in mRNA:
    #    if i > (mRNA_size - miRNA_size):
    #        break
    #    mRNA_substring = mRNA[i:int(i+miRNA_size)]
    #    num_matches = 0
    #    for i in range(0, miRNA_size):
    #        if mRNA_substring[i] == miRNA[i]:
    #            num_matches += 1
    #    current_match_score = num_matches/miRNA_size
    #    if current_match_score > maximum_match_strength:
    #        maximum_match_strength = current_match_score
    #    i += 1
    #return maximum_match_strength

def compare_miRNA_mRNA_match_strength_single_v2(miRNA,mRNA, debugLog = True, fileWrite = False, fileWrite_filepath_root = "debug_files"):
    """
    Compares the match strength of miRNA and mRNA. This function contains the same algorithm which the CUDA kernels use
    to brute-force compute the best match miRNA-mRNA match strength. Algorithm anneals entire miRNA nucleotide sequence on
    all possible mRNA starting nucleotides.

    @param miRNA: a miRNA sequence or id (a MI or HSA id)
    @param mRNA: a mRNA sequence or id (uniprot id)
    @param debugLog: if True, will print the forward and backward comparison scores.
    @param fileWrite: writes the operations to the file
    @param fileWrite_filepath_root: the root folder where the output will be saved

    @return: the best match strength between miRNA and mRNA
    """
    def _compare_sequences(miRNA_seq, mRNA_substring):
        j = 0
        successful_matches = 0
        all_characters = len(miRNA)
        for character in miRNA:
            if (miRNA_seq[j] == "A" and mRNA_substring[j] == "T") or (miRNA_seq[j] == "T" and mRNA_substring[j] == "A") or (miRNA_seq[j] == "C" and mRNA_substring[j] == "G") or (miRNA_seq[j] == "G" and mRNA_substring[j] == "C"):
                successful_matches += 1
            j+=1
            constants._d_total_opcounts += 1
        return successful_matches/float(all_characters)
    
    def _compare_sequences_debug(miRNA_seq, mRNA_substring, anneal_start_index, filepath):
        j = 0
        successful_matches = 0
        all_characters = len(miRNA)
        _d_matches_str = ""
        with open(filepath, "a") as f:
            for character in miRNA:
                if ((miRNA_seq[j] == 'A' and mRNA_substring[j] == 'T') or (miRNA_seq[j] == 'T' and mRNA_substring[j] == 'A') or (miRNA_seq[j] == 'C' and mRNA_substring[j] == 'G') or (miRNA_seq[j] == 'G' and mRNA_substring[j] == 'C')):
                    _d_matches_str += "T"
                    successful_matches += 1
                else:
                    _d_matches_str += " "
                j+=1
                constants._d_total_opcounts += 1
            # sequence_op is an operation performed on a single sequence
            f.write(f"seq_op {int(constants._d_total_opcounts/len(miRNA))}, op {constants._d_total_opcounts}: {successful_matches}/{all_characters} = {successful_matches/float(all_characters)}\n")
            f.write(f"miRNA: {miRNA_seq}\n")
            f.write(f"mRNA:  {mRNA_substring}\n")
            f.write(f"match: {_d_matches_str}\n")
            f.write("\n")
        
        return successful_matches/float(all_characters)
    
    # check if miRNA is id
    _d_miRNA_id = ""
    _d_mRNA_id = ""
    if "hsa" in miRNA or "MI" in miRNA:
        _d_miRNA_id = miRNA
        miRNA = find_miRNA_sequence(miRNA)

    # check if mRNA is id
    if "Prot" in mRNA:
        _d_mRNA_id = mRNA
        mRNA = find_mRNA_sequence(mRNA)

    _d_mRNA_id_rep = _d_mRNA_id.replace(":","") 
    fileWrite_filepath_p = os.path.join(fileWrite_filepath_root, f"{_d_miRNA_id}-{_d_mRNA_id_rep}_comparison_log.txt") # full filepath to outfile for log
    miRNA_size = len(miRNA)
    mRNA_size = len(mRNA)
    constants._d_total_opcounts = 0

    if(fileWrite == True):
        with open(fileWrite_filepath_p, "a") as f:
            f.write(f"# Comparison log for {_d_miRNA_id} - {_d_mRNA_id}\n")
            f.write(f"#\n")
            f.write(f"# miRNA_seq: {miRNA}\n")
            f.write(f"# mRNA_seq: {mRNA}\n")
            f.write("\n")

    i = 0
    max_match_strength = 0.0
    for character in mRNA:
        if i > (mRNA_size - miRNA_size):
            break
        
        mRNA_substring = mRNA[i:int(i+miRNA_size)]
        if fileWrite == False:
            match_strength_straight = _compare_sequences(miRNA, mRNA_substring)
            match_strength_reversed = _compare_sequences(miRNA, mRNA_substring[::-1])
        else: # debug to file output
            match_strength_straight = _compare_sequences_debug(miRNA, mRNA_substring, i, fileWrite_filepath_p)
            match_strength_reversed = _compare_sequences_debug(miRNA, mRNA_substring[::-1], i, fileWrite_filepath_p)
        
        match_strength = 0
        if(match_strength_straight > match_strength_reversed):
            match_strength = match_strength_straight
        else:
            match_strength = match_strength_reversed
        if match_strength > max_match_strength:
            max_match_strength = match_strength
        i += 1
    
    if debugLog == True:
        logger.debug(f"{_d_miRNA_id} :: {_d_mRNA_id} match score: {max_match_strength}")
        logger.debug(f"miRNA_len = {miRNA_size}, mRNA_len = {mRNA_size} (diff = {mRNA_size - miRNA_size}): opcount_total = {constants._d_total_opcounts/2}") # divide total opcounts by 2, because we do a forward and reverse comparison!

    return max_match_strength

def init_CUDA_sequence_comparison_results_json(cuda_analysis_filepath = "test_run_2/sequence_comparison_results.json"):
    """
    This function should be called to initialise constants.CUDA_sequence_comparison_results_json! 
    """
    if constants.CUDA_sequence_comparison_results_isinit == True:
        return
    if not os.path.exists(cuda_analysis_filepath):
        logger.info(f"ERROR! The filepath {cuda_analysis_filepath} does not exist!")
        return -1
    constants.CUDA_sequence_comparison_results_json = read_file_as_json(cuda_analysis_filepath)
    if(len(constants.CUDA_sequence_comparison_results_json) > 0):
        logger.info("CUDA sequence comparison results json init success.")
        constants.CUDA_sequence_comparison_results_isinit = True

def init_MI_HSA_miRNA_mapping(mi_hsa_miRNA_mapping_filepath = "src_data_files/miRNAdbs/mirbase_miRNA_hsa-only.txt"):
    """
    Populates constants.MI_HSA_miRNA_id_mapping dictionary with MI miRNA ids (as keys) and respective HSA miRNA ids (as values).
    """
    if constants.MI_HSA_miRNA_id_mapping_isinit == True:
        return
    if not os.path.exists(mi_hsa_miRNA_mapping_filepath):
        logger.info(f"ERROR! The filepath {mi_hsa_miRNA_mapping_filepath} does not exist!")
        return -1
    MI_HSA_miRNA_file_readlines = readlines(mi_hsa_miRNA_mapping_filepath)
    for line in MI_HSA_miRNA_file_readlines:
        splitline = line.split("\t")
        MI_id = splitline[0]
        HSA_id = splitline[1]
        constants.MI_HSA_miRNA_id_mapping[MI_id] = HSA_id
    if len(constants.MI_HSA_miRNA_id_mapping) > 0: 
        logger.info(f"Initialised {len(constants.MI_HSA_miRNA_id_mapping)} elements in constants.MI_HSA_miRNA_id_mapping")
        constants.MI_HSA_miRNA_id_mapping_isinit = True
    else:
        logger.info("Error! constants.MI_HSA_miRNA_id_mapping has 0 elements.")

def map_mi_hsa(miRNA_mi_id, miRNA_hsa_id, reverse = 0):
    """
    Finds a hsa-xxx mRNA id from a mi-xxx mRNA id.
    If reverse = 1, then it finds a mi-xxx mRNA id from a hsa-xxx mRNA id

    If you want to find a HSA id from a MI id:
        map_mi_hsa("MIXXXXXX", "", 0)
    
    If you want to find a MI id from a HSA id:
        map_mi_hsa("", "hsa-let-xxx", 1)
    """
    if constants.MI_HSA_miRNA_id_mapping_isinit == False:
        logger.info(f"MI_HSA_miRNA_mappings haven't been initialised. Attempting init.")
        init_MI_HSA_miRNA_mapping()
    if reverse == 0:
        result_hsa_id = constants.MI_HSA_miRNA_id_mapping[miRNA_mi_id]
        if result_hsa_id != None:
            return result_hsa_id
        else:
            logger.debug(f"Error querying {miRNA_mi_id} for hsa id!")
    else: # reverse == 1, should map hsa to mi
        result_mi_id = get_dict_key_at_value(constants.MI_HSA_miRNA_id_mapping, miRNA_hsa_id)
        if  result_mi_id != None:
            return result_mi_id
        else:
            logger.debug(f"Error querying {miRNA_hsa_id} for mi id!")


def find_CUDA_miRNA_mRNA_match_strength(miRNA_id, mRNA_id, cuda_analysis_filepath = "test_run_2/sequence_comparison_results.json"):
    """
    Finds the match strength (from the brute-force CUDA process) between a mRNA and a miRNA.

    @param miRNA_id: Either a MI (miRDB) id (eg. MI0000060) or a HSA id (eg. hsa-let-7f-1)
    @param mRNA_id: A full UniprotKB mRNA id (eg. UniProtKB:Q0VGL1)
    @cuda_analysis_filepath: a filepath to the .json file which was generated using util.load_sequence_comparison_results and subsequently saved using util.save_json

    @return (float) brute-force match strength (based on sequence comparison) between the specified mRNA and miRNA
    """

    if not check_file_exists(cuda_analysis_filepath):
        return -1
    if constants.CUDA_sequence_comparison_results_isinit == False:
        logger.info(f"CUDA sequence comparison results haven't been initialised using util.init_CUDA_sequence_comparison_results_json. Attempting init.")
        init_CUDA_sequence_comparison_results_json()
    if constants.MI_HSA_miRNA_id_mapping_isinit == False:
        logger.info(f"MI_HSA_miRNA_mappings haven't been initialised. Attempting init.")
        init_MI_HSA_miRNA_mapping()
    if constants.mRNA_miRNA_id_sequence_dicts_isinit == False:
        logger.info(f"miRNA-mRNA id to sequence dictionary not initialised. Attempting init.")
        init_mRNA_and_miRNA_sequence_dicts()

    # if a miRNA sequence is supplied, get the ID
    if "hsa" not in miRNA_id and "MI" not in miRNA_id:
        _miRNA_id = get_dict_key_at_value(constants.miRNA_id_sequence_dict, miRNA_id)
        if _miRNA_id == None:
            logger.info(f"ERROR! Couldn't find miRNA id from {miRNA_id}")
        miRNA_id = _miRNA_id
    
    # if a mRNA sequence is supplied, get the ID
    if "Prot" not in mRNA_id:
        _mRNA_id = get_dict_key_at_value(constants.mRNA_id_sequence_dict, mRNA_id)
        if _mRNA_id == None:
            logger.info(f"ERROR! Couldn't find miRNA id from {mRNA_id}")
        mRNA_id = _mRNA_id
    
    # if miRNA_id is HSA, convert it to MI (implications in the querying of elements in constants.CUDA_sequence_comparison_results_json, where keys are saved as MI identifiers)
    if "hsa" in miRNA_id:
        miRNA_id = map_mi_hsa("", miRNA_id, 1)
    
    # main logic -> return the CUDA-computed match strength
    for list_element in constants.CUDA_sequence_comparison_results_json: # each list element is actually a dicionary
        if list_element.get(miRNA_id) != None:
            match_strengths_dict = list_element[miRNA_id][1] # returns the list of mRNAs and their match strength for this specific miRNA, eg {'UniProtKB:Q0VGL1': 0.2125, 'UniProtKB:Q96GR2': 0.2125, ...}
            return match_strengths_dict[mRNA_id]

def init_mRNA_and_miRNA_sequence_dicts(mRNA_sequences_filepath = os.path.join(constants.TARGET_FOLDER, "product_mRNAs_cpp.txt"), miRNA_sequences_filepath = "src_data_files/miRNAdbs/mirbase_miRNA_hsa-only.txt"):
    """
    Uses mRNA_sequences_filepath (should point to test_run_X/product_mRNAs_cpp.txt, generated by util.save_mRNAs_for_cpp) and
    miRNA_sequences_filepath (should point to src_data_files/miRNAdbs/mirbase_miRNA_hsa-only.txt, generated by util.save_mirbase_hsap_miRNAs_for_cpp)
    to populate constants.mRNA_id_sequence_dict and constants.miRNA_id_sequence_dict, with keys as mRNA UniprotKB ids or miRNA MI ids
    and values as respective sequences.
    """
    if constants.mRNA_miRNA_id_sequence_dicts_isinit == True:
        return
    if not check_file_exists(mRNA_sequences_filepath):
        return -1
    if not check_file_exists(miRNA_sequences_filepath):
        return -1
    
    mRNA_readlines = readlines(mRNA_sequences_filepath)
    miRNA_readlines = readlines(miRNA_sequences_filepath)

    for line in mRNA_readlines:
        line = line.split("\t")
        if len(line) > 1: # prevent empty lines from being read
            # logger.debug(f"line = {line}")
            if line[1] != "None": # check if the sequence isnt None (eg. UniProtKB:Q96EY9 has a faulty sequence)
                constants.mRNA_id_sequence_dict[line[0]] = line[1] # line[0] = UniprotKB id, line[1] = sequence
    
    for line in miRNA_readlines:
        line = line.split("\t")
        if len(line) > 1: # prevent empty lines from being read
            constants.miRNA_id_sequence_dict[line[0]] = line[2] # line[0] = MI id, line[1] = HSA id, line[2] = sequnece
    
    if len(constants.mRNA_id_sequence_dict) > 0 and len(constants.miRNA_id_sequence_dict) > 0:
        logger.info(f"Initialised {len(constants.mRNA_id_sequence_dict)} mRNA and {len(constants.miRNA_id_sequence_dict)} miRNA ids and their respective sequences.")
        constants.mRNA_miRNA_id_sequence_dicts_isinit = True
    else:
        logger.info(f"Error initialising mRNA and miRNA id sequence dictionaries. len constants.mRNA_id_sequence_dict = {len(constants.mRNA_id_sequence_dict)}, len constrants.miRNA_id_sequence_dict = {len(constants.miRNA_id_sequence_dict)}")

def find_mRNA_sequence(mRNA_id):
    """
    Finds the mRNA sequence given the mRNA_id (should be a UniprotKB:xxxxx id).
    You should call util.init_mRNA_and_miRNA_sequence_dicts before calling this function.
    """
    if constants.mRNA_miRNA_id_sequence_dicts_isinit == False:
        logger.info(f"mRNA and miRNA id to sequence mappings not initialised. Attempting init.")
        init_mRNA_and_miRNA_sequence_dicts()
    if constants.mRNA_id_sequence_dict.get(mRNA_id) != None:
        return constants.mRNA_id_sequence_dict[mRNA_id]
    else:
        logger.info(f"ERROR! Couldn't find {mRNA_id} in constants.mRNA_id_sequence_dict.")
        return None

def find_miRNA_sequence(miRNA_id):
    """
    Finds the miRNA sequence given the miRNA_id (should be either MI or HSA type id).
    You should call util.init_mRNA_and_miRNA_sequence_dicts before calling this function.
    """
    if constants.mRNA_miRNA_id_sequence_dicts_isinit == False:
        logger.info(f"mRNA and miRNA id to sequence mappings not initialised. Attempting init.")
        init_mRNA_and_miRNA_sequence_dicts()
    
    # convert to MI type id
    if "hsa" in miRNA_id:
        miRNA_id = map_mi_hsa("", miRNA_id, reverse=1)
    
    if constants.miRNA_id_sequence_dict.get(miRNA_id) != None:
        return constants.miRNA_id_sequence_dict[miRNA_id]
    else:
        logger.info(f"ERROR! Couldn't find {miRNA_id} in constants.miRNA_id_sequence_dict.")
        return None

def compare_python_CUDA_miRNA_mRNA_match_strength(miRNA,mRNA, debugLog=True):
    """
    Does a python comparison of the mRNA to miRNA sequence (using util.compare_miRNA_mRNA_match_strength_single)
    against a CUDA-processed mRNA to miRNA match strength (using util.find_CUDA_miRNA_mRNA_match_strength).
    
    @param miRNA: a miRNA id (either HSA or MI type) or a miRNA sequence
    @param mRNA: a mRNA id (of uniprotkb type) or a mRNA sequence
    @param debugLog: if True, prints out both match strengths (one computed via python and one via CUDA)

    @return: python's match strength divided by CUDA's match strength. If both are equal (this is desired behaviour),
    the return value should be 1 (True).
    """
    CUDA_match_strength = find_CUDA_miRNA_mRNA_match_strength(miRNA,mRNA)
    python_match_strength = compare_miRNA_mRNA_match_strength_single(miRNA, mRNA)
    result = python_match_strength / CUDA_match_strength
    if debugLog == True:
        logger.debug(f"Match strength {miRNA} :: {mRNA} = {result}")
    return result


def check_file_exists(filepath, debug_message=""):
    if os.path.exists(filepath):
        return True
    else:
        if debug_message == "":
            logger.info(f"File at path {filepath} does not exist!")
        else:
            logger.info(debug_message)
        return False

def convert_thymine_uracil(sequence, T_to_U = 1):
    """
    Converts thymine (T) to uracil (U), if 'T_to_U' is 1.
    Converts uracil (U) to thymine (T), if 'T_to_U' is 0.

    Returns the converted sequence
    """
    if T_to_U == 1:
        return sequence.replace("T", "U")
    else:
        return sequence.replace("U", "T")

def replace_in_file(src, dst, src_string, replace_string, ignore_comments=True, comment_delimiter="#"):
    """
    Replaces all occurences of 'src_string' with 'replace_string' from 'src' and writes
    the output to 'dst' file.

    If 'ignore_comments' is True, it ignores all lines delimited by the 'comment_delimiter'
    """
    readlines = readlines(src)
    result_lines = []
    for line in readlines:
        if "#" in line:
            result_lines.append(line)
            continue
        res = line.replace(src_string, replace_string)
        result_lines.append(res)
        
    dst_file = open(dst, "w")
    dst_file.writelines(result_lines)
    dst_file.close()

def generate_go_names_in_input_file(dst_file):
    """
    Appends a Gene Ontology process name in the lines where terms are noted.
    File = constants.TARGET_FOLDER/input.txt
    """
    inputs_file_readlines = readlines(os.path.join(constants.TARGET_FOLDER, "input.txt"))
    terms_angio_pos = convert_linear_list_to_dictionary(constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY, flip=True)
    terms_angio_neg = convert_linear_list_to_dictionary(constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY, flip=True)
    terms_angio_general = convert_linear_list_to_dictionary(constants.TERMS_ANGIOGENESIS_GENERAL, flip=True)
    terms_dia_pos = convert_linear_list_to_dictionary(constants.TERMS_DIABETES_POSITIVE_ARRAY, flip=True)
    terms_dia_neg = convert_linear_list_to_dictionary(constants.TERMS_DIABETES_NEGATIVE_ARRAY, flip=True)
    terms_dia_general = convert_linear_list_to_dictionary(constants.TERMS_DIABETES_GENERAL, flip=True)

    # compute dict between all terms and their represetnative names
    all_terms_dict = {**terms_angio_pos, **terms_angio_neg, **terms_angio_general, **terms_dia_pos, **terms_dia_neg, **terms_dia_general}

    result_lines = []
    segment_terms_pointer = 0
    i = 0
    for line in inputs_file_readlines:
        if segment_terms_pointer == 1: # we are now in the GO_terms segment
            current_GO_term = line.split("\t")[0]
            current_GO_label = all_terms_dict[current_GO_term]
            line = line.replace("\n","")
            line = f"{line}\t{current_GO_label}\n"
        if segment_terms_pointer == 0 and "GO_terms" in line:
            segment_terms_pointer = 1
        i+=1

        result_lines.append(line)
    
    writelines(os.path.join(constants.TARGET_FOLDER, "input1.txt"),result_lines)
        

def convert_linear_list_to_dictionary(linear_list, flip=False):
    """
    Converts a linear list of [a1,b1,a2,b2,a3,b3] to dictionary between ax and bx
    
    If flip == True, then final dictionary is between bx (keys) and ax (values)
    """
    i = 0
    temp_pair = []
    result_dict = {}
    for element in linear_list:
        temp_pair.append(element)
        if i % 2 == 1 and i != 0:
            if flip == False:
                result_dict[temp_pair[0]] = temp_pair[1]
            else:
                result_dict[temp_pair[1]] = temp_pair[0]
            temp_pair = []
        i+=1
    return result_dict
        


""" An older and recursive implementation (new is get_uniprotId_from_geneName_new). Would cause me too much pain to delete.
def get_uniprotId_from_geneName(gene_name, recursion=0, prefix="UniProtKB:", trust_genes=True):
    #
    # Retrieves uniprot identifier from a gene symbol/name; e.g. UniProtKB:Q86SQ4, if gene_name=adgrg6
    #
    # Parameters:
    #  - gene_name: A gene name or symbol e.g. ADGRG6
    #  - recursion: An internal function parameter to iterate through an array of UniProtKB ids supplied by the response json
    #  - prefix: A database prefix, is prefixed before the specifi uniprot gene id
    #  - trust_genes: If True, all trusted genes inside genes_trusted.txt will override this function. If false, it
    #    notifies you that the gene was found among trusted genes and if you wish to proceed with the function.

    next_recursive_step = recursion + 1
    if gene_name in constants.TRUSTED_GENES: # trusted genes loaded at program startup
        if trust_genes == True:
            # return the element ahead of the gene_name, which is the previously found uniprot_id
            return "UniProtKB:" + constants.TRUSTED_GENES[constants.TRUSTED_GENES.index(gene_name)+1]
        else:
            trust_genes_user_response = int(input(f"Gene {gene_name} is found among trusted genes. Press 1 to continue with UniProt Id query, or 0 to use the trusted gene."))
            if trust_genes_user_response == 0:
                return "UniProtKB:" + constants.TRUSTED_GENES[constants.TRUSTED_GENES.index(gene_name)+1]

    global _uniprot_identifier_query_result
    uniprot_gene_identifier = ""
    if recursion == 0:
        _uniprot_identifier_query_result = requests.get(f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene_name}+AND+organism_id:9606&format=json&fields=accession,gene_names,organism_name,reviewed")
        if _uniprot_identifier_query_result.text == "{'results': []}":
            # empty response, may be because ortholog was found, but human ortholog has different name than the zebrafish ortholog
            raise Exception(f"No uniprot identifier query result for {gene_name} found.")
        _uniprot_identifier_query_result = _uniprot_identifier_query_result.json()
        logger.debug(_uniprot_identifier_query_result)

    uniprot_gene_identifier = _uniprot_identifier_query_result["results"][recursion]["primaryAccession"]
    results_arr_len = len(_uniprot_identifier_query_result["results"])

    # as per 'uniprot-mapping-multiple-results-issue': implemented functionality to dismiss results that have different
    # genes as the speicifed gene (idk why uniprot returns such genes nevertheless)
    genes_all = []
    genes_synonyms = "" # init to prevent errors
    genes_primary = _uniprot_identifier_query_result["results"][recursion]["genes"][0]["geneName"]["value"]
    if "synonyms" in _uniprot_identifier_query_result["results"][recursion]["genes"][0]:
        genes_synonyms = _uniprot_identifier_query_result["results"][recursion]["genes"][0]["synonyms"][0]["value"]
    logger.debug(f"genes_primary = {genes_primary}, genes_synonyms = {genes_synonyms}")
    if isinstance(genes_primary, list): 
        for g in genes_primary: genes_all.append(g)
    else: genes_all.append(genes_primary)
    if isinstance(genes_synonyms, list):
        for g in genes_synonyms: genes_all.append(g)
    else: genes_all.append(genes_synonyms)
    logger.debug(f"genes_all = {genes_all}")

    if results_arr_len == 1: #if only one result, auto accept it
        logger.info(f"[get_uniprot_identifier]: Auto translated {gene_name} -> {uniprot_gene_identifier}. Reason: Only 1 result.")
        if prefix != "":
            return prefix+uniprot_gene_identifier
        else:
            return uniprot_gene_identifier

    is_reviewed = False if "TrEMBL" in _uniprot_identifier_query_result["results"][recursion]["entryType"] else True

    logger.info(f"Gene name {gene_name} found to correspond to {uniprot_gene_identifier} (Reviewed: {is_reviewed}). Displaying response {recursion+1}/{results_arr_len}: {_get_uniprot_identifier_json_nth_response(_uniprot_identifier_query_result,recursion)}")
    logger.info(get_uniprotId_description(uniprot_gene_identifier))
    # check if current gene name is found among genes_all (final check to mitigate 'uniprot-mapping-multiple-results-issue')
    if gene_name not in genes_all:
        logger.info(f"WARNING: Gene name {gene_name} was not found among gene fields of uniprot query. Fields = {genes_all}. Skipping this result.")
        get_uniprotId_from_geneName(gene_name, recursion=next_recursive_step)

    user_logic = int(input("Press 1 to confirm current result, 2 to cycle another result or 0 to continue the program and discard all options."))
    if user_logic == 1:
        # result is confirmed, save so in TRUSTED_GENES
        logger.info(f"[get_uniprot_identifier]: Confirmed {gene_name} -> {uniprot_gene_identifier}")
        with open("src_data_files/genes_trusted.txt", "a+") as f:
            f.seek(0) # a+ has file read pointer at bottom, f.seek(0) returns to top
            logger.debug(f"Opened genes_trusted file.")
            if gene_name not in f.read():
                f.seek(0,2) # move file pointer back to the end of the file
                f.write(f"{gene_name} {uniprot_gene_identifier}\n")
                f.close()
                logger.debug(f"Writing to genes_trusted file done!")
    elif user_logic == 2:
        # cycle another result
        if next_recursive_step > (results_arr_len-1): 
            # cycled through all options, return error
            logger.info("Cycled out of options!")
            return f"[get_uniprot_identifier_new]: CycleOutOfBoundsError: Cycled through all the uniprot gene IDs for {gene_name} without confirming any"
            # todo: major bug here, it should return error, but it returns a gene id instead
            # look at: https://stackoverflow.com/questions/11356168/return-in-recursive-function
            # and https://stackoverflow.com/questions/23543485/python-recursive-function-executing-return-statement-incorrectly 
            # todo: handle this error in main script file
        get_uniprotId_from_geneName(gene_name, recursion=next_recursive_step)
    elif user_logic == 0:
        return f"[get_uniprot_identifier_new]: No uniprot gene IDs for {gene_name} found."
    else:
        # wrong input, recurse again
        logger.debug("[get_uniprot_identifier_new]: Wrong input! Must be 0, 1 or 2.")
        get_uniprotId_from_geneName(gene_name, recursion)

    if prefix != "":
        return prefix+uniprot_gene_identifier
    else:
        return uniprot_gene_identifier
"""






