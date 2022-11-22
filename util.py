# Utility functions file
import requests
import constants
import json
import os
import collections
import shutil

import logging
logger = logging.getLogger(__name__)

import sys

_response_cycle_counter = 0
_uniprot_query_result = ""

_zfin_ortholog_readlines = ""
_xenbase_ortholog_readlines = ""
_mgi_ortholog_readlines = ""
_rgd_ortholog_readlines = ""


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
        json.dump(json_object, f)
    logger.info(f"Stored json to file {filepath}")

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

def store_json_dictionaries(filepath, dictionaries):
    """
    Writes the json dictionaries to file at filepath
    """
    if json.dumps(dictionaries) != "[]":
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        file = open(filepath, "w+")
        file.write(json.dumps(dictionaries)+"\n")
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
    keys_list = list(dictionary) # call list(dict) on a dictionary to return a list of its keys
    key_at_index = keys_list[index]
    return key_at_index

def _return_ensembl_from_id_and_uniprot_query(uniprotId, query):
    logger.debug(f"Starting retrival of ensemblId for uniprotId {uniprotId}")
    index = next((index for (index, d) in enumerate(query["results"]) if d["primaryAccession"] == uniprotId), None)
    ensembl_index_list=[]
    xref_arr_length = len(query["results"][index]["uniProtKBCrossReferences"])
    for i in range(xref_arr_length):
        if query["results"][index]["uniProtKBCrossReferences"][i]["database"] == "Ensembl":
            ensembl_index_list.append(i)

    if len(ensembl_index_list) == 0:
        enId = None
    elif len(ensembl_index_list) == 1:
        enId=query["results"][index]["uniProtKBCrossReferences"][ensembl_index_list[0]]["id"].split(".")[0]
    elif len(ensembl_index_list) > 1:
        if any("idoformId" in query["results"][index]["uniProtKBCrossReferences"][i] for i in ensembl_index_list):
            for i in ensembl_index_list:
                if "-1" in query["results"][index]["uniProtKBCrossReferences"][i]["idoformId"]:
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
            _uniprot_identifier_query_result = requests.get(f"https://rest.uniprot.org/uniprotkb/search?query={gene_name}+AND+organism_id:9606&format=json&fields=accession,gene_names,organism_name,reviewed,xref_ensembl")
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
        return None, None
    
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
            raise Exception(f"[get_uniprot_identifier_new]: Wrong input!")
    logger.info("No uniprot geneIds selected")
    return f"No uniprot gene Ids for {gene_name} selected."

def zfin_find_human_ortholog(gene_id, ortholog_file_path="src_data_files/zfin_human_ortholog_mapping.txt"):
    """
    If gene_id is from the ZFIN database, searches through the zebrafish-human orthologs and returns the name of the
    symbol of the human gene ortholog.
    """
    def _zfin_get_human_gene_symbol_from_line(line, improved_algorithm=True):
        """
        Splits zfin line and retrieves human gene symbol (full caps of zebrafish gene symbol)
        """
        if improved_algorithm == True:
            # better, as zfin human ortholog sometimes has different name than the zebrafish gene
            return line.split("\t")[3] # look at zfin orthologs txt file (in src_data_files) -> when you higlight a row, you see a TAB as '->' and a SPACEBAR as '.' -> splitting at \t means every [3] linesplit element is the human gene name
        else: 
            return str(line.split("\t")[1]).upper() # split lines at tabs (spacebar is not ok!)

    gene_id=gene_id.split(":")[1] # eliminates 'ZFIN:' 
    for line in _zfin_ortholog_readlines:
        if gene_id in line:
            human_symbol = _zfin_get_human_gene_symbol_from_line(line)
            logger.info(f"[zfin_find_human_ortholog]: Returning human symbol {human_symbol}")
            return human_symbol
    return f"[ZfinError_No-human-ortholog-found:gene_id={gene_id}"

def xenbase_find_human_ortholog(gene_id, ortholog_file_path="src_data_files/xenbase_human_ortholog_mapping.txt"):
    """
    Attempts to find a human ortholog from the xenbase database.
    Parameters:
      - gene_id: eg. Xenbase:XB-GENE-495335 or XB-GENE-495335
    Returns: symbol of the human ortholog gene (eg. rsu1) or 'XenbaseError_no-human-ortholog-found'
    """
    def _xenbase_get_human_symbol_from_line(line):
        """Splits xenbase line at tabs and gets human gene symbol (in full caps)"""
        return str(line.split("\t")[2]).upper()

    gene_id_short = ""
    if ":" in gene_id: gene_id_short = gene_id.split(":")[1]
    else: gene_id_short = gene_id
    
    for line in _xenbase_ortholog_readlines:
        if gene_id_short in line:
            human_symbol = _xenbase_get_human_symbol_from_line(line)
            logger.info(f"Found human ortholog {human_symbol} for xenbase gene {gene_id}")
            return human_symbol
    return f"[XenbaseError_No-human-ortholog-found:gene_id={gene_id}"

def mgi_find_human_ortholog(gene_id):
    """
    Attempts to find a human ortholog from the mgi database.
    Parameters: gene-id eg. MGI:MGI:98480
    Returns: symbol of the human ortholog gene or "MgiError_no-human-ortholog-found".
    """
    def _mgi_get_human_symbol_from_line(line, line_index):
        """
        Splits mgi line at tabs and gets human gene symbol
        """
        split = line.split("\t")
        if split[1] != "human":
            # try i+2 to check one line further down
            line = _mgi_ortholog_readlines[line_index+2]
            split = line.split("\t")
            if split[1] == "human":
                logger.debug(f"Found keyword 'human' on secondpass line querying.")
                return split[3]
            else:
                # this still means no human ortholog!
                # example: MGI:2660935 (Prl3d2) contains no "human" (neither i+1 nor i+2), also checked uniprot and no human gene for prl3d2 exists
                return f"[MgiError_No-human-ortholog-found:gene_id={gene_id}"
        return split[3]

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
    """ Attempts to find a human ortholog from the RGD (rat genome database) """
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




