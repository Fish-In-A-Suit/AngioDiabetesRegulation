# Utility functions file
import requests
import constants
import json
import os

import logging
logger = logging.getLogger(__name__)

import sys

_response_cycle_counter = 0
_uniprot_identifier_query_result = ""

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
    if array_name == 'ALL':
        if term_shrink: return shrink_term_list(constants.TERMS_ANGIOGENESIS_GENERAL + constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY + constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY + constants.TERMS_DIABETES_GENERAL + constants.TERMS_DIABETES_NEGATIVE_ARRAY + constants.TERMS_DIABETES_POSITIVE_ARRAY)
        else: return constants.TERMS_ANGIOGENESIS_GENERAL + constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY + constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY + constants.TERMS_DIABETES_GENERAL + constants.TERMS_DIABETES_NEGATIVE_ARRAY + constants.TERMS_DIABETES_POSITIVE_ARRAY
    elif array_name == 'ANGIOGENESIS':
        if term_shrink: return shrink_term_list(constants.TERMS_ANGIOGENESIS_GENERAL + constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY + constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY)
        else: return constants.TERMS_ANGIOGENESIS_GENERAL + constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY + constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY
    elif array_name == 'ANGIOGENESIS-NEGATIVE':
        if term_shrink: return shrink_term_list(constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY)
        else: return constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY
    elif array_name == 'ANGIOGENESIS-POSITIVE':
        if term_shrink: return shrink_term_list(constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY)
        else: return constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY
    elif array_name == 'ANGIOGENESIS-GENERAL':
        if term_shrink: return shrink_term_list(constants.TERMS_ANGIOGENESIS_GENERAL)
        else: return constants.TERMS_ANGIOGENESIS_GENERAL
    elif array_name == 'DIABETES':
        if term_shrink: return shrink_term_list(constants.TERMS_DIABETES_GENERAL + constants.TERMS_DIABETES_NEGATIVE_ARRAY + constants.TERMS_DIABETES_POSITIVE_ARRAY)
        else: return constants.TERMS_DIABETES_GENERAL + constants.TERMS_DIABETES_NEGATIVE_ARRAY + constants.TERMS_DIABETES_POSITIVE_ARRAY
    elif array_name == 'DIABETES-NEGATIVE':
        if term_shrink: return shrink_term_list(constants.TERMS_DIABETES_NEGATIVE_ARRAY)
        else: return constants.TERMS_DIABETES_NEGATIVE_ARRAY
    elif array_name == 'DIABETES-POSITIVE':
        if term_shrink: return shrink_term_list(constants.TERMS_DIABETES_POSITIVE_ARRAY)
        else: return constants.TERMS_DIABETES_POSITIVE_ARRAY
    elif array_name == 'DIABETES-GENERAL':
        if term_shrink: return shrink_term_list(constants.TERMS_DIABETES_GENERAL)
        else: return constants.TERMS_DIABETES_GENERAL
    else:
        print(array_name + " could not be found! Returning empty array.")
        empty = []
        return empty

def read_file_as_json(filepath):
    """
    Reads the file into json
    """
    with open(filepath, "r") as read_content:
        return json.load(read_content)

def store_json_dictionaries(filepath, dictionaries):
    """
    Writes the json dictionaries to file at filepath
    """
    file = open(filepath, "w+")
    file.write(json.dumps(dictionaries)+"\n")
    file.close()

def readlines(filepath):
    """
    Reads the lines of the specified filepath
    """
    with open(filepath, "r") as read_content:
        return read_content.readlines()

def get_last_geneId_in_crash_json(json):
    """
    Gets the Id of the last gene in supplied json. Used in the crash recovery algorithm
    """
    logger.debug(f"[get_last_geneId_in_crash_json]: json_len = {len(json)}, json = {json}")
    geneId = json[len(json)-1]["product"] # get last item in array (len(json[0])-1) and then "product", which is geneId
    logger.debug(f"[get_last_geneId_in_crash_json]: geneId = {geneId}")
    return geneId

def shrink_term_list(list):
    i = 0
    result_list = []
    for element in list:
        if i%2: # 0 = generic term name, 1 = GO:xxxx, 2 = generic term name, 3 = GO:xxxx -> modulo op is 1 at 1, 3, 5 etc
            result_list.append(element)
        i = i+1
    return result_list

def list_directionshrink(list, reference_element, forward=True):
    """
    If forward = True, keeps all elements in list starting after the reference_element
    If forward = False, keeps all elements in list starting before the reference_element

    Example: 
    ls = [0,1,2,3,4,5,6]
    list_directionshrink(ls, 2, forward=True) --> [3,4,5,6]
    """
    result_list = []
    start_append_forward = False
    stop_append_backward = False
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


def zfin_find_human_ortholog(gene_id, ortholog_file_path="src_data_files/zfin_human_ortholog_mapping.txt"):
    """
    If gene_id is from the ZFIN database, searches through the zebrafish-human orthologs and returns the name of the
    symbol of the human gene ortholog.
    """
    #file = open(ortholog_file_path, "r") # TODO: make ortholog files global, init at runtime
    #lines = file.readlines()
    gene_id=gene_id.split(":")[1] # eliminates ZFIN: 
    for line in _zfin_ortholog_readlines:
        if gene_id in line:
            human_symbol = _zfin_get_human_gene_symbol_from_line(line)
            logger.info(f"[zfin_find_human_ortholog]: Returning human symbol {human_symbol}")
            #file.close()
            return human_symbol
    #file.close()
    return f"[ZfinError_No-human-ortholog-found:gene_id={gene_id}"

def _zfin_get_human_gene_symbol_from_line(line, improved_algorithm=True):
    """
    Splits zfin line and retrieves human gene symbol (full caps of zebrafish gene symbol)
    """
    if improved_algorithm == True:
        # better, as zfin human ortholog sometimes has different name than the zebrafish gene
        return line.split("\t")[3] # look at zfin orthologs txt file (in src_data_files) -> when you higlight a row, you see a TAB as '->' and a SPACEBAR as '.' -> splitting at \t means every [3] linesplit element is the human gene name
    else: 
        return str(line.split("\t")[1]).upper() # split lines at tabs (spacebar is not ok!)

def get_uniprotId_from_geneName_new(gene_name, prefix="UniProtKB:", trust_genes=True):
    """
    Retrieves UniProt Identifier from a gene symbol/name; e.g. UniProtKB:Q86SQ4, if gene_name=adgrg6. 
    This function no longer uses recursion and also checks for verification status of the returned Ids.

    Parameters:
      - gene_name: A gene name or symbol e.g. ADGRG6
      - prefix: A database prefix, is prefixed before the specifi uniprot gene id
      - trust_genes: If True, all trusted genes inside genes_trusted.txt will override this function. If false, it
        notifies you that the gene was found among trusted genes and if you wish to proceed with the function.
    """
    # check if gene exists in trusted genes
    if gene_name in constants.TRUSTED_GENES:
        if trust_genes == True:
            # return the element ahead of the gene_name, which is the previously found uniprot_id
            return "UniProtKB:" + constants.TRUSTED_GENES[constants.TRUSTED_GENES.index(gene_name)+1]
        else:
            trust_genes_user_response = int(input(f"Gene {gene_name} is found among trusted genes. Press 1 to continue with UniProt Id query, or 0 to use the trusted gene."))
            if trust_genes_user_response == 0:
                return "UniProtKB:" + constants.TRUSTED_GENES[constants.TRUSTED_GENES.index(gene_name)+1]
    
    # load all of the genes
    _uniprot_identifier_query_result = requests.get(f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene_name}+AND+organism_id:9606&format=json&fields=accession,gene_names,organism_name,reviewed")
    if _uniprot_identifier_query_result.text == "{'results': []}":
        # empty response, may be because ortholog was found, but human ortholog has different name than the gene from the file - denotes an error in file parsing
        raise Exception(f"No uniprot identifier query result for {gene_name} found.")
    _uniprot_identifier_query_result = _uniprot_identifier_query_result.json()

    uniprot_geneIds_dictionary = {} # stores uniprotId : reviewed_status pairs
    results_arr_len = len(_uniprot_identifier_query_result["results"])

    # if only one result, auto accept it; TODO: check if this 1 result is verified; if unverified, give user option what to do && ALSO CHECK IF GENE_NAME IS IN GENENAMES REQUEST FIELD, possibly move this down
    if results_arr_len == 1:
        uniprotId = _uniprot_identifier_query_result["results"][0]["primaryAccession"]
        logger.info(f"[get_uniprot_identifier]: Auto accepted {gene_name} -> {uniprotId}. Reason: Only 1 result.")
        if prefix != "": return prefix + uniprotId
        else: return uniprotId

    for i in range(results_arr_len):
        # reviewed is either TrEMBL or SwissProt. TrEMBL should be set to False and SwissProt reviews should be set to True
        is_reviewed = False if "TrEMBL" in _uniprot_identifier_query_result["results"][i]["entryType"] else True
        geneId = _uniprot_identifier_query_result["results"][i]["primaryAccession"]
        uniprot_geneIds_dictionary[geneId] = is_reviewed # add to dict
    logger.info(f"[get_uniprotId_from_geneName_new]: uniprot_geneIds_dictionary = {uniprot_geneIds_dictionary}")

    # Used to check the gene_names field of the request against the gene_name parameter supplied by the function. If gene_name (parameter) is not
    # among gene_names, then such element should not be analysed further.
    # "As per 'uniprot-mapping-multiple-results-issue': implemented functionality to dismiss results that have different genes as the specified gene (idk why uniprot returns such genes nevertheless)"
    geneIds_impactGenes_dictionary = {} # stores <str> uniprotId : <list> impact_genes; used later to check that uniprotId contains relevant gene
    for j in range(results_arr_len):
        impactGenes_all = []
        impactGenes_synonyms = "" # init to prevent errors
        impactGenes_primary = _uniprot_identifier_query_result["results"][j]["genes"][0]["geneName"]["value"]
        if "synonyms" in _uniprot_identifier_query_result["results"][j]["genes"][0]:
            impactGenes_synonyms = _uniprot_identifier_query_result["results"][j]["genes"][0]["synonyms"][0]["value"]
        if isinstance(impactGenes_primary, list):
            for g in impactGenes_primary: impactGenes_all.append(g)
        else: impactGenes_all.append(impactGenes_primary)
        if isinstance(impactGenes_synonyms, list):
            for g in impactGenes_synonyms: impactGenes_all.append(g)
        else: impactGenes_all.append(impactGenes_synonyms)
        geneIds_impactGenes_dictionary[get_dict_key_at_index(uniprot_geneIds_dictionary,j)] = impactGenes_all
    logger.debug(f"[get_uniprotId_from_geneName_new]: geneIds_impactGenes_dictionary = {geneIds_impactGenes_dictionary}")

    # Eliminate entries where gene_name is not among "impact genes" - final check to mitigate 'uniprot-mapping-multiple-results-issue'
    # Also checks for autoaccept (if onyl one of the UniprotIds is reviewed) to prevent multiple for loop reiterations
    pop_keys = [] # indices of items to eliminate
    pop_keys_indices = []
    reviewedId_single = ""
    NO_reviewed_Ids = 0
    for key, value in geneIds_impactGenes_dictionary.items():
        if gene_name not in value: 
            pop_keys.append(key)
            pop_keys_indices.append(i)
        else: # gene_name is found among geneNames request return field -
            if uniprot_geneIds_dictionary.get(key) == True: # check if UniprotId is reviewed
                if reviewedId_single == "": # check if no other reviewed id was added before
                    reviewedId_single = key
                    NO_reviewed_Ids += 1
    # eliminate faulty elements    
    for pop_key in pop_keys:
        uniprot_geneIds_dictionary.pop(pop_key)
        geneIds_impactGenes_dictionary.pop(pop_key)
    logger.debug(f"[get_uniprotId_from_geneName_new]: removed {len(pop_keys)} faulty ids (gene_name did not match with geneNames field). New uniprot_geneIds_dictionary len = {len(uniprot_geneIds_dictionary)}")

    # Autoaccept if only one of the UniprotIds is reviewed
    if NO_reviewed_Ids == 1:
        logger.info(f"[get_uniprotId_from_geneName_new]: Found a single reviewed UniProt Id for gene_name {gene_name}: {reviewedId_single}")
        return reviewedId_single
    
    # Either 2+ or 0 verified UniprotIds -> begin user cycling options
    results_arr_len = len(uniprot_geneIds_dictionary) # update after eliminations
    i = 0
    next_step = 1
    for uprId, isReviewed in uniprot_geneIds_dictionary: # used camelcase for differentiating inter-for-loop variables (from function variables)
        logger.info(f"Gene name {gene_name} found to correspond to {uprId} (reviewed = {isReviewed}). Displaying response {i}/{results_arr_len}: {_get_uniprot_identifier_json_nth_response(_uniprot_identifier_query_result,i)}")
        logger.info(get_uniprotId_description(uprId))
        
        user_logic = int(input("Press 1 to confirm current result, 2 to cycle another result, or 0 to continue the program and discard all options."))
        if user_logic == 1: # result is confirmed, save so in TRUSTED_GENES
            logger.info(f"[get_uniprotId_from_geneName_new]: Confirmed {gene_name} -> {uprId}")
            with open("src_data_files/genes_trusted.txt", "a+") as f:
                f.seek(0) # a+ has file read pointer at bottom, f.seek(0) returns to top
                logger.debug(f"Opened genes_trusted file.")
                if gene_name not in f.read():
                    f.seek(0,2) # move file pointer back to the end of the file
                    f.write(f"{gene_name} {uprId}\n")
                    f.close()
                    logger.debug(f"Writing to genes_trusted file done!")
            # return
            # TODO: handle a case if a user wants to confirm multiple reviewed proteins -> store them in a list and return when next_step > (results_arr_len-1)
            if prefix != "": return prefix + uprId
            else: return uprId
        elif user_logic == 2: #cycle next result
            if next_step > (results_arr_len - 1):
                logger.info("Cycled out of options!")
                return f"[get_uniprot_identifier_new]: CycleOutOfBoundsError: Cycled through all the uniprot gene IDs for {gene_name} without confirming any"
            i += 1 
            next_step += 1
            continue # skip to the next for loop element
        elif user_logic == 0: #TODO: debate the use of this, add return / redo / start loop over functionality
            logger.info("No uniprot geneIds selected.")
            return f"[get_uniprot_identifier_new]: No uniprot gene Ids for {gene_name} selected."
        else:
            # wrong input
            raise Exception(f"[get_uniprot_identifier_new]: Wrong input!")
    logger.info("No uniprot geneIds selected")
    return f"[get_uniprot_identifier_new]: No uniprot gene Ids for {gene_name} selected."

def get_uniprotId_from_geneName(gene_name, recursion=0, prefix="UniProtKB:", trust_genes=True):
    """
    Retrieves uniprot identifier from a gene symbol/name; e.g. UniProtKB:Q86SQ4, if gene_name=adgrg6

    Parameters:
      - gene_name: A gene name or symbol e.g. ADGRG6
      - recursion: An internal function parameter to iterate through an array of UniProtKB ids supplied by the response json
      - prefix: A database prefix, is prefixed before the specifi uniprot gene id
      - trust_genes: If True, all trusted genes inside genes_trusted.txt will override this function. If false, it
        notifies you that the gene was found among trusted genes and if you wish to proceed with the function.
    """
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

    # TODO:Reviewed proteins have priority. If only one in the list is reviewed, auto accept it.
    # reviewed either SwissProt or Trembl -> if Trembl, then false
    # poglej, če je samo 1 reviewed gen v setu genov -> če ja: autoselect | če ne: ročno select
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
            # TODO: major bug here, it should return error, but it returns a gene id instead
            # look at: https://stackoverflow.com/questions/11356168/return-in-recursive-function
            # and https://stackoverflow.com/questions/23543485/python-recursive-function-executing-return-statement-incorrectly 
            # TODO: handle this error in main script file
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

def _get_uniprot_identifier_json_nth_response(json, nth_response):
    "Gets the entire nth element in 'results' array inside the json retrieved by get_uniprot_identifier function"
    return json["results"][nth_response]

def load_trusted_genes(trusted_genes_file_path):
    """
    Loads constants.TRUSTED_GENES list with genes from genes_trusted.txt
    """
    file = open(trusted_genes_file_path, "r")
    lines = file.readlines()
    for line in lines:
        splitlist = line.split(" ")
        for element in splitlist:
            constants.TRUSTED_GENES.append(element)

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

def xenbase_find_human_ortholog(gene_id, ortholog_file_path="src_data_files/xenbase_human_ortholog_mapping.txt"):
    """
    Attempts to find a human ortholog from the xenbase database.
    Parameters:
      - gene_id: eg. Xenbase:XB-GENE-495335 or XB-GENE-495335
    Returns: symbol of the human ortholog gene (eg. rsu1) or 'XenbaseError_no-human-ortholog-found'
    """
    # file = open(ortholog_file_path, "r") # TODO: make ortholog files global, init at runtime
    # lines = file.readlines()
    gene_id_short = ""
    if ":" in gene_id: gene_id_short = gene_id.split(":")[1]
    else: gene_id_short = gene_id
    
    for line in _xenbase_ortholog_readlines:
        if gene_id_short in line:
            human_symbol = _xenbase_get_human_symbol_from_line(line)
            logger.info(f"Found human ortholog {human_symbol} for xenbase gene {gene_id}")
            # file.close()
            return human_symbol
    # file.close()
    return f"[XenbaseError_No-human-ortholog-found:gene_id={gene_id}"

def _xenbase_get_human_symbol_from_line(line):
    """
    Splits xenbase line at tabs and gets human gene symbol (in full caps)
    """
    return str(line.split("\t")[2]).upper()

def mgi_find_human_ortholog(gene_id):
    """
    Attempts to find a human ortholog from the mgi database.
    Parameters: gene-id eg. MGI:MGI:98480
    Returns: symbol of the human ortholog gene or "MgiError_no-human-ortholog-found".
    """
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
            human_symbol = _mgi_get_human_symbol_from_line(_mgi_ortholog_readlines[i+1])
            logger.info(f"Found human ortholog {human_symbol} for mgi gene {gene_id}")
            return human_symbol # return here doesnt affect line counter 'i', since if gene is found i is no longer needed
        i += 1
    return f"[MgiError_No-human-ortholog-found:gene_id={gene_id}"

def _mgi_get_human_symbol_from_line(line):
    """
    Splits mgi line at tabs and gets human gene symbol
    """
    split = line.split("\t")
    if split[1] != "human":
        raise Exception(f"MGI line {line} doesn't contain keyword 'human'!")
    return split[3]

def rgd_find_human_ortholog(gene_id):
    """
    Attempts to find a human ortholog from the RGD (rat genome database)
    """
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

def _rgd_get_human_symbol_from_line(line):
    """
    Splits rgd line at tabs and gets human gene smybol
    """
    # also clears whitespace from linesplit (which is split at tab). Some lines in RGD db text file had whitespace instead of \t -> clear whitespace from array to resolve
    # example: linesplit = ['Ang2', '1359373', '497229', '', '', '', '', 'Ang2', '1624110', '11731', 'MGI:104984', 'RGD', '\n']
    linesplit = line.split("\t")
    result_list = [] 
    for element in linesplit: 
        if element != "":
            result_list.append(element)
    return result_list[3]

def sort_list_of_dictionaries(input, field, direction_reversed = True):
    """Sorts the list of dictionaries by the key "field", default direction is reversed (descending)"""
    return sorted(input, key=lambda d: d[field], reverse=direction_reversed)

def append_to_file(srcfilepath, append):
    """
    Appends the 'append' to the file before the filetype.
    Example usage: util.append_to_file("term_genes/GO-0001525.json", ";params=homosapiens_only,v1")
    --> result: term_genes/GO-0001525;params=homosapiens_only,v1.json
    """
    split=srcfilepath.split(".")
    dstfilepath = split[0] + append + "." + split[1]
    os.rename(srcfilepath, dstfilepath)

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

def get_dict_key_at_index(dictionary, index):
    """
    Returns the value of the dictionary key at specified index.
    """
    keys_list = list(dictionary) # call list(dict) on a dictionary to return a list of its keys
    key_at_index = keys_list[index]
    return key_at_index




