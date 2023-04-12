import json
import requests
import time
import os
from typing import List, Dict, Set, Optional
from tqdm import trange, tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
from tabulate import tabulate, SEPARATING_LINE
import gzip
import time
from datetime import timedelta

import logging
# Initialize logger
logger = logging.getLogger(__name__)


class GOApi:
    def __init__(self):
        self.base_url = "http://api.geneontology.org/api/"
        self.timeout = 5
        self.max_retries = 3

    def fetch_term_data(self, term_id):
        """
        Fetches term data for a given term ID from the Gene Ontology API
        """
        url = f"{self.base_url}ontology/term/{term_id}"
        params = {}
        retries = self.max_retries
        while retries > 0:
            try:
                response = requests.get(url, params=params, timeout=self.timeout)
                if response.ok:
                    data = response.json()
                    return data
                else:
                    logger.warning(f"Error: {response.status_code} - {response.reason}")
            except requests.exceptions.Timeout:
                logger.warning("Timeout error. Retrying...")
            except requests.exceptions.RequestException as e:
                logger.warning(f"Error: {e}")
            retries -= 1
            time.sleep(2)
        return None

    def fetch_term_products(self, term_id):
        """
        Fetches product IDs associated with a given term ID from the Gene Ontology API
        """
        APPROVED_DATABASES = [["UniProtKB", ["NCBITaxon:9606"]],
                      ["ZFIN", ["NCBITaxon:7955"]],
                      #["RNAcentral", ["NCBITaxon:9606"]],
                      ["Xenbase", ["NCBITaxon:8364"]],
                      ["MGI", ["NCBITaxon:10090"]],
                      ["RGD", ["NCBITaxon:10116"]]]
        url = f"{self.base_url}bioentity/function/{term_id}/genes"
        params = {"rows": 10000000}
        retries = self.max_retries
        products_set = set()
        while retries > 0:
            try:
                response = requests.get(url, params=params, timeout=self.timeout)
                response.raise_for_status()
                json = response.json()
                for assoc in json['associations']:
                    if assoc['object']['id'] == term_id and any((database[0] in assoc['subject']['id'] and any(taxon in assoc['subject']['taxon']['id'] for taxon in database[1])) for database in APPROVED_DATABASES):
                        product_id = assoc['subject']['id']
                        products_set.add(product_id)
                products = list(products_set)
                logger.info(f"Fetched products for GO term {term_id}")
                return products
            except (requests.exceptions.RequestException, ValueError):
                retries -= 1
                time.sleep(2)
        return None

class UniProtAPI:
    def __init__(self):
        self.base_url = "https://rest.uniprot.org/uniprotkb/"
    
    def get_uniprot_id(self, gene_name, retries=3, timeout=5):
        """
        Given a gene name, returns the corresponding UniProt ID using the UniProt API.

        Parameters:
        - gene_name (str): name of the gene to search for.
        - retries (int): maximum number of times to retry the request in case of network errors.
        - timeout (int): timeout in seconds for the request.

        Returns:
        - str: UniProt ID if found, None otherwise.
        """

        # Define the URL to query the UniProt API
        url = f"{self.base_url}search?query=gene:{gene_name}+AND+organism_id:9606&format=json&fields=accession,gene_names,organism_name,reviewed,xref_ensembl"

        # Try the request up to `retries` times
        for i in range(retries):
            try:
                # Make the request and raise an exception if the response status is not 200 OK
                response = requests.get(url, timeout=timeout)
                response.raise_for_status()

                # Parse the response JSON and get the list of results
                results = response.json()["results"]

                # If no results were found, return None
                if len(results) == 0:
                    return None

                # If only one result was found, accept it automatically
                elif len(results) == 1:
                    uniprot_id = results[0]["primaryAccession"]
                    logger.info(f"Auto accepted {gene_name} -> {uniprot_id}. Reason: Only 1 result.")
                    return "UniProtKB:" + uniprot_id

                # If multiple results were found, filter out the non-reviewed ones
                reviewed_ids = []
                for result in results:
                    # Skip the result if the gene name is not a match
                    if gene_name not in result["genes"][0]["geneName"]["value"]:
                        continue
                    # Skip the result if it is not reviewed
                    if "TrEMBL" not in result["entryType"]:
                        reviewed_ids.append(result)

                # If no reviewed result was found, return None
                if len(reviewed_ids) == 0:
                    return None

                # If only one reviewed result was found, accept it automatically
                elif len(reviewed_ids) == 1:
                    uniprot_id = reviewed_ids[0]["primaryAccession"]
                    logger.info(f"Auto accepted {gene_name} -> {uniprot_id}. Reason: Only 1 reviewed result.")
                    return "UniProtKB:" + uniprot_id

                # If multiple reviewed results were found, ask the user to choose one
                logger.info(f"Multiple reviewed results found for {gene_name}. Please choose the correct UniProt ID from the following list:")
                for i, result in enumerate(reviewed_ids):
                    genes = result["genes"]
                    impact_genes = set()
                    for gene in genes:
                        impact_genes.add(gene["geneName"]["value"])
                        if "synonyms" in gene:
                            for synonym in gene["synonyms"]:
                                impact_genes.add(synonym["value"])
                    print(f"{i + 1}. {result['primaryAccession']} ({', '.join(impact_genes)})")
                # Get the user's choice and return the corresponding UniProt ID
                # choice = input("> ")  # prompt the user for input, but commented out for now
                choice = "1"  # for testing purposes, use "1" as the user's choice
                if choice.isdigit() and 1 <= int(choice) <= len(reviewed_ids):  # check if the user's choice is valid
                    # get the UniProt ID of the chosen result and return it
                    uniprot_id = reviewed_ids[int(choice) - 1]["primaryAccession"]
                    return "UniProtKB:" + uniprot_id
                else:
                    # raise an error if the user's choice is not valid
                    raise ValueError(f"Invalid choice: {choice}")
            except requests.exceptions.RequestException:
                # if there was an error with the HTTP request, log a warning
                logger.warning(f"Failed to fetch UniProt data for {gene_name}")
            # wait for the specified timeout before retrying the request
            time.sleep(timeout)
            # if all retries fail, return None
        return None
  
    def get_uniprot_info(self, uniprot_id, retries=3, timeout=5):
        """
        Given a UniProt ID, returns a dictionary containing various information about the corresponding protein using the UniProt API.
        """
        def _return_ensembl_from_id_and_uniprot_query(uniprotId, query):
            logger.debug(f"Starting retrival of ensemblId for uniprotId {uniprotId}")
            index = next((index for (index, d) in enumerate(query) if d["primaryAccession"] == uniprotId), None)
            ensembl_index_list=[]
            xref_arr_length = len(query[index]["uniProtKBCrossReferences"]) # the count of cross-referenced databases
            for i in range(xref_arr_length):
                if query[index]["uniProtKBCrossReferences"][i]["database"] == "Ensembl":
                    ensembl_index_list.append(i)

            if len(ensembl_index_list) == 0:
                enId = None
            elif len(ensembl_index_list) == 1:
                enId=query[index]["uniProtKBCrossReferences"][ensembl_index_list[0]]["id"].split(".")[0]
            elif len(ensembl_index_list) > 1:
                if any("isoformId" in query[index]["uniProtKBCrossReferences"][i].keys() for i in ensembl_index_list):
                    for i in ensembl_index_list:
                        if "-1" in query[index]["uniProtKBCrossReferences"][i].get("isoformId", ""):
                            enId=query[index]["uniProtKBCrossReferences"][i]["id"].split(".")[0]
                try: enId
                except NameError:
                    enId=query[index]["uniProtKBCrossReferences"][ensembl_index_list[0]]["id"].split(".")[0] 
                        
            logger.info(f"uniprotId {uniprotId} -> ensemblId {enId}")
            return enId
    
        def _return_refseqnt_from_id_and_uniprot_query(uniprotId, query):
            logger.debug(f"Starting retrival of refseqNT for uniprotId {uniprotId}")
            index = next((index for (index, d) in enumerate(query) if d["primaryAccession"] == uniprotId), None)
            refseqnt_index_list=[]
            xref_arr_length = len(query[index]["uniProtKBCrossReferences"]) # the count of cross-referenced databases
            for i in range(xref_arr_length):
                if query[index]["uniProtKBCrossReferences"][i]["database"] == "RefSeq":
                    refseqnt_index_list.append(i)

            if len(refseqnt_index_list) == 0:
                rsntId = None
            elif len(refseqnt_index_list) == 1:
                if "properties" in query[index]["uniProtKBCrossReferences"][refseqnt_index_list[0]]:
                    rsntId=query[index]["uniProtKBCrossReferences"][refseqnt_index_list[0]]["properties"][0]["value"]
            elif len(refseqnt_index_list) > 1:
                if any("isoformId" in query[index]["uniProtKBCrossReferences"][i].keys() for i in refseqnt_index_list):
                    for i in refseqnt_index_list:
                        if "-1" in query[index]["uniProtKBCrossReferences"][i].get("isoformId", ""):
                            rsntId=query[index]["uniProtKBCrossReferences"][i]["properties"][0]["value"]
                try: rsntId
                except NameError:
                    if "properties" in query[index]["uniProtKBCrossReferences"][refseqnt_index_list[0]]:
                        rsntId=query[index]["uniProtKBCrossReferences"][refseqnt_index_list[0]]["properties"][0]["value"]
                    else:
                        rsntId = None
                        
            logger.info(f"uniprotId {uniprotId} -> RefSeqNTID {rsntId}")
            return rsntId

        uniprot_id = uniprot_id.split(":")[1]
        url = f"{self.base_url}search?query={uniprot_id}+AND+organism_id:9606&format=json&fields=accession,gene_names,organism_name,reviewed,xref_ensembl,xref_refseq,protein_name"
        for i in range(retries):
            try:
                response = requests.get(url, timeout=timeout)
                response.raise_for_status()
                results = response.json()["results"]
                if len(results) == 0:
                    return {}
                else:
                    #get values!
                    full_name = results[0]["proteinDescription"]["recommendedName"]["fullName"]["value"]
                    ensembl_id = _return_ensembl_from_id_and_uniprot_query(uniprot_id, results)
                    refseq_nt_id = _return_refseqnt_from_id_and_uniprot_query(uniprot_id, results)
                    return {"full_name" : full_name, "ensembl_id" : ensembl_id, "refseq_nt_id": refseq_nt_id}
            except requests.exceptions.RequestException:
                logger.warning(f"Failed to fetch UniProt data for {uniprot_id}")
            time.sleep(timeout)
        return {}

class EnsemblAPI:
    def __init__(self, retries=3, timeout=5):
        self.retries = retries
        self.timeout = timeout

    def get_sequence(self, ensembl_id, sequence_type="cdna"):
        """
        Given an Ensembl ID, returns the corresponding nucleotide sequence using the Ensembl API.
        """
        url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?object_type=transcript;type={sequence_type}"
        for i in range(self.retries):
            try:
                response = requests.get(url, headers={"Content-Type": "text/plain"}, timeout=self.timeout)
                response.raise_for_status()
                sequence = response.text
                logger.info(f"Received sequence for id {ensembl_id}.")
                return sequence
            except requests.exceptions.RequestException:
                logger.warning(f"Failed to fetch Ensembl sequence for {ensembl_id}. Retrying in {self.timeout} seconds.")
                time.sleep(self.timeout)
        logger.info(f"Failed to get sequence for id {ensembl_id} after {self.retries} retries.")
        return None
    
class HumanOrthologFinder:
    def __init__(self):
        self.zfin = ZFINHumanOrthologFinder()
        self.xenbase = XenbaseHumanOrthologFinder()
        self.mgi = MGIHumanOrthologFinder()
        self.rgd = RGDHumanOrthologFinder()

    def find_human_ortholog(self, product):
        """
        Finds the human ortholog for the given product.

        Args:
            product (str): The product (id) for which to find the human ortholog.

        Returns:
            The human gene symbol or None if no human ortholog was found.
        """
        if "ZFIN" in product:
            human_gene_symbol = self.zfin.find_human_ortholog(product)[0]
            return None if "Error" in human_gene_symbol else human_gene_symbol
        elif "Xenbase" in product:
            human_gene_symbol = self.xenbase.find_human_ortholog(product)[0]
            return None if "Error" in human_gene_symbol else human_gene_symbol
        elif "MGI" in product:
            human_gene_symbol = self.mgi.find_human_ortholog(product)
            return None if "Error" in human_gene_symbol else human_gene_symbol
        elif "RGD" in product:
            human_gene_symbol = self.rgd.find_human_ortholog(product)
            return None if "Error" in human_gene_symbol else human_gene_symbol
        else:
            logger.info(f"No database found for {product}")
            return None

class ZFINHumanOrthologFinder(HumanOrthologFinder):
    def __init__(self):
        self._filepath = "src_data_files/zfin_human_ortholog_mapping.txt"
        self._check_file()
        with open(self._filepath, "r") as read_content:
            self._readlines = read_content.readlines()
    
    def _check_file(self):
        os.makedirs(os.path.dirname(self._filepath), exist_ok=True)
        if not os.path.exists(self._filepath):
            url = "https://zfin.org/downloads/human_orthos.txt"
            response = requests.get(url)
            with open(self._filepath, "wb") as f:
                f.write(response.content)
            logger.info(f"Downloaded zfin_human_ortholog_mapping.txt to {self._filepath}")

    def find_human_ortholog(self, product_id):
        """
        If product_id is from the ZFIN database, searches through the zebrafish-human orthologs and returns the name of the
        symbol of the human gene ortholog.

        Returns:
        - [0]: gene symbol
        - [1]: long name of the gene
        """
        def _zfin_get_human_gene_symbol_from_line(line):
            """
            Splits zfin line and retrieves human gene symbol (full caps of zebrafish gene symbol)
            """
            # better, as zfin human ortholog sometimes has different name than the zebrafish gene
            human_symbol = line.split("\t")[3]
            human_gene_name = line.split("\t")[4]
            return human_symbol, human_gene_name # look at zfin orthologs txt file (in src_data_files) -> when you higlight a row, you see a TAB as '->' and a SPACEBAR as '.' -> splitting at \t means every [3] linesplit element is the human gene name

        product_id=product_id.split(":")[1] # eliminates 'ZFIN:' 
        for line in self._readlines:
            if product_id in line:
                e = _zfin_get_human_gene_symbol_from_line(line)
                human_symbol = e[0]
                human_gene_name = e[1]
                logger.info(f"[ Returning human symbol {human_symbol} and {human_gene_name}")
                return human_symbol, human_gene_name
        return [f"ZfinError_No-human-ortholog-found:product_id={product_id}"]

class XenbaseHumanOrthologFinder(HumanOrthologFinder):
    def __init__(self):
        self._filepath = "src_data_files/xenbase_human_ortholog_mapping.txt"
        self._check_file()
        with open(self._filepath, "r") as read_content:
            self._readlines = read_content.readlines()
    
    def _check_file(self):
        os.makedirs(os.path.dirname(self._filepath), exist_ok=True)
        if not os.path.exists(self._filepath):
            url = "https://download.xenbase.org/xenbase/GenePageReports/XenbaseGeneHumanOrthologMapping.txt"
            response = requests.get(url)
            with open(self._filepath, "wb") as f:
                f.write(response.content)
            logger.info(f"Downloaded xenbase_human_ortholog_mapping.txt to {self._filepath}")

    def find_human_ortholog(self, product_id):
        """
        Attempts to find a human ortholog from the xenbase database.
        Parameters:
        - product_id: eg. Xenbase:XB-GENE-495335 or XB-GENE-495335
        Returns: 
        - [0]: symbol of the human ortholog gene (eg. rsu1) or 'XenbaseError_no-human-ortholog-found'
        - [1]: long name of the gene
        """
        def _xenbase_get_human_symbol_from_line(line):
            """Splits xenbase line at tabs and gets human gene symbol (in full caps)"""
            symbol = str(line.split("\t")[2]).upper()
            name = str(line.split("\t")[3])
            return symbol, name

        product_id_short = ""
        if ":" in product_id: product_id_short = product_id.split(":")[1]
        else: product_id_short = product_id
        
        for line in self._readlines:
            if product_id_short in line:
                e = _xenbase_get_human_symbol_from_line(line)
                human_symbol = e[0]
                human_gene_name = e[1]
                logger.info(f"Found human ortholog {human_symbol}, name = {human_gene_name} for xenbase gene {product_id}")
                return human_symbol, human_gene_name
        return [f"[XenbaseError_No-human-ortholog-found:product_id={product_id}"]

class MGIHumanOrthologFinder(HumanOrthologFinder):
    def __init__(self):
        self._filepath = "src_data_files/mgi_human_ortholog_mapping.txt"
        self._check_file()
        with open(self._filepath, "r") as read_content:
            self._readlines = read_content.readlines()
    
    def _check_file(self):
        os.makedirs(os.path.dirname(self._filepath), exist_ok=True)
        if not os.path.exists(self._filepath):
            url = "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"
            response = requests.get(url)
            with open(self._filepath, "wb") as f:
                f.write(response.content)
            logger.info(f"Downloaded mgi_human_ortholog_mapping.txt to {self._filepath}")

    def find_human_ortholog(self, product_id):
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
                line = self._readlines[line_index+2]
                second_split = line.split("\t")
                if second_split[1] == "human":
                    logger.debug(f"Found keyword 'human' on secondpass line querying.")
                    return second_split[3]
                else:
                    # this still means no human ortholog!
                    # example: MGI:2660935 (Prl3d2) contains no "human" (neither i+1 nor i+2), also checked uniprot and no human gene for prl3d2 exists
                    return f"[MgiError_No-human-ortholog-found:product_id={product_id}"
            else: return split[3]

        logger.debug(f"Starting MGI search for {product_id}")
        product_id_short = ""
        if ":" in product_id:
            split = product_id.split(":")
            if len(split) == 3: product_id_short = split[2] # in case of MGI:xxx:xxxxx
            elif len(split) == 2: product_id_short = split[1] # in case of MGI:xxxxx
        else: product_id_short = product_id

        i = 0
        for line in self._readlines:
            if product_id_short in line:
                # if "mouse" gene smybol is found at line i, then human gene symbol will be found at line i+1
                logger.debug(f"i = {i}, product_id_short = {product_id_short}, line = {line}")
                human_symbol = _mgi_get_human_symbol_from_line(self._readlines[i+1], i)
                logger.info(f"Found human ortholog {human_symbol} for mgi gene {product_id}")
                return human_symbol # return here doesnt affect line counter 'i', since if gene is found i is no longer needed
            i += 1
        return f"[MgiError_No-human-ortholog-found:product_id={product_id}"

class RGDHumanOrthologFinder(HumanOrthologFinder):
    def __init__(self):
        self._filepath = "src_data_files/rgd_human_ortholog_mapping.txt"
        self._check_file()
        with open(self._filepath, "r") as read_content:
            self._readlines = read_content.readlines()
    

    def _check_file(self):
        os.makedirs(os.path.dirname(self._filepath), exist_ok=True)
        if not os.path.exists(self._filepath):
            url = "https://download.rgd.mcw.edu/pub/data_release/orthologs/RGD_ORTHOLOGS_Ortholog.txt"
            response = requests.get(url)
            with open(self._filepath, "wb") as f:
                f.write(response.content)
            logger.info(f"Downloaded rgd_human_ortholog_mapping.txt to {self._filepath}")

    def find_human_ortholog(self, product_id):
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

        product_id_short = ""
        if ":" in product_id: product_id_short = product_id.split(":")[1]
        else: product_id_short = product_id

        i = 0
        for line in self._readlines:
            if product_id_short in line:
                splitline_debug = line.split("\t")
                human_symbol = _rgd_get_human_symbol_from_line(line)
                logger.info(f"Found human ortholog {human_symbol} for RGD gene {product_id}")
                return human_symbol
        return f"[RgdError_No-human-ortholog-found:product_id={product_id}"

class GOTerm:
    def __init__(self, id: str, process: str, direction: str, name: Optional[str] = None, description: Optional[str] = None, weight: float = 1.0, products: List[str] = []):
        """
        A class representing a Gene Ontology term.

        Args:
            id (str): The ID of the GO term.
            process (str): The name of the biological process associated with the GO term.
            direction (str): The type of regulation associated with the GO term (e.g. positive or negative or general).
            name (str): Name (optional).
            description (str): A description of the GO term (optional).
            weight (float): The weight of the GO term.
            products (list): Products associated with the term (optional).
        """
        self.id = id
        self.process = process
        self.direction = direction
        self.name = name
        self.description = description
        self.weight = float(weight)
        self.products = products

    def fetch_name_description(self, api: GOApi):
        data = api.fetch_term_data(self.id)
        if data:
            self.name = data['label']
            self.description = data['definition']
            logger.info(f"Fetched name and description for GO term {self.id}")

    def fetch_products(self, api: GOApi):
        products = api.fetch_term_products(self.id)
        if products:
            self.products = products
            logger.info(f"Fetched products for GO term {self.id}")

    @classmethod
    def from_dict(cls, d: dict):
        """
        Creates a GOTerm object from a dictionary.

        Args:
            d (dict): A dictionary containing the GO term data.

        Returns:
            A new instance of the GOTerm class.
        """
        goterm = cls(d['id'], d['process'], d['direction'], d.get('name'), d.get('description'), d.get('weight', 1.0), d.get('products', []))
        return goterm

class Product:
    def __init__(self, id_synonyms: List[str], uniprot_id: str = None, description: str = None, ensembl_id: str= None, refseq_nt_id: str = None, mRNA: str = None, scores: dict = None):
        """
        A class representing a product (e.g. a gene or protein).

        Args:
            id_synonyms (str): The list of ID of the product and synonyms. -> after ortholog translation it turns out that some products are the same. 
            uniprot_id (str): The UniProt ID of the product.
            description (str): A description of the product.
            ensembl_id (str):
            mRNA (str): The mRNA sequence of the product.
            scores (dict): A dictionary of scores associated with the product (e.g. expression score, functional score).
        """
        self.id_synonyms = id_synonyms
        self.uniprot_id = uniprot_id
        self.description = description
        self.ensembl_id = ensembl_id
        self.refseq_nt_id = refseq_nt_id
        self.mRNA = mRNA
        self.scores = {} if scores is None else scores.copy()
    
    def fetch_UniprotID(self, human_ortolog_finder: HumanOrthologFinder, uniprot_api: UniProtAPI) -> None:
        if len(self.id_synonyms) == 1 and 'UniProtKB' in self.id_synonyms[0]:
            self.uniprot_id = self.id_synonyms[0]
        elif len(self.id_synonyms) == 1:
            human_ortholog_gene_id = human_ortolog_finder.find_human_ortholog(self.id_synonyms[0])
            if human_ortholog_gene_id is not None:
                uniprot_id = uniprot_api.get_uniprot_id(human_ortholog_gene_id)
                if uniprot_id is not None:
                    self.uniprot_id = uniprot_id
        # TODO: What if len(self.id_synonyms) > 1 ?

    def fetch_Uniprot_info(self, uniprot_api: UniProtAPI) -> None:
        """
        includes description, ensembl_id and refseq_nt_id
        """
        if self.uniprot_id == None: return
        info_dict = uniprot_api.get_uniprot_info(self.uniprot_id)
        full_name = info_dict.get("full_name")
        ensembl_id = info_dict.get("ensembl_id")
        refseq_nt_id = info_dict.get("refseq_nt_id")
        if full_name is not None:
            self.description = full_name
        if ensembl_id is not None:
            self.ensembl_id = ensembl_id
        if refseq_nt_id is not None:
            self.refseq_nt_id = refseq_nt_id
    
    def fetch_mRNA_sequence(self, ensembl_api: EnsemblAPI) -> None:
        if self.ensembl_id is None: return
        sequence = ensembl_api.get_sequence(self.ensembl_id)    
        if sequence is not None:
            self.mRNA = sequence    

    @classmethod
    def from_dict(cls, d: dict) -> 'Product':
        """
        Class method to create a new Product instance from a dictionary.

        Args:
            d (dict): The dictionary containing the data to create the Product instance.

        Returns:
            Product: A new Product instance created from the input dictionary.
        """
        return cls(d['id_synonyms'], d.get('uniprot_id'), d.get('description'), d.get('ensembl_id'), d.get('refseq_nt_id'), d.get('mRNA'), d.get('scores'))

class miRNA:
    def __init__(self, id: str, sequence: str = None, mRNA_overlaps: Dict[str, float] = None, scores: Dict[str, float] = None) -> None:
        """
        Initializes an instance of miRNA class.

        Args:
        - id: a string that uniquely identifies this miRNA.
        - sequence: an optional string that represents the sequence of this miRNA.
        - mRNA_overlaps: an optional dictionary that represents the overlaps of this miRNA with mRNA sequences.
        - scores: an optional dictionary that represents the scores of this miRNA.

        Returns: None
        """
        self.id = id
        self.sequence = sequence
        self.mRNA_overlaps = {} if mRNA_overlaps is None else mRNA_overlaps.copy()
        self.scores = {} if scores is None else scores.copy()

    @classmethod
    def from_dict(cls, d: dict) -> 'miRNA':
        """
        Creates a new instance of miRNA class based on the values in the input dictionary.

        Args:
        - d: a dictionary that represents the values of the miRNA.

        Returns: 
        - A new instance of miRNA class based on the values in the input dictionary.
        """
        return cls(d['id'], d.get('sequence'), d.get('mRNA_overlaps'), d.get('scores'))

class ReverseLookup:
    def __init__(self, goterms: List[GOTerm], target_processes: List[Dict[str, str]], products: List[Product] = [], miRNAs: List[miRNA] = [], miRNA_overlap_treshold: float = 0.6, model_name: str = "model", execution_times: dict = {}):
        """
        A class representing a reverse lookup for gene products and their associated Gene Ontology terms.

        Args:
            goterms (set): A set of GOTerm objects.
            target_processes (list): A list of dictionaries containing process names and directions.
            products (set, optional): A set of Product objects. Defaults to an empty set.
        """
        self.goterms = goterms
        self.products = products
        self.target_processes = target_processes
        self.miRNAs = miRNAs
        self.miRNA_overlap_treshold = miRNA_overlap_treshold
        self.model_name = model_name # arbitrary model name
        self.execution_times = execution_times # dict of execution times, logs of runtime for functions

        self.timer = Timer()
        
        self.save_path_root = os.path.join("program_data_files", "models", model_name) # path where model files are saved
        if not os.path.exists(self.save_path_root):
            logger.debug(f"Created directory {self.save_path_root} for model {self.model_name}.")
            os.makedirs(self.save_path_root)
        else:
            logger.debug(f"Directory {self.save_path_root} already exists. Continuing.")

    def fetch_all_go_term_names_descriptions(self):
        """
        Iterates over all GOTerm objects in the go_term set and calls the fetch_name_description method for each object.
        """
        self.timer.set_start_time()

        api = GOApi()
        with logging_redirect_tqdm():
            for goterm in tqdm(self.goterms):
                if goterm.name == "":
                    goterm.fetch_name_description(api)

        if "fetch_all_go_term_names_descriptions" not in self.execution_times: # to prevent overwriting on additional runs of the same model name
            self.execution_times["fetch_all_go_term_names_descriptions"] = self.timer.get_elapsed_time()
        self.timer.print_elapsed_time()

    def fetch_all_go_term_products(self):
        """
        Iterates over all GOTerm objects in the go_term set and calls the fetch_products method for each object.
        """
        self.timer.set_start_time()

        api = GOApi()
        with logging_redirect_tqdm():
            for goterm in tqdm(self.goterms):
                goterm.fetch_products(api)
        
        if "fetch_all_go_term_products" not in self.execution_times:
            self.execution_times["fetch_all_go_term_products"] = self.timer.get_elapsed_time()
        self.timer.print_elapsed_time()

    def create_products_from_goterms(self) -> None:
        """
        This method creates Product objects from the set of products contained in each GOTerm object and
        adds them to the ReverseLookup object's products list.

        The method iterates over each GOTerm object in the goterms set and retrieves the set of products associated
        with that GOTerm object. It then adds these products to a products_set, which is a set object that ensures
        that no duplicate products are added.

        Finally, the method iterates over each product in the products_set and creates a new Product object from the
        product ID using the Product.from_dict() classmethod. The resulting Product objects are added to the
        ReverseLookup object's products list.

        Args:
            None

        Returns:
            None
        """
        self.timer.set_start_time()

        # Create an empty set to store unique products
        products_set = set()
        
        # Iterate over each GOTerm object in the go_term set and retrieve the set of products associated with that GOTerm
        # object. Add these products to the products_set.
        for term in self.goterms:
            products_set.update(term.products)
        
        # Iterate over each product in the products_set and create a new Product object from the product ID using the
        # Product.from_dict() classmethod. Add the resulting Product objects to the ReverseLookup object's products list.
        for product in products_set:
            self.products.append(Product.from_dict({'id_synonyms':[product]}))
        logger.info(f"Created Product objects from GOTerm object definitions")

        if "create_products_from_goterms" not in self.execution_times:
            self.execution_times["create_products_from_goterms"] = self.timer.get_elapsed_time()
        self.timer.print_elapsed_time()

    def fetch_UniprotID_products(self) -> None:
        self.timer.set_start_time()

        try:   
            human_ortolog_finder = HumanOrthologFinder()
            uniprot_api = UniProtAPI()
            # Iterate over each Product object in the ReverseLookup object.
            with logging_redirect_tqdm():
                for product in tqdm(self.products):
                    # Check if the Product object doesn't have a UniProt ID.
                    if product.uniprot_id == None:
                        # If it doesn't, fetch UniProt data for the Product object.
                        product.fetch_UniprotID(human_ortolog_finder, uniprot_api)
        except Exception as e:
            # If there was an exception while fetching UniProt data, save all the Product objects to a JSON file.
            self.save_products_to_datafile('crash_products.json')
            # Re-raise the exception so that the caller of the method can handle it.
            raise e
        
        if "fetch_uniprotid_products" not in self.execution_times:
            self.execution_times["fetch_uniprotid_products"] = self.timer.get_elapsed_time()
        self.timer.print_elapsed_time()

    def prune_products(self) -> None:
        self.timer.set_start_time()

        # Create a dictionary that maps UniProt ID to a list of products
        reverse_uniprotid_products = {}
        for product in self.products:
            if product.uniprot_id is not None:
                reverse_uniprotid_products.setdefault(product.uniprot_id, []).append(product)

        # For each UniProt ID that has more than one product associated with it, create a new product with all the synonyms
        # and remove the individual products from the list
        for uniprot_id, product_list in reverse_uniprotid_products.items():
            if len(product_list) > 1:
                id_synonyms = []
                for product in product_list:
                    self.products.remove(product)
                    id_synonyms.extend(product.id_synonyms)
                # Create a new product with the collected information and add it to the product list
                self.products.append(Product(id_synonyms, uniprot_id, product_list[0].description, product_list[0].ensembl_id, product_list[0].refseq_nt_id, product_list[0].mRNA, {}))

        if "prune_products" not in self.execution_times:
            self.execution_times["prune_products"] = self.timer.get_elapsed_time()
        self.timer.print_elapsed_time()

    def fetch_Uniprot_infos(self) -> None:
        self.timer.set_start_time()

        try:   
            uniprot_api = UniProtAPI()
            # Iterate over each Product object in the ReverseLookup object.
            with logging_redirect_tqdm():
                for product in tqdm(self.products):
                    # Check if the Product object doesn't have a UniProt ID.
                    if (product.description is None or product.ensembl_id is None or product.refseq_nt_id is None) and product.uniprot_id is not None:
                        # If it doesn't, fetch UniProt data for the Product object.
                        product.fetch_Uniprot_info(uniprot_api)
        except Exception as e:
            raise e

        if "fetch_uniprot_infos" not in self.execution_times:
            self.execution_times["fetch_uniprot_infos"] = self.timer.get_elapsed_time()
        self.timer.print_elapsed_time()

    def score_products(self) -> None:
        self.timer.set_start_time()

        # create a Scoring object with the current instance of miRGadget as the argument
        scoring = Scoring(self)
        # redirect the tqdm logging output to the logging module to avoid interfering with the normal output
        with logging_redirect_tqdm():
            # iterate over each Product object in self.products and score them using the Scoring object
            for product in tqdm(self.products):
                # calculate the advanced product score for the current Product object using the Scoring object
                adv_score_result = scoring.adv_product_score(product)
                # calculate the number of GO terms per process for the current Product object using the Scoring object
                dict_of_terms_per_process = scoring.nterms(product)
                # add the advanced product score and GO terms per process to the Product object's scores attribute
                product.scores = {"adv_score": adv_score_result, "nterms": dict_of_terms_per_process}
        
        if "score_products" not in self.execution_times:
            self.execution_times["score_products"] = self.timer.get_elapsed_time()
        self.timer.print_elapsed_time()

    def fetch_mRNA_sequences(self) -> None:
        self.timer.set_start_time()
        
        try:   
            ensembl_api = EnsemblAPI()
            # Iterate over each Product object in the ReverseLookup object.
            with logging_redirect_tqdm():
                for product in tqdm(self.products):
                    # Check if the Product object doesn't have a EnsemblID
                    if product.mRNA == None and product.ensembl_id is not None:
                        # If it has, fetch mRNA sequence data for the Product object.
                        product.fetch_mRNA_sequence(ensembl_api)
        except Exception as e:
            raise e
        
        if "fetch_mRNA_sequences" not in self.execution_times:
            self.execution_times["fetch_mRNA_sequences"] = self.timer.get_elapsed_time()
        self.timer.print_elapsed_time()

    def predict_miRNAs(self, prediction_type: str = 'miRDB') -> None:
        self.timer.set_start_time()

        # check the prediction type
        if prediction_type == 'miRDB':
            # use the miRDB60predictor to predict miRNAs
            predictor = miRDB60predictor()
            # iterate through each product and predict miRNAs
            with logging_redirect_tqdm():
                for product in tqdm(self.products):
                    match_dict = predictor.predict_from_product(product)
                    # if there are matches, add them to the corresponding miRNA objects
                    if match_dict is not None:
                        for miRNA_id, match in match_dict.items():
                            # check if the miRNA already exists in the list of miRNAs
                            for mirna in self.miRNAs:
                                if mirna.id == miRNA_id:
                                    mirna.mRNA_overlaps[product.uniprot_id] = match
                                    break
                            # if the miRNA doesn't exist in the list, create a new miRNA object
                            else:
                                self.miRNAs.append(miRNA(miRNA_id, mRNA_overlaps={product.uniprot_id : match}))
        elif prediction_type == 'other_type':
            # do something else
            pass
        else:
            # raise an error if the prediction type is invalid
            raise ValueError("Invalid prediction type")
        
        if "predict_miRNAs" not in self.execution_times:
            self.execution_times["predict_miRNAs"] = self.timer.get_elapsed_time()
        self.timer.print_elapsed_time()

    def change_miRNA_overlap_treshold(self, treshold: float, yes: bool = False) -> None:
        self.miRNA_overlap_treshold = treshold
        logger.warning(f"Sorry, but changing the treshold will delete all the calculated miRNA scores. You will have to calculate them again!")
        if not yes:
            confirmation = input(f"Are you sure you want to proceed? (y/n) ")
            if confirmation.lower() != 'y':
                print("Aborting operation.")
                return
        for _miRNA in self.miRNAs:
            _miRNA.scores = {}

    def score_miRNAs(self) -> None:
        self.timer.set_start_time()

        # create a Scoring object to compute scores
        scoring = Scoring(self)
        
        with logging_redirect_tqdm():
            # iterate over miRNAs using tqdm for progress tracking
            for mirna in tqdm(self.miRNAs):
                # if there is no overlap, skip the miRNA
                if not mirna.mRNA_overlaps:
                    continue
                
                # compute the basic score for the miRNA
                basic_score_result = scoring.basic_mirna_score(mirna, self.miRNA_overlap_treshold)
                
                # store the score result in the miRNA object
                mirna.scores = {"basic_score": basic_score_result}
        
        if "score_miRNAs" not in self.execution_times:
            self.execution_times["score_miRNAs"] = self.timer.get_elapsed_time()
        self.timer.print_elapsed_time()
    
    def compute_all(self, report_filepath: str = "") -> None:
        # Fetch all GO term names and descriptions
        self.fetch_all_go_term_names_descriptions() # c1-TN
        self.save_model()
        # Fetch all GO term products
        self.fetch_all_go_term_products() # c2-TP
        self.save_model()
        # Create products from GO terms
        self.create_products_from_goterms() # c3-P
        self.save_model()
        # Fetch UniProt ID products
        self.fetch_UniprotID_products() # c4-UidP
        self.save_model()
        # Prune products
        self.prune_products() # c5-Pr
        self.save_model()
        # Fetch UniProt information
        self.fetch_Uniprot_infos() # c6-Uinf
        self.save_model()
        # Score products
        self.score_products() # c7-SC
        self.save_model()
         # Optional: Fetch mRNA sequences
        self.fetch_mRNA_sequences() # c8-mRNA
        self.save_model()
        # Predict miRNAs
        self.predict_miRNAs() # c9-miRNA
        self.save_model()
        # Score miRNAs
        self.score_miRNAs() # c10-miRNASC
        self.save_model()
        # Generate report
        if report_filepath == "": report_filepath = os.path.join(self.save_path_root, "report.txt")
        report = ReverseLookup.ReportGenerator(self, verbosity=3)
        report.general_report(report_filepath)

#housekeeping functions

    def get_all_goterms_for_product(self, product: Product | str) -> List[GOTerm]:
        if isinstance(product, str):
            for prod in self.products:
                if prod.uniprot_id == product:
                    product = prod
                    break
        
        goterms_list = []
        for goterm in self.goterms:
            if any(product_id in goterm.products for product_id in product.id_synonyms):
                goterms_list.append(goterm)
        return goterms_list

    def list_goterms_id(self) -> List[str]:
        """
        Returns a list of all GO term IDs in the GO ontology.
        """
        # Use a list comprehension to extract IDs from the GO terms and return the resulting list
        return [goterm.id for goterm in self.goterms]
    
    @classmethod
    def analyze_saves(self, model_name: str) -> tuple[int,List]:
        """
        Analyzes the amount of unique model.save jsons inside program_data_files/models/model_name/

        Returns:
          - [0]: (int) count of model save files
          - [1]: [List] relative paths to model save files
        
        Returns [-1, []] if program_data_files/models/model_name/ doesn't exist
        """
        dirpath = os.path.join("program_data_files","models",model_name)
        count = 0
        saved_model_files = []
        if os.path.exists(dirpath):
            files = os.listdir(dirpath)
            for filename in files:
                if "data" in filename:
                    count += 1
                    saved_model_files.append(os.path.join("program_data_files", "models", model_name, filename))
            return count, saved_model_files
        else:
            return -1,[]
        
    
    def save_model(self, filepath: str = "") -> None:
        """
        Saves the model as a JSON. If input filepath is not provided, model is saved to self.save_path_root, else the model
        is saved to the input filepath.

        The output JSON is saved as data.json
        """
        data = {}
        #save options - currently not used

        #save target_process
        data['model_name'] = self.model_name
        data['execution_times'] = self.execution_times
        data['target_processes'] = self.target_processes
        data['miRNA_overlap_treshold'] = self.miRNA_overlap_treshold
        #save goterms
        for goterm in self.goterms:
            data.setdefault('goterms', []).append(goterm.__dict__)
        #save products
        for product in self.products:
            data.setdefault('products', []).append(product.__dict__)
        #save miRNAs
        for miRNA in self.miRNAs:
            data.setdefault('miRNAs', []).append(miRNA.__dict__)
        #write to file
        if filepath == "": filepath = os.path.join(self.save_path_root, "data.json") # if filepath not provided, use self.save_path_root to compute destination json location
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=4)

    @classmethod
    def load_model(cls, filepath: str = "", model_name: str = "") -> 'ReverseLookup':
        """
        Loads the model from filepath. If filepath is specified (eg. "diabetes_angio_1/data.json"), the model_name is ignored.
        If filepath is left empty ("") and model_name is specified, then program_data_files/models/model_name/data.json is loaded.

        If there is more than one model save file located at program_data_files/models/model_name/, then the user is given an option
        to choose among the model save files. This function displays only the filenames with the keyword "data" as valid load files,
        therefore make sure, to save all of the model json files with the "data" keyword in the filename.
        """
        save_files_analysis = cls.analyze_saves(model_name)
        if filepath == "" and model_name != "" and save_files_analysis[0] == 1: # if there is exactly one saved version in save_files_analysis
            filepath = os.path.join("program_data_files", "models", model_name, "data.json")
        elif filepath == "" and model_name != "" and save_files_analysis[0] > 1: # if there are more than one saved versions
            logger.info(f"Model {model_name} has multiple save files. Type the number of the save file you wish to choose.")
            i=0
            for filepath in save_files_analysis[1]:
                filename = filepath.split(os.sep)[-1]
                logger.info(f"[{i}]: {filename}")
                i+=1
            user_input = -1
            while True:
                user_input = input("Type the number of the save file you wish to choose: ") 
                try:
                    user_input = int(user_input)
                except ValueError:
                    logger.info("Invalid input! Please enter an integer.")
                if isinstance(user_input, int) and user_input < len(save_files_analysis[1]):
                    break
                else:
                    logger.info(f"Invalid input. Please enter an integer no greater than {len(save_files_analysis[1])}.")
            filepath = save_files_analysis[1][user_input]
            logger.info(f"You chose: {filepath} as the model load file")

        with open(filepath, "r") as f:
            data = json.load(f)
        
        model_name = data['model_name']
        execution_times = data['execution_times']
        target_processes = data['target_processes']
        miRNA_overlap_treshold = data['miRNA_overlap_treshold']

        goterms = []
        for goterm_dict in data['goterms']:
            goterms.append(GOTerm.from_dict(goterm_dict))

        products = []
        for product_dict in data.get('products', []):
            products.append(Product.from_dict(product_dict))

        miRNAs = []
        for miRNAs_dict in data.get('miRNAs', []):
            miRNAs.append(miRNA.from_dict(miRNAs_dict))
        
        if model_name:
            return cls(goterms, target_processes, products, miRNAs, miRNA_overlap_treshold, model_name=model_name, execution_times=execution_times)
        else:
            return cls(goterms, target_processes, products, miRNAs, miRNA_overlap_treshold, execution_times=execution_times)

    @classmethod
    def from_input_file(cls, filepath: str, mod_name: str = "model") -> 'ReverseLookup':
        """
        Creates a ReverseLookup object from a JSON file.

        Args:
            filepath (str): The path to the input file.
            mode_name (str): An arbitrary name of this model, or leave default ('model').

        Returns:
            ReverseLookup: A ReverseLookup object.
        """
        # Define constants used in parsing the file
        LINE_ELEMENT_DELIMITER = '\t' # Data is tab separated
        COMMENT_DELIMITER = "#" # Character used to denote a comment
        LOGIC_LINE_DELIMITER = "###" # Special set of characters to denote a "logic line"

        def process_comment(line):
            """
            Processes a comment in the line: returns the part of the line before the comment.

            Parameters:
            - line: the line whose comment to process
            """
            if LOGIC_LINE_DELIMITER in line:
                # Logic lines should be marked with "###" at the start. For a logic line, the returned result is line without the line_keep_delimiter
                return line.replace(LOGIC_LINE_DELIMITER, "")
                
            if COMMENT_DELIMITER in line:
                return line.split(COMMENT_DELIMITER)[0]
            else:
                return line
        
        # process the input filepath
        with open(filepath, "r") as read_content:
            target_processes = []
            go_terms = []
            read_lines = read_content.read().splitlines()[2:] #skip first 2 lines
            section = "" #what is the current section i am reading
            for line in read_lines:
                line = process_comment(line)
                if line == "":
                    continue
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
                    chunks = line.split(LINE_ELEMENT_DELIMITER)
                elif section == "process":
                    chunks = line.split(LINE_ELEMENT_DELIMITER)
                    target_processes.append({"process":chunks[0], "direction":chunks[1]})
                elif section == "GO":
                    chunks = line.split(LINE_ELEMENT_DELIMITER)
                    if len(chunks) == 5:
                        d = {"id":chunks[0], "process":chunks[1], "direction":chunks[2], "weight":chunks[3], "description": chunks[4]}
                    else:
                        d = {"id":chunks[0], "process":chunks[1], "direction":chunks[2], "weight":chunks[3]}
                    go_terms.append(GOTerm.from_dict(d))
        
        # process program_data_files/models/self.model_name/data.json if it exists
        model_existing_data_filepath = os.path.join("program_data_files", "models", mod_name, "data.json")
        if os.path.exists(model_existing_data_filepath):
            with open(model_existing_data_filepath, "r") as f:
                model_existing_data_json = json.load(f)
        
        if model_existing_data_json:
            execution_times = model_existing_data_json['execution_times']
            target_processes = model_existing_data_json['target_processes']
            miRNA_overlap_treshold = model_existing_data_json['miRNA_overlap_treshold']

            goterms = []
            for goterm_dict in model_existing_data_json['goterms']:
                goterms.append(GOTerm.from_dict(goterm_dict))

            products = []
            for product_dict in model_existing_data_json.get('products', []):
                products.append(Product.from_dict(product_dict))

            miRNAs = []
            for miRNAs_dict in model_existing_data_json.get('miRNAs', []):
                miRNAs.append(miRNA.from_dict(miRNAs_dict))
            
            # data.json for this model exists, use it to load the model
            return cls(goterms, target_processes, products, miRNAs, miRNA_overlap_treshold, model_name=mod_name, execution_times=execution_times)
        
        else:
            # data.json for this model doesn't exist yet
            return cls(go_terms, target_processes, model_name=mod_name)
        
    
    @classmethod
    def from_dict(cls, data: Dict[str, List[Dict]]) -> 'ReverseLookup':
        """
        Creates a ReverseLookup object from a dictionary.

        Args:
            data (dict): A dictionary containing lists of GOTerm and target_processes.

        Returns:
            ReverseLookup: A ReverseLookup object.
        """

        goterms = [GOTerm.from_dict(d) for d in data['goterms']]
        target_processes = data['target_processes']
        return cls(goterms, target_processes)

class Timer:
    def __init__(self):
        self.start_time = time.time()
    
    def set_start_time(self):
        """
        Sets a new reference start time.
        """
        self.start_time = time.time()
    
    def get_elapsed_seconds(self) -> int:
        """
        Returns the amount of seconds unformatted (contains decimal places)
        """
        return time.time() - self.start_time
    
    def get_elapsed_time(self) -> str:
        """
        Gets elapsed time in hh mm ss format.
        """
        sec = int(self.get_elapsed_seconds())
        td = timedelta(seconds=sec)
        return str(td)
    
    def print_elapsed_time(self, useLogger: bool = True, prefix: str = "Elapsed: "):
        """
        Prints the elapsed time in hh mm ss format. 
        
        Args:
          - useLogger: if True, then logger.info is used. If false, then print is used.
          - prefix: the string you want to use as a prefix
        """
        if useLogger:
            logger.info(f"{prefix}{self.get_elapsed_time()}")
        else:
            print(f"{prefix}{self.get_elapsed_time()}")

    
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

class Scoring:
    def __init__(self, reverse_lookup: ReverseLookup) -> None:
        self.reverse_lookup = reverse_lookup
    
    def _opposite_direction(self, direction: str) -> str:
        if direction == "0":
            return "0"
        elif direction == "+":
            return "-"
        elif direction == "-":
            return "+"

    def adv_product_score(self, product: Product) -> float:
        a = 10
        b1 = 2
        b2 = 0.5
        c1 = 1
        c2 = 0.1

        goterms_list = self.reverse_lookup.get_all_goterms_for_product(product)
        score = 0.0
        # Check if all target processes are regulated in the same direction as the GOTerms in the list
        # and none of them are regulated in the opposite direction
        if (
            # Check if all processes in target_processes have a GOTerm in goterms_list that regulates it in the same direction
            all(
                any(process['direction'] == goterm.direction and process['process'] == goterm.process for goterm in goterms_list)
                for process in self.reverse_lookup.target_processes
            )
            # Check if none of the processes in target_processes have a GOTerm in goterms_list that regulates it in the opposite direction
            and not any(
                any(self._opposite_direction(process['direction']) == goterm.direction and process['process'] == goterm.process for goterm in goterms_list)
                for process in self.reverse_lookup.target_processes
            )
        ):
            # If all target processes are regulated in the same direction, add a points to the score
            score += a

        # Check if all target processes are regulated in the opposite direction as the GOTerms in the list
        # and none of them are regulated in the same direction
        if (
            # Check if all processes in target_processes have a GOTerm in goterms_list that regulates it in the opposite direction
            all(
                any(self._opposite_direction(process['direction']) == goterm.direction and process['process'] == goterm.process for goterm in goterms_list)
                for process in self.reverse_lookup.target_processes
            )
            # Check if none of the processes in target_processes have a GOTerm in goterms_list that regulates it in the same direction
            and not any(
                any(process['direction'] == goterm.direction and process['process'] == goterm.process for goterm in goterms_list)
                for process in self.reverse_lookup.target_processes
            )
        ):
            # If all target processes are regulated in the opposite direction, subtract a points from the score
            score -= a

        # Calculate the score based on the number of processes in target_processes that are regulated
        # by GOTerms in the same direction as defined in the list
        score += sum(
            (b1**(b2 * sum(
                # Check if the direction and process in the process dict matches with the direction and process in any GOTerm dict
                goterm.weight for goterm in goterms_list if process['direction'] == goterm.direction and process['process'] == goterm.process)))
                for process in self.reverse_lookup.target_processes
        )
        # Calculate the score based on the number of processes in target_processes that are regulated
        # by GOTerms in the oposite direction as defined in the list
        score -= sum(
            (b1**(b2 * sum(
                # Check if the direction and process in the process dict matches with the direction and process in any GOTerm dict
                goterm.weight for goterm in goterms_list if self._opposite_direction(process['direction']) == goterm.direction and process['process'] == goterm.process)))
                for process in self.reverse_lookup.target_processes
        )

        # Calculate the score by multiplying the current score with a factor based on the number of GOTerms with direction "0"
        score = score * (
            c1  # Start with a base factor of 1
            + (c2  # Add a factor based on a constant value c
                * sum(  # Multiply c by the sum of weights of all GOTerms with direction "0"
                    goterm.weight  # Get the weight of each GOTerm
                    for goterm in goterms_list  # Iterate over all GOTerms in the list
                    if goterm.direction == "0"  # Only consider GOTerms with direction "0"
                )
            )
        )

        return score

    def nterms(self, product: Product) -> dict:

        goterms_list = self.reverse_lookup.get_all_goterms_for_product(product)

        # Create an empty dictionary to store the count of GOTerms for each process and direction
        nterms_dict = {}
        
        # Iterate over each process in the target_processes list
        for process in self.reverse_lookup.target_processes:
            # Count the number of GOTerms that have a direction of "+" and a process matching the current process
            nterms_dict[f"{process['process']}+"] = sum(1 for goterm in goterms_list if (goterm.direction == "+" and process['process'] == goterm.process))
            
            # Count the number of GOTerms that have a direction of "-" and a process matching the current process
            nterms_dict[f"{process['process']}-"] = sum(1 for goterm in goterms_list if (goterm.direction == "-" and process['process'] == goterm.process))
            
            # Count the number of GOTerms that have a direction of "0" and a process matching the current process
            nterms_dict[f"{process['process']}0"] = sum(1 for goterm in goterms_list if (goterm.direction == "0" and process['process'] == goterm.process))
        
        # Return the dictionary containing the count of GOTerms for each process and direction
        return nterms_dict
    
    def inhibited_products_id(self, mirna: miRNA, treshold: float = 0.6) -> List[str]:
        inhibited_product_id = []
        for product_id, overlap in mirna.mRNA_overlaps.items():
            if overlap >= treshold:
                inhibited_product_id.append(product_id)
        return inhibited_product_id

    def basic_mirna_score(self, mirna: miRNA, treshold: float = 0.6) -> float:
        score = 0.0
        for product_id, overlap in mirna.mRNA_overlaps.items():
            product = next((x for x in self.reverse_lookup.products if x.uniprot_id == product_id), None)
            if product is not None:
                if overlap >= treshold: #inhibited
                    a = -1 #deduct the score, since high score indicates the products is favourable for our target processes
                else:
                    a = 1
                score += a * product.scores["adv_score"]
        return score

class ReportGenerator:
    def __init__(self, reverse_lookup: ReverseLookup, verbosity: int = 1, top_n: int = 5, width: int = 80):
        # initialize the report generator with a reverse lookup object and parameters for verbosity, top_n, and width
        self.reverse_lookup = reverse_lookup  # the reverse lookup object to use for generating the report
        self.width = width  # the width of the report output
        self.top_n = top_n  # the number of top results to display in the report
        self.verbosity = verbosity  # the level of detail to include in the report

    
    def _generate_header(self) -> str:
        header = "Gene Ontology Reverse Lookup Tool".center(self.width)+"\n"
        header += "Authors: Vladimir Smrkolj (SI), Aljosa Skorjanc (SI)".center(self.width)+"\n"
        header += "March 2023".center(self.width) + "\n"
        return header
    
    def _generate_section(self, text: str) -> str:
        string = "-"*self.width+"\n"
        string += text.center(self.width)+"\n"
        string += "-"*self.width+"\n"
        return string

    def _generate_goterm_per_process_table(self) -> str:
        # Use a dictionary comprehension for the grouped_goterms initialization
        grouped_goterms = {(target["process"], direction): [] for target in self.reverse_lookup.target_processes for direction in ("+", "-", "0")}
        
        # Use a for-loop to populate the grouped_goterms dictionary
        for goterm in self.reverse_lookup.goterms:
            key = (goterm.process, goterm.direction)
            grouped_goterms[key].append(goterm)

        string = "GO TERMS PER PROCESS" + "\n"

        if self.verbosity >= 1:
            table = [["Process", "+", "-", "0", "Total"]]
            table.extend([
                [target["process"],
                len(grouped_goterms[(target["process"], "+")]),
                len(grouped_goterms[(target["process"], "-")]),
                len(grouped_goterms[(target["process"], "0")]),
                sum([len(grouped_goterms[(target["process"], "+")]), len(grouped_goterms[(target["process"], "-")]), len(grouped_goterms[(target["process"], "0")])])
                ]
                for target in self.reverse_lookup.target_processes
            ])
            table.append(["Total", sum([a[1] for a in table[1:]]), sum([a[2] for a in table[1:]]), sum([a[3] for a in table[1:]]), sum([a[4] for a in table[1:]])])
            string += tabulate(table, headers="firstrow", tablefmt="grid").center(self.width) + "\n"

        if self.verbosity == 2:
            table = [["Process", "+", "-", "0"]]
            table.extend([
                [target["process"],
                '\n'.join(str(g.id) for g in grouped_goterms[(target["process"], "+")]),
                '\n'.join(str(g.id) for g in grouped_goterms[(target["process"], "-")]),
                '\n'.join(str(g.id) for g in grouped_goterms[(target["process"], "0")]) ]
                for target in self.reverse_lookup.target_processes
            ])
            string += tabulate(table, headers="firstrow", tablefmt="grid").center(self.width) + "\n"

        if self.verbosity == 3:
            table = [["Process", "+", "-", "0"]]
            table.extend([
                [target["process"],
                '\n'.join(str(f"{g.id} - {g.name}") for g in grouped_goterms[(target["process"], "+")]),
                '\n'.join(str(f"{g.id} - {g.name}") for g in grouped_goterms[(target["process"], "-")]),
                '\n'.join(str(f"{g.id} - {g.name}") for g in grouped_goterms[(target["process"], "0")]) ]
                for target in self.reverse_lookup.target_processes
            ])
            string += tabulate(table, headers="firstrow", tablefmt="grid") + "\n"

        return string

    def _generate_goterms_statistics(self) -> str:
        string = "GO TERMS STATISTICS\n"
        
        # Calculate statistics and format the string
        if self.verbosity >= 1:
            products_per_goterm = [len(g.products) for g in self.reverse_lookup.goterms]
            min_g, max_g, avg_g = min(products_per_goterm), max(products_per_goterm), sum(products_per_goterm) / len(products_per_goterm)
            string += f"Products per GO Term (min - avg - max): {min_g} - {avg_g:.0f} - {max_g}\n"

        return string
    
    def _generate_top_bottom_products_summary(self) -> str:
        # Initialize the summary string with the header
        string = f"TOP and BOTTOM {self.top_n} PRODUCTS\n"

        # Get the top and bottom products based on the advanced score
        sorted_products = sorted(self.reverse_lookup.products, key=lambda x: x.scores["adv_score"], reverse=True)
        top_products = sorted_products[:self.top_n]
        bottom_products = sorted_products[-self.top_n:]

        # If verbosity is at least 1, create a table of the top and bottom products with their scores and descriptions
        if self.verbosity >= 1:
            # Create the table as a list of lists and append each row
            table = [["UniProtID", "Score", "Protein name"]]
            for product in top_products:
                table.append([product.uniprot_id, f"{product.scores['adv_score']:.2f}", product.description])
            # Add a separator row and append each row for the bottom products
            table.append(["----", "----", "----"])
            for product in bottom_products:
                table.append([product.uniprot_id, f"{product.scores['adv_score']:.2f}", product.description])
            # Add the table to the summary string
            string += tabulate(table, headers="firstrow", tablefmt="grid") + "\n\n"

        # If verbosity is at least 2, add details about the GO terms for each top and bottom product
        if self.verbosity >= 2:
            # Define a function to create the table of GO terms for a product
            def create_go_table(product):
                table_go = [["GO term", "GO label", "GO description"]]
                for goterm in self.reverse_lookup.get_all_goterms_for_product(product):
                    table_go.append([goterm.id, goterm.name, goterm.description])
                return table_go

            # Add the details for the top products
            for product in top_products:
                string += " "*10+f"{product.uniprot_id} - {product.scores['adv_score']:.2f} - {product.description}".center(100)+"\n"
                string += tabulate(create_go_table(product), headers="firstrow", tablefmt="grid", maxcolwidths=100) + "\n\n"

            # Add a separator and add the details for the bottom products
            string += ("-"*30)+"\n\n"
            for product in bottom_products:
                string += " "*10+f"{product.uniprot_id} - {product.scores['adv_score']:.2f} - {product.description}".center(100)+"\n"
                string += tabulate(create_go_table(product), headers="firstrow", tablefmt="grid", maxcolwidths=100) + "\n\n"

        return string

    def _generate_top_miRNAs_summary(self) -> str:
        # Create the header string.
        string = f"TOP {self.top_n} miRNAs" + "\n"
        string += f"+ annotates Product which is in top {self.top_n}, and - annotates Product which is in bottom {self.top_n}\n\n" 

        # Get the top and bottom products based on the advanced score
        sorted_products = sorted(self.reverse_lookup.products, key=lambda x: x.scores["adv_score"], reverse=True)
        top_products = sorted_products[:self.top_n]
        bottom_products = sorted_products[-self.top_n:]

        # Get the top miRNAs.
        top_miRNAs = sorted(self.reverse_lookup.miRNAs, key=lambda x: x.scores["basic_score"], reverse=True)[:self.top_n]

        # If verbosity is set to 1, create a table with the top miRNAs and their scores.
        if self.verbosity == 1:
            table = [["miRNA", "score"]]
            for _miRNA in top_miRNAs:
                table.append([_miRNA.id, _miRNA.scores["basic_score"]])
            string += tabulate(table, headers="firstrow", tablefmt="grid") + "\n\n"

        # If verbosity is set to 2 or higher, create a table with the top miRNAs, their scores, and the products they inhibit.
        if self.verbosity >= 2:
            table = [["miRNA", "score", "suppressed products"]]

            for _miRNA in top_miRNAs:
                inhibited_product_id = []
                for product_id, overlap in _miRNA.mRNA_overlaps.items():
                    if overlap >= self.reverse_lookup.miRNA_overlap_treshold:
                        inhibited_product_id.append(product_id)
                
                temp_list=[]
                temp_t_list = []
                temp_b_list = []
                for product_id in inhibited_product_id:
                    if any(product_id in sub.uniprot_id for sub in top_products):
                        temp_t_list.append(f"+{product_id}")
                    if any(product_id in sub.uniprot_id for sub in bottom_products):
                        temp_b_list.append(f"-{product_id}")
                    else:
                        temp_list.append(f"{product_id}")
                temp_list = temp_t_list + temp_list #To nsure the top/bottom products are displyed on top.
                temp_list = temp_b_list + temp_list
                table.append([_miRNA.id, _miRNA.scores["basic_score"], "\n".join(item for item in temp_list)])
            string += tabulate(table, headers="firstrow", tablefmt="grid") + "\n\n"

        # Return the summary string.
        return string

    def general_report(self, filepath:str):
        """
        Generates the general report and writes it to a file.

        Args:
        filepath (str): The path to the output file.
        """
        # Generate header of the report
        report = self._generate_header()+"\n\n"

        # Generate section on GOTerms
        if len(self.reverse_lookup.goterms) > 0:
            report += self._generate_section("GO TERMS")
            report += self._generate_goterms_statistics() + "\n" 
            report += self._generate_goterm_per_process_table() + "\n" 
        
        # Generate section on Products
        if len(self.reverse_lookup.products) > 0:
            report += self._generate_section("PRODUCTS")
            report += self._generate_top_bottom_products_summary() + "\n" 

        # Generate section on miRNAs
        if len(self.reverse_lookup.miRNAs) > 0:
            report += self._generate_section("miRNAs")
            report += self._generate_top_miRNAs_summary() + "\n" 

        # Create directory for the report file, if it does not exist
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        
        # Write the report to the output file
        with open(filepath, 'w') as f:
            f.write(report)
