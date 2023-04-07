import json
import requests
import time
import os
from typing import List, Dict, Set, Optional
from tqdm import trange, tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
from prettytable import PrettyTable
from texttable import Texttable

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
    
        uniprot_id = uniprot_id.split(":")[1]
        url = f"{self.base_url}search?query={uniprot_id}+AND+organism_id:9606&format=json&fields=accession,gene_names,organism_name,reviewed,xref_ensembl,protein_name"
        for i in range(retries):
            try:
                response = requests.get(url, timeout=timeout)
                response.raise_for_status()
                results = response.json()["results"]
                if len(results) == 0:
                    return {}
                else:
                    full_name = results[0]["proteinDescription"]["recommendedName"]["fullName"]["value"]
                    ensembl_id = _return_ensembl_from_id_and_uniprot_query(uniprot_id, results)
                    return {"full_name" : full_name, "ensembl_id" : ensembl_id}
            except requests.exceptions.RequestException:
                logger.warning(f"Failed to fetch UniProt data for {uniprot_id}")
            time.sleep(timeout)
        return {}

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
    def __init__(self, id: str, uniprot_id: str = None, description: str = None, ensembl_id: str= None, mRNA: str = None, scores: dict = None):
        """
        A class representing a product (e.g. a gene or protein).

        Args:
            id (str): The ID of the product.
            uniprot_id (str): The UniProt ID of the product.
            description (str): A description of the product.
            ensembl_id (str):
            mRNA (str): The mRNA sequence of the product.
            scores (dict): A dictionary of scores associated with the product (e.g. expression score, functional score).
        """
        self.id = id
        self.uniprot_id = uniprot_id
        self.description = description
        self.ensembl_id = ensembl_id
        self.mRNA = mRNA
        self.scores = {} if scores is None else scores.copy()
    
    def fetch_UniprotID(self, human_ortolog_finder: HumanOrthologFinder, uniprot_api: UniProtAPI) -> None:
        if 'UniProtKB' in self.id:
            self.uniprot_id = self.id
        else:
            human_ortholog_gene_id = human_ortolog_finder.find_human_ortholog(self.id)
            if human_ortholog_gene_id is not None:
                uniprot_id = uniprot_api.get_uniprot_id(human_ortholog_gene_id)
                if uniprot_id is not None:
                    self.uniprot_id = uniprot_id

    def fetch_Uniprot_info(self, uniprot_api: UniProtAPI) -> None:
        if self.uniprot_id == None: return
        info_dict = uniprot_api.get_uniprot_info(self.uniprot_id)
        full_name = info_dict.get("full_name")
        ensembl_id = info_dict.get("ensembl_id")
        if full_name is not None:
            self.description = full_name
        if ensembl_id is not None:
            self.ensembl_id = ensembl_id
                

    @classmethod
    def from_dict(cls, d: dict) -> 'Product':
        """
        Class method to create a new Product instance from a dictionary.

        Args:
            d (dict): The dictionary containing the data to create the Product instance.

        Returns:
            Product: A new Product instance created from the input dictionary.
        """
        return cls(d['id'], d.get('uniprot_id'), d.get('description'), d.get('ensembl_id'), d.get('mRNA'), d.get('scores'))

class ReverseLookup:
    def __init__(self, goterms: List[GOTerm], target_processes: List[Dict[str, str]], products: List[Product] = []):
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

    def fetch_all_go_term_names_descriptions(self):
        """
        Iterates over all GOTerm objects in the go_term set and calls the fetch_name_description method for each object.
        """
        api = GOApi()
        with logging_redirect_tqdm():
            for goterm in tqdm(self.goterms):
                goterm.fetch_name_description(api)

    def fetch_all_go_term_products(self):
        """
        Iterates over all GOTerm objects in the go_term set and calls the fetch_products method for each object.
        """
        api = GOApi()
        with logging_redirect_tqdm():
            for goterm in tqdm(self.goterms):
                goterm.fetch_products(api)

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
        # Create an empty set to store unique products
        products_set = set()
        
        # Iterate over each GOTerm object in the go_term set and retrieve the set of products associated with that GOTerm
        # object. Add these products to the products_set.
        for term in self.goterms:
            products_set.update(term.products)
        
        # Iterate over each product in the products_set and create a new Product object from the product ID using the
        # Product.from_dict() classmethod. Add the resulting Product objects to the ReverseLookup object's products list.
        for product in products_set:
            self.products.append(Product.from_dict({'id':product}))
        logger.info(f"Created Product objects from GOTerm object definitions")

    def fetch_UniprotID_products(self) -> None:
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

    def fetch_Uniprot_infos(self) -> None:
        try:   
            uniprot_api = UniProtAPI()
            # Iterate over each Product object in the ReverseLookup object.
            with logging_redirect_tqdm():
                for product in tqdm(self.products):
                    # Check if the Product object doesn't have a UniProt ID.
                    if product.description == None and product.uniprot_id is not None:
                        # If it doesn't, fetch UniProt data for the Product object.
                        product.fetch_Uniprot_info(uniprot_api)
        except Exception as e:
            raise e 

    def score_products(self) -> None:
        def _opposite_direction(direction):
            if direction == "0":
                return "0"
            elif direction == "+":
                return "-"
            elif direction == "-":
                return "+"

        def adv_score(goterms_list: List[GOTerm]) -> float:
            score = 0.0

            a = 10
            b1 = 2
            b2 = 0.5
            c1 = 1
            c2 = 0.1

            # Check if all target processes are regulated in the same direction as the GOTerms in the list
            # and none of them are regulated in the opposite direction
            if (
                # Check if all processes in target_processes have a GOTerm in goterms_list that regulates it in the same direction
                all(
                    any(process['direction'] == goterm.direction and process['process'] == goterm.process for goterm in goterms_list)
                    for process in self.target_processes
                )
                # Check if none of the processes in target_processes have a GOTerm in goterms_list that regulates it in the opposite direction
                and not any(
                    any(_opposite_direction(process['direction']) == goterm.direction and process['process'] == goterm.process for goterm in goterms_list)
                    for process in self.target_processes
                )
            ):
                # If all target processes are regulated in the same direction, add a points to the score
                score += a

            # Check if all target processes are regulated in the opposite direction as the GOTerms in the list
            # and none of them are regulated in the same direction
            if (
                # Check if all processes in target_processes have a GOTerm in goterms_list that regulates it in the opposite direction
                all(
                    any(_opposite_direction(process['direction']) == goterm.direction and process['process'] == goterm.process for goterm in goterms_list)
                    for process in self.target_processes
                )
                # Check if none of the processes in target_processes have a GOTerm in goterms_list that regulates it in the same direction
                and not any(
                    any(process['direction'] == goterm.direction and process['process'] == goterm.process for goterm in goterms_list)
                    for process in self.target_processes
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
                    for process in self.target_processes
            )
            # Calculate the score based on the number of processes in target_processes that are regulated
            # by GOTerms in the oposite direction as defined in the list
            score -= sum(
                (b1**(b2 * sum(
                    # Check if the direction and process in the process dict matches with the direction and process in any GOTerm dict
                    goterm.weight for goterm in goterms_list if _opposite_direction(process['direction']) == goterm.direction and process['process'] == goterm.process)))
                    for process in self.target_processes
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

        def nterms_per_process(goterms_list: List[GOTerm]) -> dict:
            # Create an empty dictionary to store the count of GOTerms for each process and direction
            nterms_dict = {}
            
            # Iterate over each process in the target_processes list
            for process in self.target_processes:
                # Count the number of GOTerms that have a direction of "+" and a process matching the current process
                nterms_dict[f"{process['process']}+"] = sum(1 for goterm in goterms_list if (goterm.direction == "+" and process['process'] == goterm.process))
                
                # Count the number of GOTerms that have a direction of "-" and a process matching the current process
                nterms_dict[f"{process['process']}-"] = sum(1 for goterm in goterms_list if (goterm.direction == "-" and process['process'] == goterm.process))
                
                # Count the number of GOTerms that have a direction of "0" and a process matching the current process
                nterms_dict[f"{process['process']}0"] = sum(1 for goterm in goterms_list if (goterm.direction == "0" and process['process'] == goterm.process))
            
            # Return the dictionary containing the count of GOTerms for each process and direction
            return nterms_dict

        for product in self.products:
            product_id = product.id
            goterms_list = []
            for goterm in self.goterms:
                if product_id in goterm.products:
                    goterms_list.append(goterm)
            
            adv_score_result = adv_score(goterms_list)
            dict_of_terms_per_process = nterms_per_process(goterms_list)

            product.scores = {"adv_score":adv_score_result, "nterms":dict_of_terms_per_process}

    def load_go_term_datafile(self, filename: str) -> None:
        with open(filename, 'r') as f:
            data = json.load(f)
        for element in data['goterms']:
            for term in self.goterms:
                if term.id in element.values():
                    term.name = element.get('name')
                    term.description = element.get('description')
                    term.products = element.get('products')
        logger.info(f"Loaded goterms from {filename}")

    def save_goterms_to_datafile(self, filename: str) -> None:
        """
        Saves all GOTerm objects to a JSON file.

        Args:
            filename (str): The name of the file to save the data to.
        """
        data = {'goterms': []}
        for goterm in self.goterms:
            data['goterms'].append(goterm.__dict__)
        
        with open(filename, 'w') as f:
            json.dump(data, f, indent=4)

    def load_products_datafile(self, filename: str) -> None:
        with open(filename, 'r') as f:
            data = json.load(f)
        for element in data['products']:
            for product in self.products:
                if product.id in element.values():
                    product.uniprot_id = element.get('uniprot_id')
                    product.description = element.get('description')
                    product.ensembl_id = element.get('ensembl_id')
                    product.mRNA = element.get('mRNA')
                    product.scores = element.get('scores')

    def save_products_to_datafile(self, filename: str) -> None:
        """
        Saves all Product objects to a JSON file.

        Args:
            filename (str): The name of the file to save the data to.
        """
        data = {'products': []}
        for product in self.products:
            data['products'].append(product.__dict__)
        
        with open(filename, 'w') as f:
            json.dump(data, f, indent=4)
    
    @classmethod
    def from_file(cls, filepath: str) -> 'ReverseLookup':
        """
        Creates a ReverseLookup object from a JSON file.

        Args:
            filepath (str): The path to the input file.

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
                
        with open(filepath, "r") as read_content:
            target_processes = []
            go_terms = []
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
        return cls(go_terms, target_processes)
    
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

class ReportGenerator:
    def __init__(self, reverse_lookup: ReverseLookup):
        self.reverse_lookup = reverse_lookup
    
    def _generate_header(self, text: str) -> str:
        header = f"{'=' * 40}\n{text}\n{'=' * 40}\n"
        return header
    
    def _generate_section(self, text: str) -> str:
        header = f"{'-' * 40}\n{text}\n{'-' * 40}\n"
        return header

    def generate_detailed_design_report(self, filepath: str):
        # Group the GOTerms by process and direction
        grouped_goterms = {}
        for target in self.reverse_lookup.target_processes:
            grouped_goterms[(target["process"], "+")] = []
            grouped_goterms[(target["process"], "-")] = []
            grouped_goterms[(target["process"], "0")] = []
        for goterm in self.reverse_lookup.goterms:
            key = (goterm.process, goterm.direction)
            grouped_goterms[key].append(goterm)

        # Generate the report header
        report = self._generate_header("Detailed design report")

        # Generate the Target Processes section of the report
        target_processes_table = PrettyTable()
        target_processes_table.field_names = ["Process", "Direction"]
        for target in self.reverse_lookup.target_processes:
            target_processes_table.add_row([target["process"], target["direction"]])
        report += "\nTarget Processes\n\n" + str(target_processes_table)

        # Generate the GOTerm section of the report
        goterm_table = PrettyTable()
        goterm_table.field_names = ["Process", "Direction", "GOTerm ID", "GOTerm Name"]
        for (process, direction), goterms in grouped_goterms.items():
            goterm_table.add_row([process, direction, "", ""])
            for goterm in goterms:
                goterm_table.add_row(["", "", goterm.id, goterm.name])
        report += "\n\nGO Terms\n\n" + str(goterm_table)

        # Generate the Product section of the report
        product_table = PrettyTable()
        product_table.field_names = ["Product Name", "Product Description", "Involved in Processes"]
        for product in self.reverse_lookup.products:
            processes = ", ".join([goterm.id for goterm in self.reverse_lookup.goterms if product.id in goterm.products])
            product_table.add_row([product.id, product.description, processes])
        report += "\n\nProducts\n\n" + str(product_table)

        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        with open(filepath, 'w') as f:
            f.write(report)

    def generate_summary_report(self, filepath=None):
            
        report = ""

        # Generate header
        report += self._generate_header("MODEL SUMMARY")
        
        # Generate the Target Processes section of the report
        target_processes_table = Texttable()
        target_processes_table.header(["Process", "Direction"])
        for target in self.reverse_lookup.target_processes:
            target_processes_table.add_row([target["process"], target["direction"]])
        report += "\nTarget Processes\n" + target_processes_table.draw() + '\n\n'
        
        # Generate GOTerms summary
        goterms = self.reverse_lookup.goterms
        goterms_table = Texttable()
        goterms_table.header(["Process", "+", "-", "0", "Total"])
        total_p = 0
        total_m = 0
        total_z = 0
        total = 0
        for process in self.reverse_lookup.target_processes:
            process_goterms = [goterm for goterm in goterms if goterm.process == process['process']]
            n_plus = len([goterm for goterm in process_goterms if goterm.direction == '+'])
            n_minus = len([goterm for goterm in process_goterms if goterm.direction == '-'])
            n_zero = len([goterm for goterm in process_goterms if goterm.direction == '0'])
            n_total = n_plus + n_minus + n_zero
            total_p += n_plus
            total_m += n_minus
            total_z += n_zero
            total += n_total
            goterms_table.add_row([f"{process['process']}", n_plus, n_minus, n_zero, n_total])
        goterms_table.add_row(['Total', total_p, total_m, total_z, total])

        report += "GOTerms Summary\n"
        report += goterms_table.draw() + '\n\n'
        
        # Products
        histogram_products_table = Texttable()
        histogram_products_table.header(["Involved in N GOTerms", "N"])
        num_goterms = []
        for product in self.reverse_lookup.products:
            num_goterms.append(len([goterm for goterm in self.reverse_lookup.goterms if product.id in goterm.products]))
        for i in range (1,max(num_goterms)+1):
            histogram_products_table.add_row([i, num_goterms.count(i)])
        histogram_products_table.add_row(["Total", len(self.reverse_lookup.products)])
        report += "\nProduct distribution\n" + histogram_products_table.draw() + '\n\n'

        # Save to file
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        with open(filepath, 'w') as f:
            f.write(report)