from typing import Set, List
import requests
from requests.adapters import HTTPAdapter, Retry
import urllib.request
import gzip
import time
import os
from json import JSONDecodeError
#from typing import List, Dict, Set, Optional, Callable
from tqdm import trange, tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
from .FileUtil import FileUtil

import logging

logger = logging.getLogger(__name__)

class GOApi:
    """
    This class enables the user to interact with the Gene Ontology database via http requests.
    """
    def __init__(self):
        # Set up a retrying session
        retry_strategy = Retry(
            total=3,
            status_forcelist=[429, 500, 502, 503, 504],
            backoff_factor=0.3
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session = requests.Session()
        session.mount("http://", adapter)
        session.mount("https://", adapter)
        self.s = session

    def get_data(self, term_id):
        """
        Fetches term data for a given term ID from the Gene Ontology API using http://api.geneontology.org/api/ontology/term/{term_id}, 
        example of a term_id is GO:1903589.

        Returns:
          - (string as json) data: a json string, representing the api request response
        """
        url = f"http://api.geneontology.org/api/ontology/term/{term_id}"
        params = {}
        try:
            response = self.s.get(url, params=params, timeout=5)
            if response.ok:
                data = response.json()
                return data
            else:
                logger.warning(f"Error: {response.status_code} - {response.reason}")
        except requests.exceptions.RequestException as e:
            logger.warning(f"Error: {e}")
            return None

    def get_products(self, term_id):
        
        """
        Fetches product IDs associated with a given term ID from the Gene Ontology API. The product IDs can be of any of the following
        databases: UniProt, ZFIN, Xenbase, MGI, RGD [TODO: enable the user to specify databases himself]

        The request uses this link: http://api.geneontology.org/api/bioentity/function/{term_id}/genes

        Returns:
          - (string as json) data: a json string, representing the api request response
        """
        APPROVED_DATABASES = [["UniProtKB", ["NCBITaxon:9606"]],
                      ["ZFIN", ["NCBITaxon:7955"]],
                      #["RNAcentral", ["NCBITaxon:9606"]],
                      ["Xenbase", ["NCBITaxon:8364"]],
                      ["MGI", ["NCBITaxon:10090"]],
                      ["RGD", ["NCBITaxon:10116"]]]
        url = f"http://api.geneontology.org/api/bioentity/function/{term_id}/genes"
        params = {"rows": 10000000}
        products_set = set()

        max_retries = 5 # try api requests for max 5 times
        for i in range(max_retries):
            try:
                response = self.s.get(url, params=params, timeout=5)
                response.raise_for_status()

                json = response.json()
                for assoc in json['associations']:
                    if assoc['object']['id'] == term_id and any((database[0] in assoc['subject']['id'] and any(taxon in assoc['subject']['taxon']['id'] for taxon in database[1])) for database in APPROVED_DATABASES):
                        product_id = assoc['subject']['id']
                        products_set.add(product_id)
                products = list(products_set)
                logger.info(f"Fetched products for GO term {term_id}")
                return products
        
            except (requests.exceptions.RequestException, JSONDecodeError) as e :
                if i == (max_retries - 1): # this was the last http request, it failed
                    logger.error(f"Experienced an http exception or a JSONDecodeError while fetching products for {term_id}")
                    error_log_filepath = FileUtil.find_win_abs_filepath("log_output/error_log")
                    error_type = type(e).__name__
                    error_text = str(e)

                    logger.error(f"Exception type: {error_type}")
                    logger.error(f"Exception text: {error_text}")
                    logger.error(f"Debug report was written to: {error_log_filepath}")
                    
                    with open(error_log_filepath, "a+") as f:
                        f.write(f"Fetch products error for: {term_id}\n")
                        f.write(f"Exception: {error_type}\n")
                        f.write(f"Cause: {error_text}\n")
                        f.write(f"\n\n\n")
                        f.write(f"------------------------------\n")
                else:
                    time.sleep(500) # sleep 500ms before trying another http request
                return None



class GOAnnotiationsFile:
    """
    This class provides access to a Gene Ontology Annotations File, which stores the relations between each GO Term and it's products (genes),
    along with an evidence code, confirming the truth of the interaction. A GO Annotation comprises of a) GO Term, b) gene / gene product c) evidence code.
    
    See also:
      - http://geneontology.org/docs/download-go-annotations/ 
      - http://current.geneontology.org/products/pages/downloads.html
    """
    def __init__(self) -> None:
        self._filepath = "src_data_files/goa_human.gaf"
        self._check_file()
        if self._check_file():
            with open(self._filepath, 'r') as read_content:
                temp_content = read_content.readlines()
                self._readlines = []
                for line in temp_content:
                    if not line.startswith('!') and not line.strip() == '':
                        self._readlines.append(line.strip())
            
    def _check_file(self):
        os.makedirs(os.path.dirname(self._filepath), exist_ok=True)
        if os.path.exists(self._filepath):
            return True
        else:
            url = "http://geneontology.org/gene-associations/goa_human.gaf.gz"
            # download the gzip file and save it to a temporary file
            temp_file, _ = urllib.request.urlretrieve(url)

            # read the contents of the gzip file and save it to the txt file
            with gzip.open(temp_file, 'rt') as f_in, open(self._filepath, 'w') as f_out:
                for line in f_in:
                    f_out.write(line)
                    
            # delete the temporary file
            os.remove(temp_file)

        if os.path.exists(self._filepath):
            return True
        else:
            return False

    def get_products(self, goterm_id: str) -> List[str]:
        """This method returns all unique products associated with the GO term id

        Args:
            goterm_id (str): _description_

        Raises:
            ValueError: _description_
            Exception: _description_

        Returns:
            Set[str]: set of products' gene names
        """
        products_set = set()
        for line in self._readlines:
            chunks = line.split('\t')
            if goterm_id == chunks[4]:
                products_set.add(chunks[2])
        return list(products_set)
                
    def get_all_terms_for_product(self, product: str) -> List[str]:
        terms_set = set()
        for line in self._readlines:
            chunks = line.split('\t')
            if product == chunks[2]:
                terms_set.add(chunks[4])
        return list(terms_set)
        
        
    def get_all_terms(self) -> List[str]:
        terms_set = set()
        for line in self._readlines:
            chunks = line.split('\t')
            terms_set.add(chunks[4])
        return list(terms_set)
        
class UniProtAPI:
    """
    This class enables the user to interact with the UniProtKB database via http requests.
    """
    def __init__(self):
        # Set up a retrying session
        retry_strategy = Retry(
            total=3,
            status_forcelist=[429, 500, 502, 503, 504],
            backoff_factor=0.3
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session = requests.Session()
        session.mount("http://", adapter)
        session.mount("https://", adapter)
        self.s = session
    
    def get_uniprot_id(self, gene_name):
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
        url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene_name}+AND+organism_id:9606&format=json&fields=accession,gene_names,organism_name,reviewed,xref_ensembl"

        # Try the request up to `retries` times
        try:
            # Make the request and raise an exception if the response status is not 200 OK
            response = self.s.get(url, timeout=5)
            response.raise_for_status()
        except requests.exceptions.RequestException:
            # if there was an error with the HTTP request, log a warning
            logger.warning(f"Failed to fetch UniProt data for {gene_name}")
            return None

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
            logger.warning(f"Auto-selectd first reviewed result for {gene_name}!")
            return "UniProtKB:" + uniprot_id
        else:
            # raise an error if the user's choice is not valid
            raise ValueError(f"Invalid choice: {choice}")
  
    def get_uniprot_info(self, uniprot_id: str) -> dict:
        """
        Given a UniProt ID, returns a dictionary containing various information about the corresponding protein using the UniProt API.
        """
        def _return_mane_select_values_from_uniprot_query(result: dict) -> tuple:
            """
            Given the UniProt search result dictionary, return Ensembl gene ID, Ensembl transcript ID, and RefSeq nucleotide ID for the MANE-select transcript.
            """
            mane_indices = [index for (index, d) in enumerate(result["uniProtKBCrossReferences"]) if d["database"] == "MANE-Select"]
            if len(mane_indices) == 1:
                i = mane_indices[0]
                enst_id = result["uniProtKBCrossReferences"][i]["id"]
                refseq_nt_id = next((entry["value"] for entry in result["uniProtKBCrossReferences"][i]["properties"] if entry["key"] == "RefSeqNucleotideId"), None)
                ensg_id = next((next((sub["value"] for sub in entry["properties"] if sub["key"] == "GeneId"), None) for entry in result["uniProtKBCrossReferences"] if (entry["database"] == "Ensembl" and entry["id"] == enst_id)), None)
            else:
                return None, None, None
            return ensg_id, enst_id, refseq_nt_id

        # Extract UniProt ID if given in "database:identifier" format
        if ":" in uniprot_id:
            uniprot_id = uniprot_id.split(":")[1]

        # Construct UniProt API query URL
        url = f"https://rest.uniprot.org/uniprotkb/search?query={uniprot_id}+AND+organism_id:9606&format=json&fields=accession,gene_names,organism_name,reviewed,xref_ensembl,xref_refseq,xref_mane-select,protein_name"

        try:
            response = requests.get(url, timeout=5)
            response.raise_for_status()
        except requests.exceptions.RequestException:
            logger.warning(f"Failed to fetch UniProt data for {uniprot_id}")
            return {}
        
        results = response.json()["results"]
        if len(results) == 0:
            return {}
        else:
            # Get values from the UniProt search result
            result = next((entry for entry in results if entry["primaryAccession"] == uniprot_id), None)
            name = result["genes"][0]["geneName"]["value"]
            if "proteinDescription" in result and "recommendedName" in result["proteinDescription"] and "fullName" in result["proteinDescription"]["recommendedName"] and "value" in result["proteinDescription"]["recommendedName"]["fullName"]:
                description = result["proteinDescription"]["recommendedName"]["fullName"]["value"]
            else:
                description = "ERROR: Couldn't fetch description."
                logger.warning(f"proteinDescription, recommendedName, fullName or value not found when querying for uniprot info for the id: {uniprot_id}")
                logger.warning(f"result: {result}")
            ensg_id, enst_id, refseq_nt_id = _return_mane_select_values_from_uniprot_query(result)
            return {"genename": name, "description": description, "ensg_id": ensg_id, "enst_id": enst_id, "refseq_nt_id": refseq_nt_id}

class EnsemblAPI:
    def __init__(self):
        # Set up a retrying session
        retry_strategy = Retry(
            total=3,
            status_forcelist=[429, 500, 502, 503, 504],
            backoff_factor=0.3
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session = requests.Session()
        session.mount("http://", adapter)
        session.mount("https://", adapter)
        self.s = session

    def get_human_ortholog(self, id):
        """
        Given an source ID, detect organism and returns the corresponding human ortholog using the Ensembl API.
        """
        if "ZFIN" in id:
            species = "zebrafish"
            id_url = id.split(":")[1]
        elif "Xenbase" in id:
            species  = "xenopus_tropicalis"
            id_url = id.split(":")[1]
        elif "MGI" in id:
            species  = "mouse"
            id_url = id
        elif "RGD" in id:
            species  = "rat"
            id_url = id.split(":")[1]
        else:
            logger.info(f"No predefined organism found for {id}")
            return None
        url = f"https://rest.ensembl.org/homology/symbol/{species}/{id_url}?target_species=human;type=orthologues;sequence=none"
        try:
            response = self.s.get(url, headers={"Content-Type": "application/json"}, timeout=5)
            response.raise_for_status()
        except requests.exceptions.RequestException:
            return None
                
        response_json = response.json()["data"][0]["homologies"]
        if response_json == []:
            return None
        best_ortholog_dict = max(response_json, key=lambda x: int(x["target"]["perc_id"]))
        ortholog = best_ortholog_dict["target"].get("id")
        logger.info(f"Received ortholog for id {id} -> {ortholog}")
        return ortholog

    def get_sequence(self, ensembl_id, sequence_type="cdna"):
        """
        Given an Ensembl ID, returns the corresponding nucleotide sequence using the Ensembl API.
        """
        url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?object_type=transcript;type={sequence_type}"
        try:
            response = self.s.get(url, headers={"Content-Type": "text/plain"}, timeout=5)
            response.raise_for_status()
        except requests.exceptions.RequestException:
            logger.warning(f"Failed to fetch Ensembl sequence for {ensembl_id}")
            return None
        sequence = response.text
        logger.info(f"Received sequence for id {ensembl_id}.")
        return sequence
    
    def get_info(self, id: str) -> dict:
        """Can receive Ensembl id or symbol (human)

        Args:
            id (str): Ensembl ID or symbol

        Returns:
            dict: Information about the gene
        """
        species_mapping = {
            "ZFIN": "zebrafish/",
            "Xenbase": "xenopus_tropicalis/",
            "MGI": "mouse/MGI:",
            "RGD": "rat/",
            "UniProtKB": "human/",
        }

        # Check if the ID is an Ensembl ID or symbol
        if id.startswith("ENS"):
            endpoint = f"id/{id}"
        else:
            prefix, id_ = id.split(":") if ":" in id else (None, id)
            species = species_mapping.get(prefix, "human/") #defaults to human if not prefix "xxx:"
            endpoint = f"symbol/{species}{id_}"

        try:
            response = self.s.get(
                f"https://rest.ensembl.org/lookup/{endpoint}?mane=1;expand=1",
                headers={"Content-Type": "application/json"},
                timeout=5,
            )
            response.raise_for_status()
            response_json = response.json()
        except requests.exceptions.RequestException:
            # If the request fails, try the xrefs URL instead
            try:
                response = self.s.get(
                    f"https://rest.ensembl.org/xrefs/{endpoint}?",
                    headers={"Content-Type": "application/json"},
                    timeout=5,
                )
                response.raise_for_status()
                response_json = response.json()
                # Use the first ENS ID in the xrefs response to make a new lookup request
                ensembl_id = next((xref["id"] for xref in response_json if "ENS" in xref["id"]), None)
                if ensembl_id:
                    response = self.s.get(
                        f"https://rest.ensembl.org/lookup/id/{ensembl_id}?mane=1;expand=1",
                        headers={"Content-Type": "application/json"},
                        timeout=5,
                    )
                    response.raise_for_status()
                    response_json = response.json()
                else:
                    raise Exception("no ensembl id returned")
            except Exception as e:
                logger.warning(f"Failed to fetch Ensembl info for {id}.")
                return {}
            

        # Extract gene information from API response
        ensg_id = response_json.get("id")
        name = response_json.get("display_name")
        description = response_json.get("description", "").split(" [")[0]
        
        canonical_transcript_id = next((entry.get("id") for entry in response_json["Transcript"] if entry.get("is_canonical")), None)
        mane_transcripts = [d for d in response_json["Transcript"] if d.get("MANE")]
        if len(mane_transcripts) == 0:
            ensembl_transcript_id = canonical_transcript_id
            refseq_id = None
        elif len(mane_transcripts) == 1:
            ensembl_transcript_id = mane_transcripts[0]["MANE"][0].get("id")
            refseq_id = mane_transcripts[0]["MANE"][0].get("refseq_match")
        else:
            selected_entry = next((entry for entry in mane_transcripts if entry.get("is_canonical")), None)
            if not selected_entry:
                ensembl_transcript_id = selected_entry["MANE"][0].get("id")
                refseq_id = selected_entry["MANE"][0].get("refseq_match")
            else:
                ensembl_transcript_id = mane_transcripts[0]["MANE"][0].get("id")  # select the first canonical transcript with MANE
                refseq_id = mane_transcripts[0]["MANE"][0].get("refseq_match")
                logger.warning(f"Found non-canonical MANE transcript for {id}")

        if ensembl_transcript_id:
            try:
                response = self.s.get(
                    f"https://rest.ensembl.org/xrefs/id/{ensembl_transcript_id}?all_levels=1;external_db=UniProt%",
                    headers={"Content-Type": "application/json"},
                    timeout=5,
                )
                response.raise_for_status()
                response_json = response.json()
            except requests.exceptions.RequestException:
                pass
            uniprot_id = next((entry.get("primary_id") for entry in response_json if entry.get("dbname") =="Uniprot/SWISSPROT"), None)

        logger.debug(f"Received info data for id {id}.")
        return {
            "ensg_id": ensg_id,
            "genename": name,
            "description": description,
            "enst_id": ensembl_transcript_id,
            "refseq_nt_id": refseq_id,
            "uniprot_id": uniprot_id,
        }

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
            if len(result_list) >= 4: # bugfix
                return result_list[3]
            else:
                logger.warning(f"FAULTY LINE IN RGD, linesplit =: {linesplit}")
                return None

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
