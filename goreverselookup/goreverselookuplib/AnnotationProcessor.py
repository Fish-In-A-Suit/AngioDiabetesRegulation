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
import aiohttp, asyncio

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

    def get_data(self, term_id, get_url_only = False):
        """
        Fetches term data for a given term ID from the Gene Ontology API using http://api.geneontology.org/api/ontology/term/{term_id}, 
        example of a term_id is GO:1903589.

        If get_url_only == True, this will only return the url.

        Returns:
          - (string as json) data: a json string, representing the api request response
        """
        url = f"http://api.geneontology.org/api/ontology/term/{term_id}"
        params = {}
        if get_url_only:
            return url
        logger.debug(f"Querying: {url}")
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

    def get_products(self, term_id, get_url_only=False, request_params = {"rows": 10000000}): 
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
        params = request_params
        
        # used in async
        if get_url_only == True:
            # create a request object with the base url and params
            request = requests.Request("GET", url, params=params)
            # prepare the request
            prepared_request = self.s.prepare_request(request)
            # get the fully constructed url with parameters
            url = prepared_request.url
            return url

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
                    #time.sleep(500) # sleep 500ms before trying another http request
                    time.sleep(0.5) # time.sleep is in SECONDS !!!
                return None
    
    async def get_products_async(self, term_id):
        """
        Fetches product IDs associated with a given term ID from the Gene Ontology API. The product IDs can be of any of the following
        databases: UniProt, ZFIN, Xenbase, MGI, RGD [TODO: enable the user to specify databases himself]

        This function works asynchronously, much faster than it's synchronous 'get_products' counterpart.

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
        MAX_RETRIES = 5
        url = f"http://api.geneontology.org/api/bioentity/function/{term_id}/genes"
        params = {"rows": 100000}

        global request_iterations
        request_iterations = 0 # global variable request_iterations to keep track of the amount of requests submitted to the server (maximum is MAX_RETRIES); a harsh bugfix
        
        # create a request object with the base url and params
        #request = requests.Request("GET", url, params=params)
        # prepare the request
        #prepared_request = self.s.prepare_request(request)
        # get the fully constructed url with parameters
        #url = prepared_request.url

        """ # this code caused timeout errors
        async with aiohttp.ClientSession() as session:            
            response = await session.get(url, params=params, timeout=5)
            # DEBUG: if the above doesn't work, uncomment the request pre-process lines above and use session.get(url) only
            if response.status == 200:
                data = await response.json()
                products_set = set()

                for assoc in data['associations']:
                    if assoc['object']['id'] == term_id and any((database[0] in assoc['subject']['id'] and any(taxon in assoc['subject']['taxon']['id'] for taxon in database[1])) for database in APPROVED_DATABASES):
                        product_id = assoc['subject']['id']
                        products_set.add(product_id)

                products = list(products_set)
                logger.info(f"Fetched products for GO term {term_id}")
                return products
            else: # warn of the error; # TODO: implement error logging to file as in 'get_products' function
                logger.error(f"Couldn't fetch products for {term_id}")
        """
        #i=0
        # TODO: PRODUCT REQUESTS START TO CAUSE TIMEOUT ERRORS AFTER ~15 ASYNC REQUESTS -> BATCH? DIFFERENT IPS?
        # as per: https://stackoverflow.com/questions/51248714/aiohttp-client-exception-serverdisconnectederror-is-this-the-api-servers-issu
        connector = aiohttp.TCPConnector(limit=20) # default limit is 100
        async with aiohttp.ClientSession(connector=connector) as session:
            #for i in range(MAX_RETRIES):
            #while i < MAX_RETRIES: # due to the async nature, each iteration resets i; hence "i" is useless -> bugfix: global variable request_iterations
            while request_iterations < MAX_RETRIES:
                try:
                    request_iterations += 1
                    response = await session.get(url, params=params, timeout=7)
                    response.raise_for_status() # checks for anything other than status 200
                    data = await response.json()
                    products_set = set()
                    for assoc in data['associations']:
                        if assoc['object']['id'] == term_id and any((database[0] in assoc['subject']['id'] and any(taxon in assoc['subject']['taxon']['id'] for taxon in database[1])) for database in APPROVED_DATABASES):
                            product_id = assoc['subject']['id']
                            products_set.add(product_id)

                    products = list(products_set)
                    logger.info(f"Fetched products for GO term {term_id}")
                    request_iterations = 0 # reset
                    return products
                except (requests.exceptions.RequestException, JSONDecodeError, asyncio.exceptions.TimeoutError) as e :
                    # logger.error(f"TimoutError on retry attempt {request_iterations}. Exception: {e}")
                    # i += 1
                    # if i == (MAX_RETRIES - 1): # this was the last http request, it failed
                    # if request_iterations == (MAX_RETRIES - 1):
                    if request_iterations == MAX_RETRIES: # due to while loop logic we don't substract 1
                        error_log_filepath = FileUtil.find_win_abs_filepath("log_output/error_log")
                        error_type = type(e).__name__
                        error_text = str(e)

                        # logger.error(f"Exception type: {error_type}")
                        # logger.error(f"Exception text: {error_text}")
                        # logger.error(f"Debug report was written to: {error_log_filepath}")
                        logger.error(f"https error for {term_id}, error_type = {error_type}, error_text = {error_text}")

                        with open(error_log_filepath, "a+") as f:
                            f.write(f"Fetch products error for: {term_id}\n")
                            f.write(f"Exception: {error_type}\n")
                            f.write(f"Cause: {error_text}\n")
                            f.write(f"\n\n\n")
                            f.write(f"------------------------------\n")
                    else:
                        # time.sleep(0.5)
                        time.sleep(1) # maybe with 1s the server won't start to block?
            # reset
            request_iterations = 0
    
    async def get_products_async_notimeout(self, term_id):
        """
        A testing variant of get_products_async. Doesn't include timeout in the url request, no retries.
        """
        APPROVED_DATABASES = [["UniProtKB", ["NCBITaxon:9606"]],
                      ["ZFIN", ["NCBITaxon:7955"]],
                      #["RNAcentral", ["NCBITaxon:9606"]],
                      ["Xenbase", ["NCBITaxon:8364"]],
                      ["MGI", ["NCBITaxon:10090"]],
                      ["RGD", ["NCBITaxon:10116"]]]
        url = f"http://api.geneontology.org/api/bioentity/function/{term_id}/genes"
        params = {"rows": 20000} # 10k rows resulted in 56 mismatches for querying products for 200 goterms (compared to reference model, loaded from synchronous query data)
        # DELAY = 1 # 1 second delay between requests

        # as per: https://stackoverflow.com/questions/51248714/aiohttp-client-exception-serverdisconnectederror-is-this-the-api-servers-issu
        connector = aiohttp.TCPConnector(limit=20, limit_per_host=20) # default limit is 100
        # as per: https://stackoverflow.com/questions/64534844/python-asyncio-aiohttp-timeout; DOESNT WORK!
        # session_timeout =   aiohttp.ClientTimeout(total=None,sock_connect=10,sock_read=10) -> async with aiohttp.ClientSession(connector=connector, timeout=session_timeout) as session; 
        # https://github.com/aio-libs/aiohttp/issues/3187 -> 504 gateways are server-limited !

        ### POSSIBLE ERROR SOLUTION ### [TODO: continue from here]
        # Current algorithm creates one aiohttp.ClientSession FOR EACH GOTERM. Therefore, each ClientSession only has one connection,
        # and the checks for connection limiting aren't enforeced. During runtime, there can be as many as 200 (as many as there are goterms)
        # active ClientSessions, each with only one request. You should code in the following manner:
        #
        # async def make_requests():
        #    connector = aiohttp.TCPConnector(limit=20, limit_per_host=20)
        #    async with aiohttp.ClientSession(connector=connector) as session:
        #        urls = [...]  # List of URLs to request
        #        for url in urls:
        #            await asyncio.sleep(1)  # Introduce a 1-second delay between requests
        #            response = await session.get(url)
        #            # Process the response

        async with aiohttp.ClientSession(connector=connector) as session:
            response = await session.get(url, params=params)
            # response.raise_for_status() # checks for anything other than status 200
            if response.status != 200: # return HTTP Error if status is not 200 (not ok), parse it into goterm.http_errors -> TODO: recalculate products for goterms with http errors
                logger.warning(f"HTTP Error when parsing {term_id}. Response status = {response.status}")
                return f"HTTP Error: status = {response.status}, reason = {response.reason}"
                     
            data = await response.json()
            products_set = set()
            for assoc in data['associations']:
                if assoc['object']['id'] == term_id and any((database[0] in assoc['subject']['id'] and any(taxon in assoc['subject']['taxon']['id'] for taxon in database[1])) for database in APPROVED_DATABASES):
                    product_id = assoc['subject']['id']
                    products_set.add(product_id)

            products = list(products_set)
            logger.info(f"Fetched products for GO term {term_id}")
            return products
        
    """ THIS CAUSED CIRCULAR IMPORT DUE TO GOTerm -> was moved into the GOTerm class
    async def get_products_async_v3(self, goterm:GOTerm, session:aiohttp.ClientSession, request_params = {"rows":20000}, req_delay=0.5):
        
        #A testing variant of get_products_async. Doesn't include timeout in the url request, no retries.
        #Doesn't create own ClientSession, but relies on external ClientSession, hence doesn't overload the server as does the get_products_async_notimeout function.
        
        # Previous algorithm created one aiohttp.ClientSession FOR EACH GOTERM. Therefore, each ClientSession only had one connection,
        # and the checks for connection limiting weren't enforeced. During runtime, there could be as many as 200 (as many as there are goterms)
        # active ClientSessions, each with only one request. You should code in the following manner:
        #
        # async def make_requests():
        #    connector = aiohttp.TCPConnector(limit=20, limit_per_host=20)
        #    async with aiohttp.ClientSession(connector=connector) as session:
        #        urls = [...]  # List of URLs to request
        #        for url in urls:
        #            await asyncio.sleep(1)  # Introduce a 1-second delay between requests
        #            response = await session.get(url)
        #            # Process the response
        
        APPROVED_DATABASES = [["UniProtKB", ["NCBITaxon:9606"]],
                      ["ZFIN", ["NCBITaxon:7955"]],
                      #["RNAcentral", ["NCBITaxon:9606"]],
                      ["Xenbase", ["NCBITaxon:8364"]],
                      ["MGI", ["NCBITaxon:10090"]],
                      ["RGD", ["NCBITaxon:10116"]]]
        url = f"http://api.geneontology.org/api/bioentity/function/{goterm.id}/genes"
        params = request_params # 10k rows resulted in 56 mismatches for querying products for 200 goterms (compared to reference model, loaded from synchronous query data)
        
        asyncio.sleep(req_delay)
        response = await session.get(url, params=params)
        if response.status != 200: # return HTTP Error if status is not 200 (not ok), parse it into goterm.http_errors -> TODO: recalculate products for goterms with http errors
            logger.warning(f"HTTP Error when parsing {goterm.id}. Response status = {response.status}")
            return f"HTTP Error: status = {response.status}, reason = {response.reason}"
        
        data = await response.json()
        products_set = set()
        for assoc in data['associations']:
            if assoc['object']['id'] == goterm.id and any((database[0] in assoc['subject']['id'] and any(taxon in assoc['subject']['taxon']['id'] for taxon in database[1])) for database in APPROVED_DATABASES):
                product_id = assoc['subject']['id']
                products_set.add(product_id)

        products = list(products_set)
        goterm.products = products
        logger.info(f"Fetched products for GO term {goterm.id}")
        return products      
    """      

class GOAnnotiationsFile:
    def __init__(self, filepath:str="") -> None:
        """
        This class provides access to a Gene Ontology Annotations File, which stores the relations between each GO Term and it's products (genes),
        along with an evidence code, confirming the truth of the interaction. A GO Annotation comprises of a) GO Term, b) gene / gene product c) evidence code.

        Parameters:
          - (str) filepath: the filepath to the GO Annotations File downloaded file from http://current.geneontology.org/products/pages/downloads.html -> Homo Sapiens (EBI Gene Ontology Database) - protein = goa_human.gaf; link = http://geneontology.org/gene-associations/goa_human.gaf.gz
                            if left to default value, self._filepath will be set to 'src_data_files/goa_human.gaf'. The file should reside in root/src_data_files/ and the parameter filepath should be the file name of the downloaded file inside src_data_files/

        See also:
          - http://geneontology.org/docs/download-go-annotations/ 
          - http://current.geneontology.org/products/pages/downloads.html
        """
        if filepath == "":
            self._filepath = "src_data_files/goa_human.gaf"
        else:
            if "src_data_files/" not in filepath:
                self._filepath = f"src_data_files/{filepath}"
            else: # 'src_data_files/' is already in the filepath
                self._filepath = filepath

        self._check_file()
        if self._check_file():
            with open(self._filepath, 'r') as read_content:
                temp_content = read_content.readlines()
                self._readlines = []
                for line in temp_content:
                    if not line.startswith('!') and not line.strip() == '':
                        self._readlines.append(line.strip())
        self.terms_dict = None
        self.products_dict = None
            
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

    def get_all_products_for_goterm(self, goterm_id: str) -> List[str]:
        """
        This method returns all unique products associated with the GO term id

        Args:
            goterm_id (str): a GO Term identifier, eg. GO:0003723

        Returns:
            List[str]: a List of all product (gene/gene products) gene names, eg. ['NUDT4B', ...]
        
        Example: for 'GO:0003723' it returns ['KMT2C', 'CLNS1A', 'ZCCHC9', 'TERT', 'BUD23', 'DDX19B', 'CCAR2', 'NAP1L4', 'SAMSN1', 'ERVK-9', 'COA6', 'RTF1', 'AHCYL1', 'SMARCA4', ... (total len = 1378)]
        """
        
        if self.terms_dict is None:
            self.populate_terms_dict()
        
        return self.terms_dict.get(goterm_id, [])
               
    def populate_poducts_dict(self):
        """
        For each line in the readlines of the GO Annotations File, it creates a connection between a product gene name and an associated GO Term.

        The result is a dictionary (self.products_dict), mapping keys (product gene names, eg. NUDT4B) to a List of all associated
        GO Terms (eg. ['GO:0003723', ...])
        """
        self.products_dict = {}
        for line in self._readlines: # example line: 'UniProtKB \t A0A024RBG1 \t NUDT4B \t enables \t GO:0003723 \t GO_REF:0000043 \t IEA \t UniProtKB-KW:KW-0694 \t F \t Diphosphoinositol polyphosphate phosphohydrolase NUDT4B \t NUDT4B \t protein \t taxon:9606 \t 20230306 \t UniProt'
            chunks = line.split('\t')
            self.products_dict.setdefault(chunks[2], set()).add(chunks[4]) # create a key with the line's product gene name (if the key already exists, don't re-create the key - specified by the setdefault method) and add the associated GO Term to the value set. eg. {'NUDT4B': {'GO:0003723'}}, after first line is processed, {'NUDT4B': {'GO:0003723'}, 'NUDT4B': {'GO:0046872'}} after second line ...
        for key, values in self.products_dict.items(): # the set() above prevents the value elements (GO Terms) in dictionary to be repeated
            self.products_dict[key] = list(values) # converts the set to a List, eg. {'NUDT4B': ['GO:0003723']}     
            
    def populate_terms_dict(self):
        """
        For each line in the readlines of the GO Annotations File, it creates a connection between a GO Term and it's associated product gene name.

        The result is a dictionary (self.terms_dict), mapping keys (GO Terms, eg. GO:0003723) to a List of all
        associated product gene names (eg. ['NUDT4B', ...])
        """
        self.terms_dict = {}
        for line in self._readlines: # example line: 'UniProtKB \t A0A024RBG1 \t NUDT4B \t enables \t GO:0003723 \t GO_REF:0000043 \t IEA \t UniProtKB-KW:KW-0694 \t F \t Diphosphoinositol polyphosphate phosphohydrolase NUDT4B \t NUDT4B \t protein \t taxon:9606 \t 20230306 \t UniProt'
            chunks = line.split('\t')
            self.terms_dict.setdefault(chunks[4], set()).add(chunks[2]) # create a key with the line's GO Term (if the key already exists, don't re-create the key - specified by the setdefault method) and add the product' gene name to the value set. eg. {'GO:0003723': {'NUDT4B'}}, after first line is processed, {'GO:0003723': {'NUDT4B'}, 'GO:0046872': {'NUDT4B'}} after second line ...
        for key, values in self.terms_dict.items(): # the previous set() prevents the value elements (product gene names) in dictionary to be repeated
            self.terms_dict[key] = list(values) # converts the set to a List, eg. {'NUDT4B': ['GO:0003723']}

    def get_all_terms_for_product(self, product: str) -> List[str]:
        """
        Gets all GO Terms associated to a product gene name.

        Args:
          - (str) product: must be a gene name corresponding to a specific gene/gene product, eg. NUDT4B
        
        Returns:
          - List[str]: a List of all GO Term ids associated with the input product's gene name
        
        Example: for 'NUDT4B', it returns ['GO:1901911', 'GO:0071543', 'GO:0005737', 'GO:0000298', 'GO:0005634', 'GO:0034431', 'GO:0034432', 'GO:0046872', 'GO:0008486', 'GO:1901909', 'GO:0003723', 'GO:1901907', 'GO:0005829']
        """
        if self.products_dict is None:
            self.populate_poducts_dict()
        
        return self.products_dict.get(product, [])
        
        
    def get_all_terms(self) -> List[str]:
        """
        Returns a List of all unique GO Terms read from the GO Annotations file.
        In the current (27_04_2023) GO Annotation File, there are 18880 unique GO Terms.
        """
        if not self.terms_dict:
            self.populate_terms_dict()
        
        terms_list = [k for k,v in self.terms_dict.items()]
        return terms_list

    """ # THIS IS INVALID, orthologs cannot be queried from the GOAF !!!
    def get_all_product_orthologs(self, product_id:str):
        #""
        #Gets all orthologs in line for a specific product (gene) id. This function uses GOAF for the ortholog query.
        #TODO: check if this is actually even valid.
        #""
        # if a user sends a uniprotkb product_id here, default to get_uniprotkb_genename
        if "UniProtKB" in product_id:
            genename =  self.get_uniprotkb_genename(product_id)
            if genename != None:
                return genename
        
        possible_orthologs = {} # a dict between possible orthologs and the readlines where they are found
        for line in self._readlines:
            if product_id in line:
                line_elements = line.split("\t")
                if line_elements[1] != product_id:
                    # query the 8th line element With (or) From
                    # GOAF line elements: (1) DB (2) DB Object Id (3) DB Object Symbol (4) Qualifier (optional) (5) GO ID (6) DB:Reference (7) Evidence Code (8) With (or) from (optional) (9) Aspect ...
                    # We are interested in the 8th line element, but it is optional. Furthermore, Qualifier (4th line element) is also optional.
                    # However, we are certain that the "With or from" line element will appear after the Evidence Code (which is always a 3-character code - http://geneontology.org/docs/guide-go-evidence-codes/) and before the
                    # Aspect, which can be either "P" (biological process), "F" (biological function) or "C" (cellular component). If the difference between the evidence code index and the aspect index (index = position in the array) is
                    # greater than 1, then we are sure that the element between them is a "With or from" element.
                    evidence_code = ""
                    aspect = ""
                    evidence_code_index = -1
                    aspect_index = -1
                    i=0
                    for line_element in line_elements:
                        if len(line_element) == 3: # this is the Evidence Code
                            evidence_code = line_element
                            evidence_code_index = i
                        if len(line_element) == 1 and i > evidence_code_index: # this is the Aspect
                            aspect = line_element
                            aspect_index = i
                        i+=1
                    if aspect_index - evidence_code_index > 1:
                        # orthologs exist
                        orthologs = 
        """
    
    def get_uniprotkb_genename(self, product_id:str):
        """
        Gets the gene name (DB Object Symbol) for the supplied UniProtKB product_id.

        Parameters:
          - product_id: must be in the format UniProtKB:XXXXX 
        
        Algorithm:
            If the product_id is a UniProtKB, then the GOAF is browsed to obtain the gene name, which
            is the third line element (DB Object Symbol) in the GOAF. If the parse doesnt find the uniprot id
            (XXXXX) in the GOAF as the second line element, then the supplied UniProtKB may be an animal protein.
            In this case, the program also records any lines that don't have (XXXXX) as the second line element, but still
            contain the (XXXXX) in the line. The program reverts to these lines and attempts to find a human ortholog. If
            all lines result in the same human ortholog, then the search was successful. TODO: implement logic if multiple different
            orthologs are found.
        """
        if "UniProtKB" in product_id:
            product_id = product_id.split(":")[1] # gets the raw id; UniProtKB:XXXXX -> XXXXX
        else:
            logger.warning(f"get_uniprotkb_genename unsucessful for {product_id}. Product id must be supplied in the UniProtKB:XXXXX format!")
            return None
        
        gene_name = ""
        ortholog_lines = [] # lines which contain product_id, but not as the second element
        for line in self._readlines:
            if product_id in line:
                line_elements = line.split("\t")
                if line_elements[1] == product_id:
                    gene_name = line_elements[2]
                    return gene_name
                elif line_elements[1] != product_id:
                    ortholog_lines.append(line)
        
        if gene_name == "" and len(ortholog_lines) > 0:
            # goaf file was read, but product_id was never the second line element,
            # but was found in some lines to be an ortholog to some genes? or is maybe involved in some genes?
            # TODO: implement logic

            # for ortho_line in ortholog_lines:
            #   ...

            gene_name = "" # delete this
        
        if gene_name == "":
            # if gene_name is still not found, return None
            
            return None

        
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
        self.uniprot_query_exceptions = []
    
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
        return self._process_uniprot_info_query_results(results, uniprot_id)
        
    def _return_mane_select_values_from_uniprot_query(self, result: dict) -> tuple:
        """
        Given the UniProt search result dictionary, return Ensembl gene ID, Ensembl transcript ID, and RefSeq nucleotide ID for the MANE-select transcript.
        """
        mane_indices = [index for (index, d) in enumerate(result["uniProtKBCrossReferences"]) if d["database"] == "MANE-Select"]
        if len(mane_indices) == 1:
            i = mane_indices[0]
            enst_id = result["uniProtKBCrossReferences"][i]["id"]
            refseq_nt_id = next((entry["value"] for entry in result["uniProtKBCrossReferences"][i]["properties"] if entry["key"] == "RefSeqNucleotideId"), None)
            ensg_id = next((next((sub["value"] for sub in entry["properties"] if sub["key"] == "GeneId"), None) for entry in result["uniProtKBCrossReferences"] if (entry["database"] == "Ensembl" and entry["id"] == enst_id)), None)
            return ensg_id, enst_id, refseq_nt_id
        else:
            return None, None, None
    
    def _process_uniprot_info_query_results(self, results:str, uniprot_id:str) -> dict:
        """
        Processes the results obtained from the get_uniprot_info query.
        """
        if len(results) == 0:
            return {}
        else:
            # Get values from the UniProt search result
            result = next((entry for entry in results if entry["primaryAccession"] == uniprot_id), None)
            name = result["genes"][0]["geneName"]["value"]
            if "proteinDescription" in result and "recommendedName" in result["proteinDescription"] and "fullName" in result["proteinDescription"]["recommendedName"] and "value" in result["proteinDescription"]["recommendedName"]["fullName"]:
                description = result["proteinDescription"]["recommendedName"]["fullName"]["value"]
            elif "submissionNames" in result["proteinDescription"]:
                # some entries, such as UniProtKB:A0A0G2JMH6 don't have recommendedName in proteinDescription, but follow this pattern: result->proteinDescription->submissionNames->List[0: fullName -> value].
                # there can be multiple proposed descriptions, this code accounts for them all:
                description = ""
                submissionNames = result['proteinDescription']['submissionNames']
                for i in range(len(submissionNames)):
                    if i == 0:
                        description = submissionNames[i]['fullName']['value']
                    else:
                        description += f", {submissionNames[i]['fullName']['value']}"
                # resulting description is the accumulation of all comma-delimited descriptions
            else:
                description = "ERROR: Couldn't fetch description."
                logger.warning(f"proteinDescription, recommendedName, fullName or value not found when querying for uniprot info for the id: {uniprot_id}")
                logger.warning(f"result: {result}")
            ensg_id, enst_id, refseq_nt_id = self._return_mane_select_values_from_uniprot_query(result)
            return {"genename": name, "description": description, "ensg_id": ensg_id, "enst_id": enst_id, "refseq_nt_id": refseq_nt_id}
    
    async def get_uniprot_info_async(self, uniprot_id: str, session: aiohttp.ClientSession) -> dict:
        # Extract UniProt ID if given in "database:identifier" format
        logger.info(f"querying info for {uniprot_id}")

        if ":" in uniprot_id:
            uniprot_id = uniprot_id.split(":")[1]
        # Construct UniProt API query URL
        url = f"https://rest.uniprot.org/uniprotkb/search?query={uniprot_id}+AND+organism_id:9606&format=json&fields=accession,gene_names,organism_name,reviewed,xref_ensembl,xref_refseq,xref_mane-select,protein_name"
        
        # async with session.get(url) as response:
        #    results = await response.json()["results"]
        #    return self._process_uniprot_info_query_results(results, uniprot_id)

        QUERY_RETRIES = 3 # TODO: make parameter
        i = 0
        for _ in range (QUERY_RETRIES):
            if i == (QUERY_RETRIES-1):
                return None
            i += 1
            try:
                response = await session.get(url, timeout=5)
            except(requests.exceptions.RequestException, TimeoutError, asyncio.CancelledError, asyncio.exceptions.TimeoutError, aiohttp.ServerDisconnectedError) as e:
                logger.warning(f"Exception when querying info for {uniprot_id}. Exception: {str(e)}")
                self.uniprot_query_exceptions.append({f"{uniprot_id}": f"{str(e)}"})
                await asyncio.sleep(2) # sleep before retrying
                continue
            
        # single query retry
        #try:
        #    response = await session.get(url, timeout=5)
        #except (requests.exceptions.RequestException, TimeoutError, asyncio.CancelledError, asyncio.exceptions.TimeoutError) as e:
        #    logger.warning(f"Exception when querying info for {uniprot_id}. Exception: {str(e)}")
        #    self.uniprot_query_exceptions.append({f"{uniprot_id}": f"{str(e)}"})
        #    return None
        
        response_json = await response.json()
        results = response_json["results"]
        return self._process_uniprot_info_query_results(results, uniprot_id)
        
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
        self.ortholog_query_exceptions = [] # the list of exceptions during the ortholog query

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

    async def get_human_ortholog_async(self, id, session: aiohttp.ClientSession):
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
            response = await session.get(url, headers={"Content-Type": "application/json"}, timeout=10)
            # response.raise_for_status()
        except (requests.exceptions.RequestException, TimeoutError, asyncio.CancelledError, asyncio.exceptions.TimeoutError) as e:
            logger.warning(f"Exception for {id_url} for request: https://rest.ensembl.org/homology/symbol/{species}/{id_url}?target_species=human;type=orthologues;sequence=none. Exception: {str(e)}")
            self.ortholog_query_exceptions.append({f"{id}": f"{str(e)}"})
            return None

        # TODO: implement this safety check, server may send text only, which causes error (content_type == "text/plain")
        #if response.content_type == "application/json":
        #    response_json = await response.json()
        response_json = await response.json()
        logger.info(f"response_json: {response_json}; expr = {response_json==[]}")
        if response_json == [] or "error" in response_json:
            return None
        elif response_json != [] or "error" not in response_json:
            response_json = response_json["data"][0]["homologies"]
            if response_json == []: # if there are no homologies, return None
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
        if "Error" in id: # this is a bugfix. Older versions had a string "[RgdError_No-human-ortholog-found:product_id=RGD:1359312" for the genename field, if no ortholog was found (for example for the genename field of "RGD:1359312"). This is to be backwards compatible with any such data.json(s). An error can also be an '[MgiError_No-human-ortholog-found:product_id=MGI:97618'
            logger.debug(f"ERROR: {id}. This means a particular RGD, Zfin, MGI or Xenbase gene does not have a human ortholog and you are safe to ignore it.")
            return {}
        
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
    
    async def get_info_async(self, id:str, session:aiohttp.ClientSession):
        """Can receive Ensembl id or symbol (human)

        Args:
            id (str): Ensembl ID or symbol

        Returns:
            dict: Information about the gene
        """
        if "Error" in id: # this is a bugfix. Older versions had a string "[RgdError_No-human-ortholog-found:product_id=RGD:1359312" for the genename field, if no ortholog was found (for example for the genename field of "RGD:1359312"). This is to be backwards compatible with any such data.json(s). An error can also be an '[MgiError_No-human-ortholog-found:product_id=MGI:97618'
            logger.debug(f"ERROR: {id}. This means a particular RGD, Zfin, MGI or Xenbase gene does not have a human ortholog and you are safe to ignore it.")
            return {}
        
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
            response = await session.get(
                f"https://rest.ensembl.org/lookup/{endpoint}?mane=1;expand=1",
                headers={"Content-Type": "application/json"},
                timeout=5
                )
            response.raise_for_status()
            response_json = await response.json()
        except (requests.exceptions.RequestException, TimeoutError, asyncio.CancelledError, asyncio.exceptions.TimeoutError):
            # If the request fails, try the xrefs URL instead
            try:
                response = await session.get(
                    f"https://rest.ensembl.org/xrefs/{endpoint}?",
                    headers={"Content-Type": "application/json"},
                    timeout=5
                )
                response.raise_for_status()
                response_json = await response.json()
                # Use the first ENS ID in the xrefs response to make a new lookup request
                ensembl_id = next((xref["id"] for xref in response_json if "ENS" in xref["id"]), None)
                if ensembl_id:
                    response = await session.get(
                        f"https://rest.ensembl.org/lookup/id/{ensembl_id}?mane=1;expand=1",
                        headers={"Content-Type": "application/json"},
                        timeout=5
                    )
                    response.raise_for_status()
                    response_json = await response.json()
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
                response = await session.get(
                    f"https://rest.ensembl.org/xrefs/id/{ensembl_transcript_id}?all_levels=1;external_db=UniProt%",
                    headers={"Content-Type": "application/json"},
                    timeout=5
                )
                response.raise_for_status() # TODO: solve Too Many Requests error (429) -> aiohttp.client_exceptions.ClientResponseError: 429, message='Too Many Requests', url=URL('https://rest.ensembl.org/xrefs/id/ENST00000301012?all_levels=1;external_db=UniProt%25')
                response_json = await response.json()
            except (requests.exceptions.RequestException, TimeoutError, asyncio.CancelledError, asyncio.exceptions.TimeoutError):
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
    def __init__(self, zfin_filepath:str = "", xenbase_filepath:str = "", mgi_filepath:str = "", rgd_filepath:str="", goaf_filepath:str = ""):
        """
        Constructs the HumanOrthologFinder, which uses file-based search on pre-downloaded 3rd party database ortholog mappings to find
        ortholog genes.

        Parameters:
          - (str) zfin_filepath: Filepath to the Zebrafish Information Network human ortholog mapping file, found at https://zfin.org/downloads -> Orthology Data -> Human and Zebrafish Orthology -> link = https://zfin.org/downloads/human_orthos.txt
          - (str) xenbase_filepath: Filepath to the Xenbase human ortholog mapping file, found at https://www.xenbase.org/ -> Download -> Data Download (https://www.xenbase.org/xenbase/static-xenbase/ftpDatafiles.jsp) -> Data Reports -> Orthology -> Xenbase genes to Human Entrez Genes -> link: https://download.xenbase.org/xenbase/GenePageReports/XenbaseGeneHumanOrthologMapping.txt
          - (str) mgi_filepath: Filepath to the Mouse Genome Informatics human ortholog mapping file, found at: https://www.informatics.jax.org/ -> Download (https://www.informatics.jax.org/downloads/reports/index.html) -> Vertebrate homology -> Human and Mouse Homology Classes with Sequence information (tab-delimited) -> link = https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt
          - (str) rgd_filepath: Filepath to the Rat Genoma Database human ortholog mapping file, found at: https://rgd.mcw.edu/ -> Data -> Download -> data/release -> RGD_ORTHOLOGS.txt -> link = https://download.rgd.mcw.edu/data_release/RGD_ORTHOLOGS.txt
                                TODO: RGD also offers various other ortholog files, of use may be the Ensembl ortholog file, which also offers some ensembl ids: RGD_ORTHOLOGS_Ensembl.txt (https://download.rgd.mcw.edu/data_release/RGD_ORTHOLOGS_Ensembl.txt)
          
        The files are expected to reside in root/src_data_files/ folder.
        """ 
        self.zfin = ZFINHumanOrthologFinder(filepath=zfin_filepath)
        self.xenbase = XenbaseHumanOrthologFinder(filepath=xenbase_filepath)
        self.mgi = MGIHumanOrthologFinder(filepath=mgi_filepath)
        self.rgd = RGDHumanOrthologFinder(filepath=rgd_filepath)
        self.goaf = GOAnnotiationsFile()

    def find_human_ortholog(self, product):
        """
        Finds the human ortholog for the given product.

        Args:
            product (str): The product (id) for which to find the human ortholog.
        
        Returns:
            The human gene symbol or None if no human ortholog was found.
        """
        if "ZFIN" in product:
            result = self.zfin.find_human_ortholog(product) # returns [0]: gene symbol, [1]: long name of the gene
            human_gene_symbol = result[0] if result != None else None
            #human_gene_symbol = self.zfin.find_human_ortholog(product)[0]
            #return None if "Error" in human_gene_symbol else human_gene_symbol
        elif "Xenbase" in product:
            result = self.xenbase.find_human_ortholog(product)
            human_gene_symbol = result[0] if result != None else None
            #human_gene_symbol = self.xenbase.find_human_ortholog(product)[0]
            #return None if "Error" in human_gene_symbol else human_gene_symbol
        elif "MGI" in product:
            human_gene_symbol = self.mgi.find_human_ortholog(product)
            human_gene_symbol = human_gene_symbol if (human_gene_symbol != None) else None
            # return None if "Error" in human_gene_symbol else human_gene_symbol
        elif "RGD" in product:
            human_gene_symbol = self.rgd.find_human_ortholog(product)
            human_gene_symbol = human_gene_symbol if (human_gene_symbol != None) else None
            # return None if "Error" in human_gene_symbol else human_gene_symbol
        else:
            logger.info(f"No database found for {product}")
        
        return human_gene_symbol
    
    async def find_human_ortholog_async(self, product):
        # TODO: START FROM HERE. Create async file browsing for all other ortholog finders.
        if "ZFIN" in product:
            result = await self.zfin.find_human_ortholog_async(product) # returns [0]: gene symbol, [1]: long name of the gene
            return result[0] if result != None else None
            #human_gene_symbol = self.zfin.find_human_ortholog(product)[0]
            #return None if "Error" in human_gene_symbol else human_gene_symbol
        elif "Xenbase" in product:
            result = await self.xenbase.find_human_ortholog_async(product)
            return result[0] if result != None else None
            #human_gene_symbol = self.xenbase.find_human_ortholog(product)[0]
            #return None if "Error" in human_gene_symbol else human_gene_symbol
        elif "MGI" in product:
            human_gene_symbol = await self.mgi.find_human_ortholog_async(product)
            return human_gene_symbol if (human_gene_symbol != None) else None
            # return None if "Error" in human_gene_symbol else human_gene_symbol
        elif "RGD" in product:
            human_gene_symbol = await self.rgd.find_human_ortholog_async(product)
            return human_gene_symbol if (human_gene_symbol != None) else None
            # return None if "Error" in human_gene_symbol else human_gene_symbol
        else:
            logger.info(f"No database found for {product}")
            return None

class ZFINHumanOrthologFinder(HumanOrthologFinder):
    def __init__(self, filepath:str=""):
        """
        This class allows the user to search Zebrafish human orthologs. The human orthologs mapping file should be downloaded
        from the ZFIN webpage: https://zfin.org/downloads -> Orthology Data -> Human and Zebrafish Orthology -> link = https://zfin.org/downloads/human_orthos.txt

        Parameters:
          - (str) filepath: if left to default value, self._filepath will be set to "src_data_files/zfin_human_ortholog_mapping.txt", else
                            self._filepath will be set to src_data_files/{filepath}
        """
        if filepath == "":
            self._filepath = "src_data_files/zfin_human_ortholog_mapping.txt"
        else:
            if "src_data_files/" not in filepath:
                self._filepath = f"src_data_files/{filepath}"
            else:
                self._filepath = filepath
        self._check_file()
        with open(self._filepath, "r") as read_content:
            self._readlines = read_content.readlines()
        logger.info(f"ZFINHumanOrthologFinder setup ok: {len(self._readlines)} readlines.")
    
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
        return None
        # return [f"ZfinError_No-human-ortholog-found:product_id={product_id}"]
    
    async def find_human_ortholog_async(self, product_id):
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
        return None
        # return [f"ZfinError_No-human-ortholog-found:product_id={product_id}"]

class XenbaseHumanOrthologFinder(HumanOrthologFinder):
    def __init__(self, filepath:str=""):
        """
        This class allows the user to search Xenbase human orthologs. The human orthologs mapping file should be downloaded
        from the Xenbase webpage: https://www.xenbase.org/ -> Download -> Data Download (https://www.xenbase.org/xenbase/static-xenbase/ftpDatafiles.jsp) -> Data Reports -> Orthology -> Xenbase genes to Human Entrez Genes -> link: https://download.xenbase.org/xenbase/GenePageReports/XenbaseGeneHumanOrthologMapping.txt
        
        Parameters:
          - (str) filepath: if left to default value, self._filepath will be set to "src_data_files/xenbase_human_ortholog_mapping.txt", else
                            self._filepath will be set to src_data_files/{filepath}
        """
        if filepath == "":
            self._filepath = "src_data_files/xenbase_human_ortholog_mapping.txt"
        else:
            if "src_data_files/" not in filepath:
                self._filepath = f"src_data_files/{filepath}"
            else:
                self._filepath = filepath

        self._check_file()
        with open(self._filepath, "r") as read_content:
            self._readlines = read_content.readlines()
        logger.info(f"XenbaseHumanOrthologFinder setup ok: {len(self._readlines)} readlines.")
    
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
        return None
        # return [f"[XenbaseError_No-human-ortholog-found:product_id={product_id}"]
    
    async def find_human_ortholog_async(self, product_id):
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
        return None
        # return [f"[XenbaseError_No-human-ortholog-found:product_id={product_id}"]

class MGIHumanOrthologFinder(HumanOrthologFinder):
    def __init__(self, filepath:str=""):
        """
        This class allows the user to search MGI human orthologs. The human orthologs mapping file should be downloaded
        from the MGI webpage: Filepath to the Mouse Genome Informatics human ortholog mapping file, found at: 
        https://www.informatics.jax.org/ -> Download (https://www.informatics.jax.org/downloads/reports/index.html) -> Vertebrate homology -> Human and Mouse Homology Classes with Sequence information (tab-delimited) -> link = https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt
        
        Parameters:
          - (str) filepath: if left to default value, self._filepath will be set to "src_data_files/mgi_human_ortholog_mapping.txt", else
                            self._filepath will be set to src_data_files/{filepath}
        """
        if filepath == "":
            self._filepath = "src_data_files/mgi_human_ortholog_mapping.txt"
        else:
            if "src_data_files/" not in filepath:
                self._filepath = f"src_data_files/{filepath}"
            else:
                self._filepath = filepath

        self._check_file()
        with open(self._filepath, "r") as read_content:
            self._readlines = read_content.readlines()
        logger.info(f"MGIHumanOrthologFinder setup ok: {len(self._readlines)} readlines.")
    
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
                if "MgiError" in human_symbol:
                    logger.info(f"Couldn't find human ortholog for mgi gene {product_id}")
                    return None
                logger.info(f"Found human ortholog {human_symbol} for mgi gene {product_id}")
                return human_symbol # return here doesnt affect line counter 'i', since if gene is found i is no longer needed
            i += 1
        return None
        # return f"[MgiError_No-human-ortholog-found:product_id={product_id}"
    
    async def find_human_ortholog_async(self, product_id):
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
                if "MgiError" in human_symbol:
                    logger.info(f"Couldn't find human ortholog for mgi gene {product_id}")
                    return None
                logger.info(f"Found human ortholog {human_symbol} for mgi gene {product_id}")
                return human_symbol # return here doesnt affect line counter 'i', since if gene is found i is no longer needed
            i += 1
        return None
        # return f"[MgiError_No-human-ortholog-found:product_id={product_id}"

class RGDHumanOrthologFinder(HumanOrthologFinder):
    def __init__(self, filepath:str=""):
        """
        This class allows the user to search RGD human orthologs. The human orthologs mapping file should be downloaded
        from the RGD webpage: https://rgd.mcw.edu/ -> Data -> Download -> data/release -> RGD_ORTHOLOGS.txt -> link = https://download.rgd.mcw.edu/data_release/RGD_ORTHOLOGS.txt
                              TODO: RGD also offers various other ortholog files, of use may be the Ensembl ortholog file, which also offers some ensembl ids: RGD_ORTHOLOGS_Ensembl.txt (https://download.rgd.mcw.edu/data_release/RGD_ORTHOLOGS_Ensembl.txt)
        
        Parameters:
          - (str) filepath: if left to default value, self._filepath will be set to "src_data_files/rgd_human_ortholog_mapping.txt", else
                            self._filepath will be set to src_data_files/{filepath}
        """
        if filepath == "":
            self._filepath = "src_data_files/rgd_human_ortholog_mapping.txt"
        else:
            if "src_data_files/" not in filepath:
                self._filepath = f"src_data_files/{filepath}"
            else:
                self._filepath = filepath

        self._check_file()
        with open(self._filepath, "r") as read_content:
            self._readlines = read_content.readlines()
        logger.info(f"RGDHumanOrthologFinder setup ok: {len(self._readlines)} readlines.")
    

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
                logger.warning(f"FAULTY LINE IN RGD while searching for {product_id}, linesplit =: {linesplit}")
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
        return None
        # return f"[RgdError_No-human-ortholog-found:product_id={product_id}"
    
    async def find_human_ortholog_async(self, product_id):
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
                logger.warning(f"FAULTY LINE IN RGD while searching for {product_id}, linesplit =: {linesplit}")
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
        return None
        # return f"[RgdError_No-human-ortholog-found:product_id={product_id}"
