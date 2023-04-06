import json
import requests
import time
import os
from typing import List, Dict, Set, Optional
from tqdm import tqdm

import logging
# Initialize logger
logger = logging.getLogger(__name__)


class GOTerm:
    def __init__(self, id: str, process: str, regulation_type: str, name: Optional[str] = None, description: Optional[str] = None, weight: float = 1.0, products: List[str] = []):
        """
        A class representing a Gene Ontology term.

        Args:
            id (str): The ID of the GO term.
            process (str): The name of the biological process associated with the GO term.
            regulation_type (str): The type of regulation associated with the GO term (e.g. positive or negative or general).
            name (str): Name (optional).
            description (str): A description of the GO term (optional).
            weight (float): The weight of the GO term.
            products (list): Products associated with the term (optional).
        """
        self.id = id
        self.process = process
        self.regulation_type = regulation_type
        self.name = name
        self.description = description
        self.weight = weight
        self.products = products

    def fetch_name_description(self):
        url = f"http://api.geneontology.org/api/ontology/term/{self.id}"
        params = {}
        retries = 3
        while retries > 0:
            try:
                response = requests.get(url, params=params, timeout=5)
                if response.ok:
                    data = response.json()
                    self.name = data['label']
                    self.description = data['definition']
                    logger.info(f"Fetched name and description for GO term {self.id}")
                    break
                else:
                    logger.warning(f"Error: {response.status_code} - {response.reason}")
            except requests.exceptions.Timeout:
                logger.warning("Timeout error. Retrying...")
            except requests.exceptions.RequestException as e:
                logger.warning(f"Error: {e}")
            retries -= 1
            time.sleep(2)

    @classmethod
    def from_dict(cls, d: dict):
        """
        Creates a GOTerm object from a dictionary.

        Args:
            d (dict): A dictionary containing the GO term data.

        Returns:
            A new instance of the GOTerm class.
        """
        goterm = cls(d['id'], d['process'], d['regulation_type'], d.get('name'), d.get('description'), d.get('weight', 1.0), d.get('products', []))
        return goterm

    def fetch_products(self):
        """
        Fetches product IDs associated with this GO term from the Gene Ontology API
        """
        url = f"http://api.geneontology.org/api/bioentity/function/{self.id}/genes"
        params = {"rows": 10000000}
        retries = 3
        self.products = []
        products_set = set()
        while retries > 0:
            try:
                response = requests.get(url, params=params, timeout=5)
                response.raise_for_status()
                json = response.json()
                for assoc in json['associations']:
                    if assoc['object']['id'] == self.id and assoc['subject']['taxon']['id'] == "NCBITaxon:9606":
                        product_id = assoc['subject']['id']
                        products_set.add(product_id)
                self.products = list(products_set)
                logger.info(f"Fetched products for GO term {self.id}")
                break
            except (requests.exceptions.RequestException, ValueError):
               retries-=1

class Product:
    def __init__(self, id: str, uniprot_id: str = None, description: str = None, mRNA: str = None, scores: dict = {}):
        """
        A class representing a product (e.g. a gene or protein).

        Args:
            id (str): The ID of the product.
            uniprot_id (str): The UniProt ID of the product.
            description (str): A description of the product.
            mRNA (str): The mRNA sequence of the product.
            scores (dict): A dictionary of scores associated with the product (e.g. expression score, functional score).
        """
        self.id = id
        self.uniprot_id = uniprot_id
        self.description = description
        self.mRNA = mRNA
        self.scores = dict(scores)

    @classmethod
    def from_dict(cls, d: dict) -> 'Product':
        """
        Class method to create a new Product instance from a dictionary.

        Args:
            d (dict): The dictionary containing the data to create the Product instance.

        Returns:
            Product: A new Product instance created from the input dictionary.
        """
        return cls(d['id'], d.get('uniprot_id'), d.get('description'), d.get('mRNA'), d.get('scores', {}))

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
        for goterm in tqdm(self.goterms):
            goterm.fetch_name_description()

    def fetch_all_go_term_products(self):
        """
        Iterates over all GOTerm objects in the go_term set and calls the fetch_products method for each object.
        """
        for goterm in tqdm(self.goterms):
            goterm.fetch_products()

    def load_go_term_datafile(self, filename: str) -> None:
        with open(filename, 'r') as f:
            data = json.load(f)
        for element in data['goterms']:
            for term in self.goterms:
                if term.id in element.values():
                    term.name = element.get('name')
                    term.description = element.get('description')
                    term.products = element.get('products')

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

    def create_products_from_goterms(self) -> None:
        products_set = set()
        for term in self.goterms:
            products_set.update(term.products)
        for product in products_set:
            self.products.append(Product.from_dict({'id':product}))
    
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
                    target_processes.append({"process":chunks[0], "regulation_goal":chunks[1]})
                elif section == "GO":
                    chunks = line.split(LINE_ELEMENT_DELIMITER)
                    if len(chunks) == 5:
                        d = {"id":chunks[0], "process":chunks[1], "regulation_type":chunks[2], "weight":chunks[3], "description": chunks[4]}
                    else:
                        d = {"id":chunks[0], "process":chunks[1], "regulation_type":chunks[2], "weight":chunks[3]}
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


class ReverseLookupIO:
    @staticmethod
    def load(filename: str) -> 'ReverseLookup':
        """
        Load a set of Gene Ontology terms and associated products from a file.

        Args:
            filename (str): The name of the file to load.

        Returns:
            ReverseLookup: A `ReverseLookup` object containing the GO terms and products.
        """
        with open(filename, 'r') as f:
            data = json.load(f)

        reverse_lookup = ReverseLookup()

        for goterm_data in data['goterms']:
            goterm = GOTerm.from_dict(goterm_data)
            reverse_lookup.add_goterm(goterm)

        for product_data in data['products']:
            product = Product.from_dict(product_data)
            reverse_lookup.add_product(product)

        return reverse_lookup

    @staticmethod
    def save(filename: str, reverse_lookup: 'ReverseLookup') -> bool:
        """
        Save a set of Gene Ontology terms and associated products to a file.

        Args:
            filename (str): The name of the file to save to.
            reverse_lookup (ReverseLookup): A `ReverseLookup` object containing the GO terms and products.

        Returns:
            bool: True if the file was saved successfully, False otherwise.
        """
        data = {
            'goterms': [goterm.to_dict() for goterm in reverse_lookup.goterms],
            'products': [product.to_dict() for product in reverse_lookup.products]
        }

        try:
            with open(filename, 'w') as f:
                json.dump(data, f, indent=4)
            return True
        except:
            return False


