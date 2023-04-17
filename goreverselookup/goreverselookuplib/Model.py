from __future__ import annotations
from .AnnotationProcessor import GOApi, EnsemblAPI, UniProtAPI, HumanOrthologFinder
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .Metrics import Metrics
import json
import os
from typing import List, Dict, Optional
from tqdm import trange, tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
import logging
import traceback

logger = logging.getLogger(__name__)


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
        goterm = cls(d['id'], d['process'], d['direction'], d.get('name'), d.get(
            'description'), d.get('weight', 1.0), d.get('products', []))
        return goterm


class Product:
    def __init__(self, id_synonyms: List[str], uniprot_id: str = None, description: str = None, ensembl_id: str = None, refseq_nt_id: str = None, mRNA: str = None, scores: dict = None):
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
            human_ortholog_gene_id = human_ortolog_finder.find_human_ortholog(
                self.id_synonyms[0])
            if human_ortholog_gene_id is not None:
                uniprot_id = uniprot_api.get_uniprot_id(human_ortholog_gene_id)
                if uniprot_id is not None:
                    self.uniprot_id = uniprot_id

    def fetch_Uniprot_info(self, uniprot_api: UniProtAPI) -> None:
        """
        includes description, ensembl_id and refseq_nt_id
        """
        if self.uniprot_id == None:
            return
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
        if self.ensembl_id is None:
            return
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

from .miRNAprediction import miRDB60predictor

class ReverseLookup:
    def __init__(self, goterms: List[GOTerm], target_processes: List[Dict[str, str]], products: List[Product] = [], miRNAs: List[miRNA] = [], miRNA_overlap_treshold: float = 0.6):
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
            self.products.append(Product.from_dict({'id_synonyms': [product]}))
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
                        product.fetch_UniprotID(
                            human_ortolog_finder, uniprot_api)
        except Exception as e:
            # If there was an exception while fetching UniProt data, save all the Product objects to a JSON file.
            self.save_products_to_datafile('crash_products.json')
            # Re-raise the exception so that the caller of the method can handle it.
            raise e

    def prune_products(self) -> None:
        # Create a dictionary that maps UniProt ID to a list of products
        reverse_uniprotid_products = {}
        for product in self.products:
            if product.uniprot_id is not None:
                reverse_uniprotid_products.setdefault(
                    product.uniprot_id, []).append(product)

        # For each UniProt ID that has more than one product associated with it, create a new product with all the synonyms
        # and remove the individual products from the list
        for uniprot_id, product_list in reverse_uniprotid_products.items():
            if len(product_list) > 1:
                id_synonyms = []
                for product in product_list:
                    self.products.remove(product)
                    id_synonyms.extend(product.id_synonyms)
                # Create a new product with the collected information and add it to the product list
                self.products.append(Product(
                    id_synonyms, uniprot_id, product_list[0].description, product_list[0].ensembl_id, product_list[0].refseq_nt_id, product_list[0].mRNA, {}))

    def fetch_Uniprot_infos(self) -> None:
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

    def score_products(self, score_class: List[Metrics]) -> None:
        if not isinstance(score_class, list):
            score_class = [score_class]
        # redirect the tqdm logging output to the logging module to avoid interfering with the normal output
        with logging_redirect_tqdm():
            # iterate over each Product object in self.products and score them using the Scoring object
            for product in tqdm(self.products):
                for _score_class in score_class:
                    product.scores[_score_class.name] = _score_class.metric(
                        product)

    def fetch_mRNA_sequences(self) -> None:
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

    def predict_miRNAs(self, prediction_type: str = 'miRDB') -> None:
        # check the prediction type
        if prediction_type == 'miRDB':
            # use the miRDB60predictor to predict miRNAs #TODO make it so that the user submitts the predictior, like metrices
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
                                self.miRNAs.append(miRNA(miRNA_id, mRNA_overlaps={
                                                   product.uniprot_id: match}))
        elif prediction_type == 'other_type':
            # do something else
            pass
        else:
            # raise an error if the prediction type is invalid
            raise ValueError("Invalid prediction type")

    def change_miRNA_overlap_treshold(self, treshold: float, yes: bool = False) -> None:
        self.miRNA_overlap_treshold = treshold
        logger.warning(
            f"Sorry, but changing the treshold will delete all the calculated miRNA scores. You will have to calculate them again!")
        if not yes:
            confirmation = input(f"Are you sure you want to proceed? (y/n) ")
            if confirmation.lower() != 'y':
                print("Aborting operation.")
                return
        for _miRNA in self.miRNAs:
            _miRNA.scores = {}

    def score_miRNAs(self, score_class: List[Metrics]) -> None:
        if not isinstance(score_class, list):
            score_class = [score_class]

        with logging_redirect_tqdm():
            # iterate over miRNAs using tqdm for progress tracking
            for mirna in tqdm(self.miRNAs):
                # if there is no overlap, skip the miRNA
                if not mirna.mRNA_overlaps:
                    continue
                for _score_class in score_class:
                    mirna.scores[_score_class.name] = _score_class.metric(
                        mirna)

# housekeeping functions

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

    def save_model(self, filepath: str) -> None:
        data = {}
        # save options - currently not used
        current_dir = os.path.dirname(os.path.abspath(
            traceback.extract_stack()[0].filename))
        # save target_process
        data['target_processes'] = self.target_processes
        data['miRNA_overlap_treshold'] = self.miRNA_overlap_treshold
        # save goterms
        for goterm in self.goterms:
            data.setdefault('goterms', []).append(goterm.__dict__)
        # save products
        for product in self.products:
            data.setdefault('products', []).append(product.__dict__)
        # save miRNAs
        for miRNA in self.miRNAs:
            data.setdefault('miRNAs', []).append(miRNA.__dict__)
        # write to file
        # Create directory for the report file, if it does not exist
        os.makedirs(os.path.dirname(os.path.join(
            current_dir, filepath)), exist_ok=True)

        with open(os.path.join(current_dir, filepath), 'w') as f:
            json.dump(data, f, indent=4)

    @classmethod
    def load_model(cls, filepath: str) -> 'ReverseLookup':
        current_dir = os.path.dirname(os.path.abspath(
            traceback.extract_stack()[0].filename))
        print(current_dir)
        with open(os.path.join(current_dir, filepath), "r") as f:
            data = json.load(f)
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

        return cls(goterms, target_processes, products, miRNAs, miRNA_overlap_treshold)

    @classmethod
    def from_input_file(cls, filepath: str) -> 'ReverseLookup':
        """
        Creates a ReverseLookup object from a JSON file.

        Args:
            filepath (str): The path to the input file.

        Returns:
            ReverseLookup: A ReverseLookup object.
        """
        # Define constants used in parsing the file
        LINE_ELEMENT_DELIMITER = '\t'  # Data is tab separated
        COMMENT_DELIMITER = "#"  # Character used to denote a comment
        LOGIC_LINE_DELIMITER = "###"  # Special set of characters to denote a "logic line"

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

        current_dir = os.path.dirname(os.path.abspath(
            traceback.extract_stack()[0].filename))
        with open(os.path.join(current_dir, filepath), "r") as read_content:
            target_processes = []
            go_terms = []
            read_lines = read_content.read().splitlines()[
                2:]  # skip first 2 lines
            section = ""  # what is the current section i am reading
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
                    target_processes.append(
                        {"process": chunks[0], "direction": chunks[1]})
                elif section == "GO":
                    chunks = line.split(LINE_ELEMENT_DELIMITER)
                    if len(chunks) == 5:
                        d = {"id": chunks[0], "process": chunks[1], "direction": chunks[2],
                             "weight": chunks[3], "description": chunks[4]}
                    else:
                        d = {"id": chunks[0], "process": chunks[1],
                             "direction": chunks[2], "weight": chunks[3]}
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
