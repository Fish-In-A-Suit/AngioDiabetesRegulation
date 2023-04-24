from __future__ import annotations
from .AnnotationProcessor import GOApi, GOAnnotiationsFile, EnsemblAPI, UniProtAPI, HumanOrthologFinder
from typing import TYPE_CHECKING, Set, List, Dict, Optional
if TYPE_CHECKING:
    from .Metrics import Metrics
import json
import os
from tqdm import trange, tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
import logging
import traceback
from .FileUtil import FileUtil
from .Timer import Timer

logger = logging.getLogger(__name__)


class GOTerm:
    def __init__(self, id: str, processes: List[Dict], name: Optional[str] = None, description: Optional[str] = None, weight: float = 1.0, products: List[str] = []):
        """
        A class representing a Gene Ontology term.

        Args:
            id (str): The ID of the GO term.
            process (List[Dict]): [{"process" : "angio", "direction" : "+"}]
            name (str): Name (optional).
            description (str): A description of the GO term (optional).
            weight (float): The weight of the GO term.
            products (list): Products associated with the term (optional).
        """
        self.id = id
        self.processes = processes if isinstance(processes, list) else [processes]
        self.name = name
        self.description = description
        self.weight = float(weight)
        self.products = products

    def fetch_name_description(self, api: GOApi):
        """
        Sets the "name" and "description" member field of the GO Term. The link used to query for the response is http://api.geneontology.org/api/ontology/term/{term_id}.
        This function sets the "name" field of the GO Term to response['label'] and the "description" field to response['definition']

        Parameters:
          - (GOApi) api: a GOApi instance
        
        Usage and calling:
            api = GOApi()
            goterms = ["GO:1903589", ...]
            for goterm in goterms:
                goterm.fetch_name_description(api)
        """
        logger.info("Fetching GO Term names (labels) and descriptions (definitions).")
        data = api.get_data(self.id)
        if data:
            self.name = data['label']
            self.description = data['definition']
            logger.info(f"Fetched name and description for GO term {self.id}")

    def fetch_products(self, source):
        """
        Fetches UniProtKB products associated with a GO Term and sets the "products" member field of the GO Term to a list of all associated products.
        The link used to query for the response is http://api.geneontology.org/api/bioentity/function/{term_id}/genes.

        The product IDs can be of any of the following databases: UniProt, ZFIN, Xenbase, MGI, RGD 
        [TODO: enable the user to specify databases himself]

        Parameters:
          - (GOApi) api: a GOApi instance
        
        Usage and calling:
            api = GOApi() or goaf = GOAnnotationFile()
            goterms = ["GO:1903589", ...]
            for goterm in goterms:
                goterm.fetch_products(api)
        """
        products = source.get_products(self.id)
        if products:
            self.products = products

    def add_process(self, process: Dict):
        if not process in self.processes:
            self.processes.append(process)

    @classmethod
    def from_dict(cls, d: dict):
        """
        Creates a GOTerm object from a dictionary.

        Args:
            d (dict): A dictionary containing the GO term data.

        Returns:
            A new instance of the GOTerm class.
        """
        goterm = cls(d['id'], d['processes'], d.get('name'), d.get(
            'description'), d.get('weight', 1.0), d.get('products', []))
        return goterm

class Product:
    def __init__(self, id_synonyms: List[str], genename: str = None, uniprot_id: str = None, description: str = None, ensg_id: str = None, enst_id: str = None, refseq_nt_id: str = None, mRNA: str = None, scores: dict = None):
        """
        A class representing a product (e.g. a gene or protein).

        Args:
            id_synonyms (str): The list of ID of the product and synonyms. -> after ortholog translation it turns out that some products are the same. 
            uniprot_id (str): The UniProt ID of the product.
            description (str): A description of the product.
            ensg_id (str): Ensembl gene ID (MAIN).
            enst_id (str): Ensembl transcript ID.
            refseq_nt_id (str): Refseq (reference sequence) transcript ID.
            mRNA (str): The mRNA sequence of the product.
            scores (dict): A dictionary of scores associated with the product (e.g. expression score, functional score).
        """
        self.id_synonyms = id_synonyms
        self.genename = genename
        self.description = description
        self.uniprot_id = uniprot_id
        self.ensg_id = ensg_id
        self.enst_id = enst_id
        self.refseq_nt_id = refseq_nt_id
        self.mRNA = mRNA
        self.scores = {} if scores is None else scores.copy()

    def fetch_ortholog(self, human_ortolog_finder: Optional[HumanOrthologFinder] = None, uniprot_api: Optional[UniProtAPI] = None, ensembl_api: Optional[EnsemblAPI] = None) -> None:
        if not human_ortolog_finder:
            human_ortolog_finder = HumanOrthologFinder()
        if not uniprot_api:
            uniprot_api = UniProtAPI()
        if not ensembl_api:
            ensembl_api = EnsemblAPI()
        if len(self.id_synonyms) == 1 and 'UniProtKB' in self.id_synonyms[0]:
            info_dict = uniprot_api.get_uniprot_info(self.uniprot_id)
            self.genename = info_dict.get("genename")
        elif len(self.id_synonyms) == 1:
            human_ortholog_gene_id = human_ortolog_finder.find_human_ortholog(
                self.id_synonyms[0])
            if human_ortholog_gene_id is None:
                logger.warning(f"file-based human ortholog finder did not find ortholog for {self.id_synonyms[0]}")
                human_ortholog_gene_ensg_id = ensembl_api.get_human_ortholog(self.id_synonyms[0]) # attempt ensembl search
                if human_ortholog_gene_ensg_id is not None:
                    enst_dict = ensembl_api.get_info(human_ortholog_gene_ensg_id)
                    self.genename = enst_dict.get("genename")
                else:
                    return # search was unsuccessful
            else:
                self.genename = human_ortholog_gene_id
                
    def fetch_info(self, uniprot_api: Optional[UniProtAPI] = None, ensembl_api: Optional[EnsemblAPI] = None) -> None:
        """
        includes description, ensg_id, enst_id and refseq_nt_id
        """
        if not (self.uniprot_id or self.genename or self.ensg_id):
            return
        if not uniprot_api:
            uniprot_api = UniProtAPI()
        if not ensembl_api:
            ensembl_api = EnsemblAPI()

        required_keys = ["genename", "description", "ensg_id", "enst_id", "refseq_nt_id"]
        #Is uniprot really necessary. If it is faster, perhaps get uniprotID from genename and then first try to get info from uniprot
        if any(getattr(self, key) is None for key in required_keys) and self.uniprot_id:
            info_dict = uniprot_api.get_uniprot_info(self.uniprot_id)
            for key, value in info_dict.items():
                if value is not None:
                    setattr(self, key, value)
        if any(getattr(self, key) is None for key in required_keys) and self.ensg_id:
            enst_dict = ensembl_api.get_info(self.ensg_id)
            for key, value in enst_dict.items():
                if value is not None:
                    setattr(self, key, value)
        if any(getattr(self, key) is None for key in required_keys) and self.genename:
            enst_dict = ensembl_api.get_info(self.genename)
            for key, value in enst_dict.items():
                if value is not None:
                    setattr(self, key, value)
        if any(getattr(self, key) is None for key in required_keys) and self.uniprot_id:
            enst_dict = ensembl_api.get_info(self.uniprot_id)
            for key, value in enst_dict.items():
                if value is not None:
                    setattr(self, key, value)
        #TODO: logger output which values are still missing
        
    @classmethod
    def from_dict(cls, d: dict) -> 'Product':
        """
        Class method to create a new Product instance from a dictionary.

        Args:
            d (dict): The dictionary containing the data to create the Product instance.

        Returns:
            Product: A new Product instance created from the input dictionary.
        """
        return cls(d.get('id_synonyms'), d.get('genename'), d.get('uniprot_id'), d.get('description'), d.get('ensg_id'), d.get('enst_id'), d.get('refseq_nt_id'), d.get('mRNA'), d.get('scores'))

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

class TargetProcess:
    def __init__(self, name: str, direction: str) -> None:
        """
        A class representing a target process. NOT USED CURRENTLY

        Args:
            name (str)
            direction (str): + or -
            goterms (set): a set of all goterms which are 
        """
        self.name = name
        self.direction = direction
    
from .miRNAprediction import miRDB60predictor

class ReverseLookup:
    def __init__(self, goterms: List[GOTerm], target_processes: List[Dict[str, str]], products: List[Product] = [], miRNAs: List[miRNA] = [], miRNA_overlap_treshold: float = 0.6, execution_times: dict = {}):
        """
        A class representing a reverse lookup for gene products and their associated Gene Ontology terms.

        Args:
            goterms (set): A set of GOTerm objects.
            target_processes (list): A list of dictionaries containing process names and directions.
            products (set, optional): A set of Product objects. Defaults to an empty set.
        """
        self.goterms = goterms #TODO make a set, it is faster
        self.products = products #TODO make a set, it is faster
        self.target_processes = target_processes
        self.miRNAs = miRNAs #TODO make a set, it is faster
        self.miRNA_overlap_treshold = miRNA_overlap_treshold

        self.execution_times = execution_times # dict of execution times, logs of runtime for functions
        self.timer = Timer()

    def fetch_all_go_term_names_descriptions(self):
        """
        Iterates over all GOTerm objects in the go_term set and calls the fetch_name_description method for each object.
        """
        logger.info(f"Fetching GO term names and their ")
        self.timer.set_start_time()

        api = GOApi()     
        with logging_redirect_tqdm():
            for goterm in tqdm(self.goterms):
                 if goterm.name == None or goterm.description == None: # if goterm.name or description don't exist, then attempt fetch
                    goterm.fetch_name_description(api)
        
        if "fetch_all_go_term_names_descriptions" not in self.execution_times: # to prevent overwriting on additional runs of the same model name
            self.execution_times["fetch_all_go_term_names_descriptions"] = self.timer.get_elapsed_time()
        self.timer.print_elapsed_time()

    def fetch_all_go_term_products(self, recalculate: bool = False):
        """
        Iterates over all GOTerm objects in the go_term set and calls the fetch_products method for each object.
        
        Args:
          - (bool) recalculate: if set to True, will recalculate (fetch again) the term's products even if they already exist (perhaps from a model loaded from data.json)
        """
        logger.info(f"Started fetching all GO Term products.")
        self.timer.set_start_time()

        goaf = GOAnnotiationsFile()
        with logging_redirect_tqdm():
            for goterm in tqdm(self.goterms):
                if goterm.products == [] or recalculate == True: # to prevent recalculation of products if they are already computed
                    goterm.fetch_products(goaf)

        # api = GOApi()
        # with logging_redirect_tqdm():
        #     for goterm in tqdm(self.goterms):
        #         goterm.fetch_products(api)

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
        logger.info(f"Creating products from GO Terms")
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
            if ':' in product:
                self.products.append(Product.from_dict({'id_synonyms': [product]}))
            else:
                self.products.append(Product.from_dict({'id_synonyms': [product], 'genename': product}))
        logger.info(f"Created Product objects from GOTerm object definitions")

        if "create_products_from_goterms" not in self.execution_times:
            self.execution_times["create_products_from_goterms"] = self.timer.get_elapsed_time()
        self.timer.print_elapsed_time()

    def fetch_ortholog_products(self) -> None:
        """
        This function tries to find the orthologs to any non-uniprot genes (products) associated with a GO Term.

        When fetching products (genes / gene products) from Gene Ontology for a specific GO Term:
            (GOTerm).fetch_products()

            api = GOApi()
            goterms = ["GO:1903589", ...]
            for goterm in goterms:
                goterm.fetch_products(api)
        
        The resulting products can be from any of the following databases: UniProtKB, ZFIN, Xenbase, MGI, RGD. For subsequent 
        Homo-Sapiens-only product analysis, it is important to find, if human ortholog genes exist for the products, fetched from a non-uniprot
        databases.

        Usage and calling:
            products = ... # define a list of Product instances
            human_ortolog_finder = HumanOrthologFinder()
            uniprot_api = UniProtAPI()
            ensembl_api = EnsemblAPI()

            for product in products:
                product.fetch_ortholog(human_ortholog_finder, uniprot_api, ensembl_api)
        
        """
        logger.info(f"Started fetching ortholog products.")
        self.timer.set_start_time()

        try:
            human_ortholog_finder = HumanOrthologFinder()
            uniprot_api = UniProtAPI()
            ensembl_api = EnsemblAPI()
            # Iterate over each Product object in the ReverseLookup object.
            with logging_redirect_tqdm():
                for product in tqdm(self.products):
                    # Check if the Product object doesn't have a UniProt ID or genename or ensg_id.
                    if product.genename == None:
                        # If it doesn't, fetch UniProt data for the Product object.
                        product.fetch_ortholog(human_ortholog_finder, uniprot_api, ensembl_api)
        except Exception as e:
            # If there was an exception while fetching UniProt data, save all the Product objects to a JSON file.
            self.save_model('crash_products.json')
            # Re-raise the exception so that the caller of the method can handle it.
            raise e
        
        if "fetch_ortholog_products" not in self.execution_times:
            self.execution_times["fetch_ortholog_products"] = self.timer.get_elapsed_time()
        self.timer.print_elapsed_time()

    def prune_products(self) -> None:
        logger.info(f"Started pruning products.")
        self.timer.set_start_time()

        # Create a dictionary that maps genename to a list of products
        reverse_genename_products = {}
        for product in self.products:
            if product.genename is not None:
                reverse_genename_products.setdefault(
                    product.genename, []).append(product)

        # For each ENSG that has more than one product associated with it, create a new product with all the synonyms
        # and remove the individual products from the list
        for genename, product_list in reverse_genename_products.items():
            if len(product_list) > 1:
                id_synonyms = []
                for product in product_list:
                    self.products.remove(product)
                    id_synonyms.extend(product.id_synonyms)
                # Create a new product with the collected information and add it to the product list
                self.products.append(Product(
                    id_synonyms, product_list[0].genename, product_list[0].uniprot_id, product_list[0].description, product_list[0].ensg_id, product_list[0].enst_id, product_list[0].refseq_nt_id, product_list[0].mRNA, {}))

        if "prune_products" not in self.execution_times:
            self.execution_times["prune_products"] = self.timer.get_elapsed_time()
        self.timer.print_elapsed_time()

    def fetch_product_infos(self) -> None:
        # TODO: ensembl support batch request

        logger.info(f"Started fetching product infos.")
        self.timer.set_start_time()

        try:
            uniprot_api = UniProtAPI()
            ensembl_api = EnsemblAPI()
            # Iterate over each Product object in the ReverseLookup object.
            with logging_redirect_tqdm():
                for product in tqdm(self.products):
                    # Check if the Product object doesn't have a UniProt ID.
                    if any(attr is None for attr in [product.genename, product.description, product.enst_id, product.ensg_id, product.refseq_nt_id]) and (product.uniprot_id or product.genename or product.ensg_id):
                        # If it doesn't, fetch UniProt data for the Product object.
                        product.fetch_info(uniprot_api, ensembl_api)
        except Exception as e:
            raise e
        
        if "fetch_product_infos" not in self.execution_times:
            self.execution_times["fetch_product_infos"] = self.timer.get_elapsed_time()
        self.timer.print_elapsed_time()

    def score_products(self, score_classes: List[Metrics]) -> None:
        """
        Scores the products of the current ReverseLookup model. This function allows you to pass a custom or a pre-defined scoring algorithm,
        which is of 'Metrics' type (look in Metrics.py), or a list of scoring algorithms. Each Product class of the current ReverseLookup instance products (self.products)
        has a member field 'scores'. For each product, score is computed and saved to the product's 'scores' dictionary as a mapping between the 
        scoring algorithm's name (eg. "adv_score") and the corresponding product's score computed with this scoring algorithm (eg. 14.6). 
        If multiple scoring algorithms are used, then the product's 'scores' dictionary will have multiple elements, each a mapping between
        the scoring algorithm's name and the corresponding score.

        Parameters:
          - score_classes: A subclass (implementation) of the Metrics superclass (interface). Current pre-defined Metrics implementations subclasses
                         are 'adv_product_score', 'nterms', 'inhibited_products_id', 'basic_mirna_score'.
        
        Calling example:
        (1) Construct a ReverseLookup model
        model = ReverseLookup.from_input_file("diabetes_angio_1/input.txt")

        (2) Create one or more Metrics scoring implementations for the model:
        adv_score = adv_product_score(model)
        nterms_score = nterms(model)

        (3) Call the score_products on the model using the Metrics scoring implementations
        model.score_products([adv_score, nterms_score])
        """
        logger.info(f"Started scoring products.")
        self.timer.set_start_time()

        if not isinstance(score_classes, list):
            score_classes = [score_classes]
        # redirect the tqdm logging output to the logging module to avoid interfering with the normal output
        with logging_redirect_tqdm():
            # iterate over each Product object in self.products and score them using the Scoring object
            for product in tqdm(self.products): # each Product has a field scores - a dictionary between a name of the scoring algorithm and it's corresponding score
                for _score_class in score_classes:
                    pass
                    product.scores[_score_class.name] = _score_class.metric(product) # create a dictionary between the scoring algorithm name and it's score for current product
            
        for _score_class in score_classes:
            i = 0
            p_values = []
            if _score_class.name == "fisher_test" or _score_class.name == "binomial_test":   

                for product in self.products:
                    for process in self.target_processes:
                        for direction in ['+', '-']:
                            p_values.append(product.scores[_score_class.name][f"{process['process']}{direction}"]["pvalue"])
                # apply Benjamini-Hochberg FDR correction
                from statsmodels.stats.multitest import multipletests
                reject, p_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
                for product in self.products:
                    for process in self.target_processes:
                        for direction in ['+', '-']:
                            product.scores[_score_class.name][f"{process['process']}{direction}"]["pvalue_corr"] = p_corrected[i]
                            i += 1
        
        if "score_products" not in self.execution_times:
            self.execution_times["score_products"] = self.timer.get_elapsed_time()
        self.timer.print_elapsed_time()
                            
    def fetch_mRNA_sequences(self) -> None:
        logger.info(f"Started fetching mRNA sequences.")
        self.timer.set_start_time()

        try:
            ensembl_api = EnsemblAPI()
            # Iterate over each Product object in the ReverseLookup object.
            with logging_redirect_tqdm():
                for product in tqdm(self.products):
                    # Check if the Product object doesn't have a EnsemblID
                    if product.mRNA == None and product.enst_id is not None:
                        # If it has, fetch mRNA sequence data for the Product object.
                        product.fetch_mRNA_sequence(ensembl_api)
        except Exception as e:
            raise e

        if "fetch_mRNA_sequences" not in self.execution_times:
            self.execution_times["fetch_mRNA_sequences"] = self.timer.get_elapsed_time()
        self.timer.print_elapsed_time()

    def predict_miRNAs(self, prediction_type: str = 'miRDB') -> None:
        logger.info(f"Started miRNA prediction analysis.")
        self.timer.set_start_time()
        
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

        if "predict_miRNAs" not in self.execution_times:
            self.execution_times["predict_miRNAs"] = self.timer.get_elapsed_time()
        self.timer.print_elapsed_time()

    def change_miRNA_overlap_treshold(self, treshold: float, safety: bool = False) -> None:
        """
        Sets the model's 'miRNA_overlap_threshold' to a new 'threshold'. The threshold should be between 0.0 and 1.0.

        WARNING: Changing the miRNA overlap threshold will delete all the calculated previous miRNA scores.

        Parameters:
          - (float) threshold: the new miRNA_overlap_threshold
          - (bool) safety: if False, will ask for the user's confirmation during runtime.
        """
        self.miRNA_overlap_treshold = treshold
        logger.warning(f"Sorry, but changing the treshold will delete all the calculated miRNA scores. You will have to calculate them again!")
        if not safety:
            confirmation = input(f"Are you sure you want to proceed? (y/n)")
            if confirmation.lower() != 'y':
                print("Aborting operation.")
                return
        for _miRNA in self.miRNAs:
            _miRNA.scores = {}

    def score_miRNAs(self, score_class: List[Metrics]) -> None:
        """
        Performs miRNA scoring on the current ReverseLookup's 'miRNAs' using the input Metrics implementation(s). This function allows the user
        to pass a custom or a pre-defined scoring algorithm, which is of the 'Metrics' type (look in Metrics.py), or a list of scoring algorithms.
        Each miRNA class of the current ReverseLookup instance has a member field 'scores'. For each miRNA instance, score is computed
        and saved to the miRNA's 'scores' dictionary as a mapping between the scoring algorithm's name (eg. "basic_miRNA_score") and the
        corresponding miRNA's float score computed with this scoring algorithm. If multiple scoring algorithms are used, then the miRNA's 
        'scores' dictionary will have multiple elements, each a mapping between the scoring algorithm's name and the corresponding score.

        Parameters:
          - score_class: A subclass (implementation) of the Metrics superclass (interface). Current pre-defined Metrics implementations subclasses
                         are 'adv_product_score', 'nterms', 'inhibited_products_id', 'basic_mirna_score'. 

                         If 'inhibited_products_id' are used, then the miRNA's 'scoring' field will have a key "inhibited products id", the
                         value at this key will be a list of all the product ids (of all the current GOTerm-associated products, which satisfy
                         the condition that the product's mRNA binding strength > miRNA_overlap_threshold)

                         If 'basic_mirna_score' is used, then [TODO]
        
        Calling example:
        (1) Construct a ReverseLookup model
        model = ReverseLookup.from_input_file("diabetes_angio_1/input.txt")

        (2) Create one or more Metrics scoring implementations for the model:
        adv_score = adv_product_score(model)
        nterms_score = nterms(model)

        (3) Call the score_products on the model using the Metrics scoring implementations
        model.score_products([adv_score, nterms_score])
        """
        logger.info(f"Started miRNA scoring.")
        self.timer.set_start_time()

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
        
        if "score_miRNAs" not in self.execution_times:
            self.execution_times["score_miRNAs"] = self.timer.get_elapsed_time()
        self.timer.print_elapsed_time()

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
    
    def get_all_goterms_for_process(self, process: str) -> List[GOTerm]:
        goterms_list = []
        for goterm in self.goterms:
            if any(proc["process"] == process for proc in goterm.processes):
                goterms_list.append(goterm)
        return goterms_list

    def list_goterms_id(self) -> List[str]:
        """
        Returns a list of all GO term IDs in the GO ontology.
        """
        # Use a list comprehension to extract IDs from the GO terms and return the resulting list
        return [goterm.id for goterm in self.goterms]

    def get_goterm(self, identifier) -> GOTerm:
        """
        Return GOTerm based on any id
        """
        goterm = next(obj for obj in self.goterms if any(getattr(obj, attr) == identifier for attr in ["id", "name", "description"]))
        return goterm
    def get_product(self, identifier) -> Product:
        """
        Return GOTerm based on any id
        """
        product = next(obj for obj in self.products if any(getattr(obj, attr) == identifier for attr in ["genename", "description", "uniprot_id", "ensg_id", "enst_id", "refseq_nt_id", "mRNA"]) or any(id == identifier for id in obj.id_synonyms))
        return product

    def save_model(self, filepath: str) -> None:
        data = {}
        # TODO: save options - currently not used

        data['target_processes'] = self.target_processes
        data['miRNA_overlap_treshold'] = self.miRNA_overlap_treshold
        data['execution_times'] = self.execution_times

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

        try: # this works on mac, not on windows
            current_dir = os.path.dirname(os.path.abspath(traceback.extract_stack()[0].filename))
            os.makedirs(os.path.dirname(os.path.join(current_dir, filepath)), exist_ok=True) # Create directory for the report file, if it does not exist
            with open(os.path.join(current_dir, filepath), 'w') as f:
                json.dump(data, f, indent=4)
        except OSError:
            # pass the error on the first attempt
            pass

        try: # if first attempt fails, try using current_dir = os.getcwd(), this works on windows
            windows_filepath = FileUtil.find_win_abs_filepath(filepath)
            with open(windows_filepath, 'w') as f:
                json.dump(data, f, indent=4)
            #current_dir = os.getcwd()
            #os.makedirs(os.path.dirname(os.path.join(current_dir, filepath)), exist_ok=True)
        except OSError:
            logger.info(f"ERROR creating filepath {filepath} at {os.getcwd()}")

    @classmethod
    def load_model(cls, filepath: str) -> 'ReverseLookup':
        if not os.path.isabs(filepath):
            fileutil = FileUtil()
            filepath = fileutil.find_file(filepath) # attempt backtrace file search
            # current_dir = os.path.dirname(os.path.abspath(traceback.extract_stack()[0].filename))
            # filepath = os.path.join(current_dir, filepath)
        with open(filepath, "r") as f:
            data = json.load(f)
        target_processes = data['target_processes']
        miRNA_overlap_treshold = data['miRNA_overlap_treshold']

        execution_times = {}
        if "execution_times" in data:
            execution_times = data['execution_times']

        goterms = []
        for goterm_dict in data['goterms']:
            goterms.append(GOTerm.from_dict(goterm_dict))

        products = []
        for product_dict in data.get('products', []):
            products.append(Product.from_dict(product_dict))

        miRNAs = []
        for miRNAs_dict in data.get('miRNAs', []):
            miRNAs.append(miRNA.from_dict(miRNAs_dict))

        return cls(goterms, target_processes, products, miRNAs, miRNA_overlap_treshold, execution_times=execution_times)

    @classmethod
    def from_input_file(cls, filepath: str) -> 'ReverseLookup':
        """
        Creates a ReverseLookup object from a text file.

        Args:
            filepath (str): The path to the input text file.

        Returns:
            ReverseLookup: A ReverseLookup object.
        """
        # Define constants used in parsing the file
        LINE_ELEMENT_DELIMITER = '\t'  # Data is tab separated
        COMMENT_DELIMITER = "#"  # Character used to denote a comment
        LOGIC_LINE_DELIMITER = "###"  # Special set of characters to denote a "logic line"

        target_processes = []
        go_terms = []

        def process_comment(line):
            """
            Processes a comment in the line: returns the part of the line before the comment. The input file should be structured to contain
            three sections - 'settings', 'processes' and 'GO_terms', annotated using the LOGIC_LINE_DELIMITER.

            For the construction of input.txt, please refer to the Readme file. [TODO]

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
        
        def process_file(filepath: str):
            with open(filepath, "r") as read_content:
                read_lines = read_content.read().splitlines()[2:]  # skip first 2 lines
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
                        target_processes.append({"process": chunks[0], "direction": chunks[1]})
                    elif section == "GO":
                        chunks = line.split(LINE_ELEMENT_DELIMITER)
                        if len(chunks) == 5:
                            d = {"id": chunks[0], "processes":{"process": chunks[1], "direction": chunks[2]},"weight": chunks[3], "description": chunks[4]}
                        else:
                            d = {"id": chunks[0], "processes": {"process": chunks[1], "direction": chunks[2]}, "weight": chunks[3]}
                        if not any(d["id"] == goterm.id for goterm in go_terms): # TODO: check this !!!!!
                            go_terms.append(GOTerm.from_dict(d))
                        else: # TODO: check this !!!!!
                            next(goterm for goterm in go_terms if d["id"] == goterm.id).add_process({"process": chunks[1], "direction": chunks[2]})


        if not os.path.isabs(filepath): # this process with traceback.extract_stack works correctly on mac, but not on windows.
            current_dir = os.path.dirname(os.path.abspath(
                traceback.extract_stack()[0].filename))
            mac_filepath = os.path.join(current_dir, filepath) 
        
        try:
            os.makedirs(os.path.dirname(mac_filepath), exist_ok=True) # this approach works on a mac computer
            process_file(mac_filepath)
        except OSError:
            # # first pass is allowed, on Windows 10 this tries to create a file at 
            # 'C:\\Program Files\\Python310\\lib\\diabetes_angio_1/general.txt'
            # which raises a permission error.
            pass
        
        # fallback if the above fails
        try:
            win_filepath = FileUtil.find_win_abs_filepath(filepath)
            os.makedirs(os.path.dirname(win_filepath), exist_ok=True)
            process_file(win_filepath)
        except OSError:
            logger.error(f"ERROR while opening win filepath {win_filepath}")
        
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


                        


        
