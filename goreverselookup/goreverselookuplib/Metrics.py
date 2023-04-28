from __future__ import annotations
from typing import TYPE_CHECKING, Dict
if TYPE_CHECKING:
    from .Model import ReverseLookup, Product, miRNA
    from .AnnotationProcessor import GOAnnotiationsFile
from typing import List

from scipy.stats import binomtest, combine_pvalues, fisher_exact
import statistics

class Metrics:
    """
    A super-class (interface) for the scoring used in the GO Reverse Lookup. It has to be implemented by a subclass, specifically you have
    to implement the 'metric' function. It is designed to measure the regulatory effect of Gene Ontology (GO) terms on a product.
    """
    def __init__(self, model: ReverseLookup):
        self.reverse_lookup = model
        self.name :str = None

    def metric(self, product: Product | miRNA):
        """
        The 'metric' function should be implemented in the subclasses of this interface.
        """
        raise NotImplementedError("Subclasses must implement metric()")


class adv_product_score(Metrics):
    """
    An advanced scoring algorithm, an implementation of the Metrics interface. It takes in  a model (in the form of a ReverseLookup object) 
    and several parameters (a, b1, b2, c1, c2) in it's constructor, which are used to tune the weights given to different aspects of the scoring algorithm.

    Parameters:
      - (ReverseLookup) model: an instance of the ReverseLookup model.
      - (float) a: is used to give a base score to a product when all target processes are regulated in the same direction as the GOTerms in the list.
      - (float) b1, b2: are used to calculate the score based on the number of processes in target_processes that are regulated by GOTerms in the same (b1) or opposite (b2) direction as defined in the list.
      - (float) c1, c2: are used to adjust the score based on the number of GOTerms with direction "0"
    
    Scoring algorithm (explanation of the metric function):
      1. Start with score = 0.0
      
      2. (a) If all GO Terms of a Product instance regulate the processes of the ReverseLookup instance (eg. angio, diabetes) positively (and none negatively), then add 'a' to score. The 'direction' value of a positive regulatory term is '+', whereas direction for a negative regulatory term is '-'.
      
      3. (b1) For each of the processes, compute sum(goterm.weight ** b2), for every GO Term of the product, which positively regulates the process.
         Final equation ADDED to the score is sum(b1 * sum(goterm.weight ** b2)). The first 'sum' is the sum of processes, whereas the second 'sum' is the sum of GO Terms, which pass the positive regulation check for the current process.
      
      4. (b2) For each of the process, compute sum(goterm.weight ** b2), for every GO Term of the product, which negatively regulates the process.
         Final equation SUBTRACTED from the score is sum(b1 * sum(goterm.weight ** b2)). The first 'sum' is the sum of processes, the second 'sum' is the sum of GO Terms, which pass the negative regulation check for the current process.
      
      5. (c1, c2): For the GO Terms of the product with "general" regulation (direction = 0), add a base score of c1 and add sum(c2 * goterm.weight) for every GO Term with direction = 0 (general regulation).
         Final equation ADDED to the score is score * (c1 + sum(c2 * goterm.weight)), the sum relating to the number of GO Terms with direction == 0. 

    Example of calling and usage:
    1. Construct a ReverseLookup model
    model = ReverseLookup.from_input_file("diabetes_angio_1/input.txt")

    2. Create the scoring object
    adv_score = adv_product_score(model)

    3. Score the products
    model.score_products(adv_score)
    """
    def __init__(self, model: ReverseLookup, a: float = 10, b1: float = 2, b2: float = 0.5, c1: float = 1, c2: float = 0.1):
        super().__init__(model)
        self.name = "adv_score"
        self.a = a
        self.b1 = b1
        self.b2 = b2
        self.c1 = c1
        self.c2 = c2

    def metric(self, product: Product)-> float:
        """
        An implementation of the scoring algorithm for an input Product instance.

        Parameters:
          - (Product) product: an instance of a Product

        Returns:
          - (float) score: a score according to this class' scoring algorithm.
        """
        # a list of GO Terms associated with the current Product
        goterms_list = self.reverse_lookup.get_all_goterms_for_product(product)
        score = 0.0

        def _opposite_direction(direction: str) -> str:
            if direction == "0":
                return "0"
            elif direction == "+":
                return "-"
            elif direction == "-":
                return "+"
            
        # Check if all target processes are regulated in the same direction as the GOTerms in the list
        # and none of them are regulated in the opposite direction
        if (
            # Check if all processes in target_processes of the ReverseLookup model
            # have a GOTerm in goterms_list that regulates it (them) in the same direction
            all(
                any(any(process['direction'] == p['direction'] and process['process'] == p['process'] for p in goterm.processes) for goterm in goterms_list)
                for process in self.reverse_lookup.target_processes
            )
            # Check if none of the processes in target_processes have a GOTerm in goterms_list that regulates it in the opposite direction
            and not any(
                any(any(_opposite_direction(process['direction']) == p['direction'] and process['process'] == p['process'] for p in goterm.processes) for goterm in goterms_list)
                for process in self.reverse_lookup.target_processes
            )
        ):
            # If all target processes are regulated in the same direction, add a points to the score
            score += self.a

        # Check if all target processes are regulated in the opposite direction as the GOTerms in the list
        # and none of them are regulated in the same direction
        if (
            # Check if all processes in target_processes have a GOTerm in goterms_list that regulates it in the opposite direction
            all(
                any(any(_opposite_direction(process['direction']) == p['direction'] and process['process'] == p['process'] for p in goterm.processes) for goterm in goterms_list)
                for process in self.reverse_lookup.target_processes
            )
            # Check if none of the processes in target_processes have a GOTerm in goterms_list that regulates it in the same direction
            and not any(
                any(any(process['direction'] == p['direction'] and process['process'] == p['process'] for p in goterm.processes) for goterm in goterms_list)
                for process in self.reverse_lookup.target_processes
            )
        ):
            # If all target processes are regulated in the opposite direction, subtract a points from the score
            score -= self.a

        # Calculate the score based on the number of processes in target_processes that are regulated
        # by GOTerms in the same direction as defined in the list
        score += sum(
            (self.b1 * sum(
                # Check if the direction and process in the process dict matches with the direction and process in any GOTerm dict
                goterm.weight for goterm in goterms_list if any(process['direction'] == p['direction'] and process['process'] == p['process'] for p in goterm.processes)
                ) ** self.b2)
                for process in self.reverse_lookup.target_processes
        )
        # Calculate the score based on the number of processes in target_processes that are regulated
        # by GOTerms in the oposite direction as defined in the list
        score -= sum(
            (self.b1 * sum(
                # Check if the direction and process in the process dict matches with the direction and process in any GOTerm dict
                goterm.weight for goterm in goterms_list if any(_opposite_direction(process['direction']) == p['direction'] and process['process'] == p['process'] for p in goterm.processes)
                ) ** self.b2)
                for process in self.reverse_lookup.target_processes
        )

        # Calculate the score by multiplying the current score with a factor based on the number of GOTerms with direction "0"
        score = score * (
            self.c1  # Start with a base factor of 1
            + (self.c2  # Add a factor based on a constant value c
                * sum(  # Multiply c by the sum of weights of all GOTerms with direction "0"
                    goterm.weight  # Get the weight of each GOTerm
                    for goterm in goterms_list  # Iterate over all GOTerms in the list
                    if any(p['direction'] == 0 for p in goterm.processes)  # Only consider GOTerms with direction "0"
                )
            )
        )

        return score

class nterms(Metrics):
    """
    An implementation of the Metrics interface, it scores the products by positive, negative or general regulation of a speciffic process.
    It takes in a model (in the form of a ReverseLookup object) in it's constructor.

    Parameters:
      - (ReverseLookup) model: an instance of the ReverseLookup model.
    
    Scoring algorithm (explanation of the metric function):
      Create an empty nterms_dict, where descriptive regulatory keys (eg. 'angio+', 'angio-', 'angio0') will be mapped to a count of terms regulating a specific process in a specific direction
      For each process in the model's (ReverseLookup) 'target_processes':
        a) create the following keys in nterms_dict: '{process}+', '{process}-', '{process}0'; if process is 'angio', then 'angio+', 'angio-', 'angio0' will be the keys in nterms_dict
        b) populate each of the keys with the count of GO Terms, which positively (direction == '+'), negatively (direction == '-') or generally (direction == '0') regulate the process
    
    Example of calling and usage:
    1. Create a ReverseLookup model
    model = ReverseLookup.from_input_file("diabetes_angio_1/input.txt")

    2. Create the scoring object
    nterms_score = nterms(model)

    3. Use the scoring object using model.score_products
    model.score_products(nterms_score)
    """
    def __init__(self, model: ReverseLookup):
        super().__init__(model) 
        self.name = "nterms"
    def metric(self, product: Product) -> dict:
        """
        An implementation of the scoring algorithm for an input Product instance.

        Parameters:
          - (Product) product: an instance of a Product

        Returns:
          - (dict) nterms_dict: a dictionary with (ReverseLookup).target_processes * 3 keys. Each process of a ReverseLookup instance has 3 keys,
                                '{process}+', '{process}-', '{process}0'. Each key has an integer count value of the amount of GO Terms of the input Product instance,
                                which positively (direction == '+'), negatively (direction == '-') or generally (direction == '0') regulate a speciffic process.

                                For a ReverseLookup model with defined processed 'angio' and 'diabetes', the returned dictionary would have 6 keys:
                                angio+, angio-, angio0, diabetes+, diabetes-, diabetes0
        """
        # A list of GO Terms associated with the current Product
        goterms_list = self.reverse_lookup.get_all_goterms_for_product(product)
        # An empty dictionary to store the count of GOTerms for each process and direction
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

class binomial_test(Metrics):
    def __init__(self, model: ReverseLookup, goaf: GOAnnotiationsFile):
        super().__init__(model) 
        self.goaf = goaf
        self.name = "binomial_test"
        self._num_all_goterms = 0
    
    def metric(self, product: Product) -> Dict:
        
        if self._num_all_goterms == 0:
            self._num_all_goterms = len(self.goaf.get_all_terms())
        
        results_dict = {}
        
        for process in self.reverse_lookup.target_processes:
            process_goterms_list = self.reverse_lookup.get_all_goterms_for_process(process["process"])
            num_goterms_product_general = len(self.goaf.get_all_terms_for_product(product.genename))
            num_goterms_all_general = self._num_all_goterms
            for direction in ['+', '-']:
                num_goterms_product_process = sum(1 for goterm in process_goterms_list if (any(i['direction'] == direction for i in goterm.processes) and (any(product_id in goterm.products for product_id in product.id_synonyms) or product.genename in goterm.products)))
                num_goterms_all_process = sum(1 for goterm in process_goterms_list if any(i['direction'] == direction for i in goterm.processes))
                
                #time for Binomial test and "risk ratio"
                #binom = binomtest(num_goterms_product_process, num_goterms_all_process, 
                #                  (num_goterms_product_general/num_goterms_all_general), alternative='greater')
                binom = binomtest(num_goterms_product_process, num_goterms_all_process, 
                                  (num_goterms_product_general/num_goterms_all_general if num_goterms_all_general != 0 else 0), alternative='greater') # bugfix: ZeroDivisionError
                binom_pvalue = binom.pvalue
                
                if num_goterms_product_general != 0 and num_goterms_all_general != 0: # bugfix: ZeroDivisionError
                    risk_ratio = (num_goterms_product_process/num_goterms_all_process if num_goterms_all_process != 0 else 0) / (num_goterms_product_general/num_goterms_all_general)
                else:
                    risk_ratio = 0

                fold_enrichment_score = 0
                if num_goterms_all_process != 0 and num_goterms_product_general != 0 and num_goterms_all_general != 0:
                    fold_enrichment_score = num_goterms_product_process / (num_goterms_all_process * (num_goterms_product_general / num_goterms_all_general))
                    
                results_dict[f"{process['process']}{direction}"] = {
                    #"n_prod_process" : num_goterms_product_process,
                    #"n_all_process" : num_goterms_all_process,
                    #"n_prod_general" : num_goterms_product_general,
                    #"n_all_general" : num_goterms_all_general,
                    "num" : num_goterms_product_process,
                    "expected" : num_goterms_all_process * (num_goterms_product_general / num_goterms_all_general if num_goterms_all_general != 0 else 0),
                    "fold_enrichment" : fold_enrichment_score, # bugfix: ZeroDivisionError
                    "pvalue" : binom_pvalue,
                    "risk_ratio" : risk_ratio,
                }
        
        #all_target_pvalues = [results_dict[f"{process['process']}{process['direction']}"]['pvalue'] for process in self.reverse_lookup.target_processes]
        #combined_p = combine_pvalues(all_target_pvalues)
        #results_dict["comb_binom_pvalue"] = combined_p.pvalue
        #combined_rr = statistics.mean([results_dict[f"{process['process']}{process['direction']}"]['risk_ratio'] for process in self.reverse_lookup.target_processes])
        #results_dict["comb_risk_ratio"] = combined_rr
        
        return results_dict
            
class fisher_exact_test(Metrics):
    def __init__(self, model: ReverseLookup, goaf: GOAnnotiationsFile):
        super().__init__(model) 
        self.goaf = goaf
        self.name = "fisher_test"
        self._num_all_goterms = 0
    
    def metric(self, product: Product) -> Dict:
        
        if self._num_all_goterms == 0:
            self._num_all_goterms = len(self.goaf.get_all_terms())
        
        results_dict = {}
        
        for process in self.reverse_lookup.target_processes:
            process_goterms_list = self.reverse_lookup.get_all_goterms_for_process(process["process"])
            num_goterms_product_general = len(self.goaf.get_all_terms_for_product(product.genename))
            num_goterms_all_general = self._num_all_goterms
            for direction in ['+', '-']:
                num_goterms_product_process = sum(1 for goterm in process_goterms_list if (any(i['direction'] == direction for i in goterm.processes) and (any(product_id in goterm.products for product_id in product.id_synonyms) or product.genename in goterm.products)))
                num_goterms_all_process = sum(1 for goterm in process_goterms_list if any(i['direction'] == direction for i in goterm.processes))
                
                #time for Binomial test and "risk ratio"
                cont_table = [[num_goterms_product_process, num_goterms_all_process-num_goterms_product_process],
                              [num_goterms_product_general-num_goterms_product_process, num_goterms_all_general-(num_goterms_all_process-num_goterms_product_process)]]
                fisher = fisher_exact(cont_table)

                fisher_pvalue = fisher.pvalue
                
                odds_ratio = fisher.statistic
                
                folder_enrichment_score = 0
                if num_goterms_all_process != 0 and num_goterms_product_general != 0 and num_goterms_all_general != 0:
                    folder_enrichment_score = num_goterms_product_process / (num_goterms_all_process * (num_goterms_product_general / num_goterms_all_general))

                results_dict[f"{process['process']}{direction}"] = {
                    "n_prod_process" : num_goterms_product_process,
                    "n_all_process" : num_goterms_all_process,
                    "n_prod_general" : num_goterms_product_general,
                    "n_all_general" : num_goterms_all_general,
                    "num" : num_goterms_product_process,
                    "expected" : num_goterms_all_process * (num_goterms_product_general / num_goterms_all_general if num_goterms_all_general != 0 else 0),
                    "fold_enrichment" : folder_enrichment_score, # BUGFIX: ZeroDivisionError
                    "pvalue" : fisher_pvalue,
                    "odds_ratio" : odds_ratio,
                }
        
        #all_target_pvalues = [results_dict[f"{process['process']}{process['direction']}"]['pvalue'] for process in self.reverse_lookup.target_processes]
        #combined_p = combine_pvalues(all_target_pvalues)
        #results_dict["comb_fisher_pvalue"] = combined_p.pvalue
        #combined_rr = statistics.mean([results_dict[f"{process['process']}{process['direction']}"]['odds_ratio'] for process in self.reverse_lookup.target_processes])
        #results_dict["comb_odds_ratio"] = combined_rr
        
        return results_dict

class inhibited_products_id(Metrics):
    """
    An implementation of the Metrics interface to return a list of all product ids inhibited by a specific miRNA, if the binding strength
    between a product id and a specific miRNA is greater than (ReverseLookup).miRNA_overlap_threshold.
    
    WARNING: field 'miRNA_overlap_threshold' must be defined in an instance of the ReverseLookup model passed to this constructor.
    
    Parameters:
      - (ReverseLookup) model: an instance of ReverseLookup

    Algorithm:
      Create an empty list.
      For the input miRNA, loop over all miRNA-mRNA binding strengths (stored in (miRNA).mRNA_overlaps)
        If binding strength > miRNA_overlap_threshold: append product id to list
      Return a list of product ids
    """
    def __init__(self, model: ReverseLookup):
        super().__init__(model)
        self.name = "inhibited_products_id"
        self.treshold = self.reverse_lookup.miRNA_overlap_treshold # it should be defined in the model, otherwise strange things happen when one mixes scores with different treshold
    
    def metric(self, mirna: miRNA) -> List[str]:
        """
        An implementation of the scoring algorithm for a specific miRNA instance. It loops over all miRNA-mRNA binding strengths in (miRNA).mRNA_overlaps
        and returns a list of mRNA product ids, whose binding strengths to this miRNA are greater than miRNA_overlap_threshold.
        """
        inhibited_product_ids = []
        for product_id, overlap in mirna.mRNA_overlaps.items():
            if overlap >= self.treshold:
                inhibited_product_ids.append(product_id)
        return inhibited_product_ids

class basic_mirna_score(Metrics):
    """
    Score calculated from adv_score and overlap

    Scoring algorithm:
        Initialise score = 0.0
        For each product_id and it's float overlap (binding strength) value in (miRNA).mRNA_overlaps
        [TODO]: explain why the miRNA score is decreased (a = -1) if it binds to a product with a good threshold? 
        if miRNA binds to a product well, then it's score should be increased!
    
    WARNING: [TODO] need to resolve scoring issue
    """
    def __init__(self, model: ReverseLookup):
        super().__init__(model)
        self.name = "basic_score"
        self.treshold = self.reverse_lookup.miRNA_overlap_treshold # it should be defined in the model, otherwise strange things happen when one mixes scores with different treshold
    
    def metric(self, mirna: miRNA) -> float:
        """
        An implementation of the scoring algorithm for a specific miRNA instance. [TODO] explain more after the scoring issue is solved
        """
        score = 0.0
        for product_id, overlap in mirna.mRNA_overlaps.items():
            product = next((x for x in self.reverse_lookup.products if x.uniprot_id == product_id), None) #  this line of code is looking through a sequence of products and finding the first product whose uniprot_id matches the value of product_id. If such a product is found, it is assigned to the variable product; otherwise, product is set to None
            if product is not None: # each miRNA can have many products in it's 'mRNA_overlaps' field, this is a check that we are only analysing the products, which are also present in (ReverseLookup).products
                if overlap >= self.treshold: # inhibited
                    a = -1 # deduct the score, since high score indicates the products is favourable for our target processes
                else:
                    a = 1
                score += a * product.scores["adv_score"]
        return score
