from .Model import ReverseLookup, Product, miRNA
from typing import List, Dict, Set, Optional, Callable

class Metrics:
    def __init__(self, model: ReverseLookup):
        self.reverse_lookup = model
        self.name :str = None

    def metric(self, product: Product):
        raise NotImplementedError("Subclasses must implement metric()")

class adv_product_score(Metrics):
    def __init__(self, model: ReverseLookup, a: float = 10, b1: float = 2, b2: float = 0.5, c1: float = 1, c2: float = 0.1):
        self.name = "adv_score"
        super().__init__(model) 
        self.a = a
        self.b1 = b1
        self.b2 = b2
        self.c1 = c1
        self.c2 = c2

    def metric(self, product: Product)-> float:
        goterms_list = self.reverse_lookup.get_all_goterms_for_product(product)
        score = 0.0

        def _opposite_direction(self, direction: str) -> str:
            if direction == "0":
                return "0"
            elif direction == "+":
                return "-"
            elif direction == "-":
                return "+"
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
                any(_opposite_direction(process['direction']) == goterm.direction and process['process'] == goterm.process for goterm in goterms_list)
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
                any(_opposite_direction(process['direction']) == goterm.direction and process['process'] == goterm.process for goterm in goterms_list)
                for process in self.reverse_lookup.target_processes
            )
            # Check if none of the processes in target_processes have a GOTerm in goterms_list that regulates it in the same direction
            and not any(
                any(process['direction'] == goterm.direction and process['process'] == goterm.process for goterm in goterms_list)
                for process in self.reverse_lookup.target_processes
            )
        ):
            # If all target processes are regulated in the opposite direction, subtract a points from the score
            score -= self.a

        # Calculate the score based on the number of processes in target_processes that are regulated
        # by GOTerms in the same direction as defined in the list
        score += sum(
            (self.b1**(self.b2 * sum(
                # Check if the direction and process in the process dict matches with the direction and process in any GOTerm dict
                goterm.weight for goterm in goterms_list if process['direction'] == goterm.direction and process['process'] == goterm.process)))
                for process in self.reverse_lookup.target_processes
        )
        # Calculate the score based on the number of processes in target_processes that are regulated
        # by GOTerms in the oposite direction as defined in the list
        score -= sum(
            (self.b1**(self.b2 * sum(
                # Check if the direction and process in the process dict matches with the direction and process in any GOTerm dict
                goterm.weight for goterm in goterms_list if _opposite_direction(process['direction']) == goterm.direction and process['process'] == goterm.process)))
                for process in self.reverse_lookup.target_processes
        )

        # Calculate the score by multiplying the current score with a factor based on the number of GOTerms with direction "0"
        score = score * (
            self.c1  # Start with a base factor of 1
            + (self.c2  # Add a factor based on a constant value c
                * sum(  # Multiply c by the sum of weights of all GOTerms with direction "0"
                    goterm.weight  # Get the weight of each GOTerm
                    for goterm in goterms_list  # Iterate over all GOTerms in the list
                    if goterm.direction == "0"  # Only consider GOTerms with direction "0"
                )
            )
        )

        return score

class nterms(Metrics):
    def __init__(self, model: ReverseLookup):
        super().__init__(model) 
        self.name = "nterms"
    def metric(self, product: Product) -> dict:

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

class inhibited_products_id(Metrics):
    def __init__(self, model: ReverseLookup):
        super().__init__(model)
        self.name = "inhibited_products_id"
        self.treshold = self.reverse_lookup.miRNA_overlap_treshold # it should be defined in the model, otherwise strange things happen when one mixes scores with different treshold
    def metric(self, mirna: miRNA) -> List[str]:
        inhibited_product_id = []
        for product_id, overlap in mirna.mRNA_overlaps.items():
            if overlap >= self.treshold:
                inhibited_product_id.append(product_id)
        return inhibited_product_id

class basic_mirna_score(Metrics):
    def __init__(self, model: ReverseLookup):
        super().__init__(model)
        self.name = "basic_score"
        self.treshold = self.reverse_lookup.miRNA_overlap_treshold # it should be defined in the model, otherwise strange things happen when one mixes scores with different treshold
    def metric(self, mirna: miRNA) -> float:
        score = 0.0
        for product_id, overlap in mirna.mRNA_overlaps.items():
            product = next((x for x in self.reverse_lookup.products if x.uniprot_id == product_id), None)
            if product is not None:
                if overlap >= self.treshold: #inhibited
                    a = -1 #deduct the score, since high score indicates the products is favourable for our target processes
                else:
                    a = 1
                score += a * product.scores["adv_score"]
        return score
