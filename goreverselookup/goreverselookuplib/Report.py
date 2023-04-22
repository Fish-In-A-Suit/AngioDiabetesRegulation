from __future__ import annotations
from typing import TYPE_CHECKING, Optional
if TYPE_CHECKING:
    from .Model import ReverseLookup
    from .Metrics import Metrics
import os
import tabulate
from tabulate import tabulate, SEPARATING_LINE
import traceback
import logging

logger = logging.getLogger(__name__)

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
    
    def _generate_top_bottom_products_summary(self, score_key) -> str:
        # Initialize the summary string with the header
        string = f"TOP and BOTTOM {self.top_n} PRODUCTS\n"

        # Get the top and bottom products based on the advanced score
        sorted_products = sorted(self.reverse_lookup.products, key=lambda x: x.scores[score_key], reverse=True)
        top_products = sorted_products[:self.top_n]
        bottom_products = sorted_products[-self.top_n:]

        # If verbosity is at least 1, create a table of the top and bottom products with their scores and descriptions
        if self.verbosity >= 1:
            # Create the table as a list of lists and append each row
            table = [["Gene Name", "Score", "Description"]]
            for product in top_products:
                table.append([product.genename, f"{product.scores[score_key]:.2f}", product.description])
            # Add a separator row and append each row for the bottom products
            table.append(["----", "----", "----"])
            for product in bottom_products:
                table.append([product.genename, f"{product.scores[score_key]:.2f}", product.description])
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
                string += " "*10+f"{product.genename} - {product.scores[score_key]:.2f} - {product.description}".center(100)+"\n"
                string += tabulate(create_go_table(product), headers="firstrow", tablefmt="grid", maxcolwidths=100) + "\n\n"

            # Add a separator and add the details for the bottom products
            string += ("-"*30)+"\n\n"
            for product in bottom_products:
                string += " "*10+f"{product.genename} - {product.scores[score_key]:.2f} - {product.description}".center(100)+"\n"
                string += tabulate(create_go_table(product), headers="firstrow", tablefmt="grid", maxcolwidths=100) + "\n\n"

        return string

    def _generate_top_miRNAs_summary(self, products_score_key, mirna_score_key) -> str:
        # Create the header string.
        string = f"TOP {self.top_n} miRNAs" + "\n"
        string += f"+ annotates Product which is in top {self.top_n}, and - annotates Product which is in bottom {self.top_n}\n\n" 

        # Get the top and bottom products based on the advanced score
        sorted_products = sorted(self.reverse_lookup.products, key=lambda x: x.scores[products_score_key], reverse=True)
        top_products = sorted_products[:self.top_n]
        bottom_products = sorted_products[-self.top_n:]

        # Get the top miRNAs.
        top_miRNAs = sorted(self.reverse_lookup.miRNAs, key=lambda x: x.scores[mirna_score_key], reverse=True)[:self.top_n]

        # If verbosity is set to 1, create a table with the top miRNAs and their scores.
        if self.verbosity == 1:
            table = [["miRNA", "mirna_score_key"]]
            for _miRNA in top_miRNAs:
                table.append([_miRNA.id, _miRNA.scores[mirna_score_key]])
            string += tabulate(table, headers="firstrow", tablefmt="grid") + "\n\n"

        # If verbosity is set to 2 or higher, create a table with the top miRNAs, their scores, and the products they inhibit.
        if self.verbosity >= 2:
            table = [["miRNA", "mirna_score_key", "suppressed products"]]

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
                table.append([_miRNA.id, _miRNA.scores[mirna_score_key], "\n".join(item for item in temp_list)])
            string += tabulate(table, headers="firstrow", tablefmt="grid") + "\n\n"

        # Return the summary string.
        return string

    def general_report(self, filepath:str, product_score: Optional[Metrics] = None, miRNA_score: Optional[Metrics] = None):
        """
        Generates the general report and writes it to a file.

        Args:
        filepath (str): The path to the output file.
        """
        filepath = filepath.replace("/", os.sep) # # replace any '/' to avoid having both \\ and / in a single filepath
        
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
            # TODO: bottom line results in error, if no product_score is supplied !!!
            report += self._generate_top_bottom_products_summary(product_score.name) + "\n"

        # Generate section on miRNAs
        if len(self.reverse_lookup.miRNAs) > 0:
            report += self._generate_section("miRNAs")
            report += self._generate_top_miRNAs_summary(miRNA_score.name) + "\n" 

        if not os.path.isabs(filepath):
            current_dir = os.path.dirname(os.path.abspath(traceback.extract_stack()[0].filename))
            mac_filepath = os.path.join(current_dir, filepath) # mac_filepath, since this approach works on a mac computer

        # Create directory for the report file, if it does not exist
        try:
            os.makedirs(os.path.dirname(mac_filepath), exist_ok=True) # this approach works on a mac computer
        
            # Write the report to the output file
            with open(mac_filepath, 'w') as f:
                f.write(report)
        except OSError:
            # TODO
            # first pass is allowed, on Windows 10 this tries to create a file at 
            # 'C:\\Program Files\\Python310\\lib\\diabetes_angio_1/general.txt'
            # which raises a permission error.
            pass  

        # fallback if the above fails
        try:
            current_dir = os.getcwd()
            win_filepath = os.path.join(current_dir, filepath) # reassign filepath to absolute path
            os.makedirs(os.path.dirname(win_filepath), exist_ok=True)

            with open(win_filepath, 'w') as f:
                f.write(report)
        except OSError:
            logger.info(f"ERROR! Cannot make directory at {os.path.dirname(filepath)}")
