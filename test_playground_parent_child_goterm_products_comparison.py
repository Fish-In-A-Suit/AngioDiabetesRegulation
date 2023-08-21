# This file compares the gene products associated with parent GO Terms and
# their children. The two parent terms selected are GO:0002698 (negative regulation of immune effector process)
# and GO:0002699 (positive regulation of immune effector process). We browsed the Gene Ontology
# terms tree and took all the terms "under" the aforementioned parent terms into the list of child terms.

# TLDR: Parent GO terms DO NOT store all of their children's gene products. 

from goreverselookuplib.GOTerm import GOTerm
from goreverselookuplib.Model import ReverseLookup
import logging

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s [%(levelname)s] %(funcName)s: %(message)s')
# Create a file handler
file_handler = logging.FileHandler('./log_output/test_json_dump.log')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(logging.Formatter('%(asctime)s [%(filename)s:%(lineno)s] - [%(funcName)s] %(message)s'))
logger = logging.getLogger(__name__)

parent_go_terms_list = ["GO:0002698", "GO:0002699"]
children_go_terms_list = [
    "GO:0032833",
    "GO:0032834",
    "GO:1905450",
    "GO:1905451",
    "GO:2001189",
    "GO:2001190",
    "GO:0045623",
    "GO:0045624",
    "GO:0045916",
    "GO:0045917",
    "GO:2001192",
    "GO:2001193",
    "GO:0002632",
    "GO:0002633",
    "GO:0043301",
    "GO:0043302",
    "GO:0002704",
    "GO:0002705",
    "GO:0033007",
    "GO:0033008",
    "GO:0043381",
    "GO:0043382",
    "GO:0032827",
    "GO:0032828",
    "GO:0032821",
    "GO:0032822",
    "GO:0002701",
    "GO:0002702"
]

parent_go_terms = []
child_go_terms = []

for go_id in parent_go_terms_list:
    goterm = GOTerm(id=go_id, processes=[])
    parent_go_terms.append(goterm)

for go_id in children_go_terms_list:
    goterm = GOTerm(id=go_id, processes=[])
    child_go_terms.append(goterm)

model = ReverseLookup(goterms=(parent_go_terms+child_go_terms), target_processes=[])
model.fetch_all_go_term_names_descriptions()
model.fetch_all_go_term_products(max_connections=50)
model.save_model("debug_files/parent-child-products-test.json")

parent_genes_list = []
child_genes_list = []
for goterm in model.goterms:
    is_parent = False
    if goterm.id in parent_go_terms_list:
        # goterm is parent
        parent_genes_list += goterm.products
    else:
        # goterm is child
        child_genes_list += goterm.products

# find out if there are more genes in child_genes_list than in parent_genes_list
if len(child_genes_list) > len(parent_genes_list):
    logger.info(f"There are more genes in child GO Terms than in parent GO Terms.")
else:
    logger.info(f"There are more genes in parent GO Terms than in child GO Terms.")
logger.info(f"  - child GO Terms genes count: {len(child_genes_list)}")
logger.info(f"  - parent GO Terms genes count: {len(parent_genes_list)}")

diff = []
for child_gene in child_genes_list:
    if child_gene not in parent_genes_list:
        diff.append(child_gene)
logger.info(f"Printing genes of child GO Terms not present in parent GO Terms: {diff}")
