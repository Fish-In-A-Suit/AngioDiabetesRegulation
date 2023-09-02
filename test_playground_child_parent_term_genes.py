# According to "Use and misuse of the gene ontology annotations" (https://www.nature.com/articles/nrg2363):
#	
#	"""
# 	Because GO annotations to a term inherit all the properties of the ancestors of those terms, every path from any term back to its root(s)
# 	must be biologically accurate or the ontology must be revised. For example, if a gene is known to be specifically
# 	involved in ‘vesicle fusion’, it will be annotated directly to that term, and it is implicitly annotated (indirectly)
# 	to all of its parents’ terms, including ‘membrane fusion’, ‘membrane organization and biogenesis’, ‘vesicle-
# 	mediated transport’, ‘transport’ and so on, back to the root node (FIG. 1c). Thus, a gene annotated to vesicle fusion can
# 	be retrieved not only with this term, but also with all of its parent terms, increasing flexibility and power when
# 	searching for and making inferences about genes.
#	"""
#
# However, the current API query requests contradict the claim that "a gene annotated to vesicle fusion can
# be retrieved not only with this term, but also with all of its parent terms". For example, on 02.09.2023 a query
# was made to retrieve the genes of parent term GO:0050678 (regulation of epithelial cell proliferation) and of it's
# child term GO:0001936 (regulation of endothelial cell proliferation). It was found that the child term contains
# 24 extra genes that are not found inside the parent term. Therefore, parent term does not contain all of the genes found
# in the child terms.
#
# Parent term genes (80): ['RGD:69342', 'UniProtKB:P24530', 'MGI:104642', 'UniProtKB:Q5TA89', 'RGD:1359498', 'MGI:1298218', 'MGI:109579', 'RGD:1308353', 'MGI:102780', 'RGD:3187', 'RGD:3032', 'MGI:105382', 'RGD:1595855', 'UniProtKB:Q15475', 'UniProtKB:P55957', 'UniProtKB:Q14469', 'RGD:1308521', 'MGI:104853', 'MGI:108093', 'UniProtKB:P78504', 'RGD:1304749', 'MGI:1097716', 'RGD:1307761', 'RGD:2937', 'MGI:96522', 'MGI:95522', 'RGD:2536', 'UniProtKB:P05114', 'MGI:105923', 'UniProtKB:P78527', 'RGD:621023', 'UniProtKB:Q96QS3', 'MGI:106034', 'MGI:101884', 'RGD:3317', 'UniProtKB:P84022', 'MGI:98961', 'MGI:1341800', 'RGD:620713', 'RGD:2535', 'ZFIN:ZDB-GENE-021115-2', 'RGD:1306726', 'UniProtKB:Q8WU20', 'MGI:1100860', 'MGI:102720', 'MGI:97567', 'UniProtKB:O15230', 'RGD:2611', 'MGI:104876', 'MGI:1095416', 'RGD:1560732', 'UniProtKB:P06401', 'RGD:3851', 'MGI:97348', 'RGD:62081', 'RGD:1308982', 'RGD:1593096', 'UniProtKB:Q96T88', 'RGD:620160', 'UniProtKB:Q9H082', 'UniProtKB:Q8TAU0', 'RGD:620857', 'MGI:96120', 'MGI:1861606', 'MGI:1338889', 'UniProtKB:Q9UIU6', 'RGD:69079', 'MGI:97363', 'UniProtKB:P48059', 'MGI:104779', 'RGD:620906', 'RGD:621403', 'RGD:1308201', 'UniProtKB:P09758', 'RGD:621340', 'RGD:3370', 'UniProtKB:P36952', 'MGI:95523', 'MGI:1201674', 'RGD:1562672']
# Child term genes (26): ['RGD:621642', 'RGD:71052', 'RGD:620857', 'MGI:104663', 'MGI:104642', 'MGI:99906', 'UniProtKB:P37023', 'ZFIN:ZDB-MIRNAG-041217-16', 'RGD:628896', 'MGI:1338938', 'UniProtKB:O76093', 'MGI:96683', 'UniProtKB:Q10981', 'RGD:2639', 'UniProtKB:Q15389', 'MGI:109375', 'UniProtKB:P41159', 'UniProtKB:P35590', 'RGD:2608', 'RGD:70989', 'RGD:3000', 'UniProtKB:P19526', 'RGD:2638', 'MGI:88180', 'UniProtKB:O43320', 'MGI:98664']
# Child term extra genes (24): ['RGD:621642', 'RGD:71052', 'MGI:104663', 'MGI:99906', 'UniProtKB:P37023', 'ZFIN:ZDB-MIRNAG-041217-16', 'RGD:628896', 'MGI:1338938', 'UniProtKB:O76093', 'MGI:96683', 'UniProtKB:Q10981', 'RGD:2639', 'UniProtKB:Q15389', 'MGI:109375', 'UniProtKB:P41159', 'UniProtKB:P35590', 'RGD:2608', 'RGD:70989', 'RGD:3000', 'UniProtKB:P19526', 'RGD:2638', 'MGI:88180', 'UniProtKB:O43320', 'MGI:98664']
#

from goreverselookuplib import ReverseLookup
from goreverselookuplib.GOTerm import GOTerm
from goreverselookuplib.AnnotationProcessor import GOAnnotiationsFile
from goreverselookuplib.Workflows import WorkflowOne, WorkflowTwo
from goreverselookuplib.AnnotationProcessor import GOApi
import asyncio
import aiohttp
from goreverselookuplib.JsonUtil import SimpleNamespaceUtil, JsonToClass, JsonUtil
from goreverselookuplib.AnnotationProcessor import HumanOrthologFinder, UniProtAPI, EnsemblAPI, GOAnnotiationsFile
from goreverselookuplib.CacheUtils import ConnectionCacher, Cacher
from goreverselookuplib.Metrics import fisher_exact_test, adv_product_score, nterms, binomial_test

import logging

from math import log10

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s [%(levelname)s] %(funcName)s: %(message)s')
# Create a file handler
file_handler = logging.FileHandler('./log_output/test_json_dump.log')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(logging.Formatter('%(asctime)s [%(filename)s:%(lineno)s] - [%(funcName)s] %(message)s'))
logger = logging.getLogger(__name__)


# Algorithm
term1 = GOTerm(id="GO:0050678", description="regulation of epithelial cell proliferation", processes="")
term2 = GOTerm(id="GO:0001936", description="regulation of endothelial cell proliferation", processes="")
api = GOApi()
term1.fetch_products(api)
term2.fetch_products(api)
parent_term_contains_all_child_genes = True
child_term_extra_genes = []
for gene in term2.products:
	if gene not in term1.products: # term 1 is parent
		parent_term_contains_all_child_genes = False
		child_term_extra_genes.append(gene)
		
logger.info(f"Parent gene contains all child term genes: {parent_term_contains_all_child_genes}")
logger.info(f"  - parent term genes ({len(term1.products)}): {term1.products}")
logger.info(f"  - child term genes ({len(term2.products)}): {term2.products}")
logger.info(f"  ---")
logger.info(f"  - child term extra genes ({len(child_term_extra_genes)}): {child_term_extra_genes}")
