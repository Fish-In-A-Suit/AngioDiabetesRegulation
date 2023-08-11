"""
This test file demonstrates how to compare two ReverseLookup instances using the
(ReverseLookup).compare_to function.
"""

from goreverselookuplib import ReverseLookup
from goreverselookuplib.AnnotationProcessor import GOAnnotiationsFile
from goreverselookuplib.Workflows import WorkflowOne
from goreverselookuplib.AnnotationProcessor import GOApi
import asyncio
import aiohttp
from goreverselookuplib.JsonUtil import SimpleNamespaceUtil, JsonToClass
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s [%(levelname)s] %(funcName)s: %(message)s')
logger = logging.getLogger(__name__)

# load two models, change some data
model1 = ReverseLookup.load_model("diabetes_angio_2/data_06-06-2023.json")
model2 = ReverseLookup.load_model("diabetes_angio_2/data_06-06-2023.json")
model2.goterms[0].name = "faulty-name"
model2.goterms[3].description = "faulty-description"
model2.goterms[32].products = ["faulty-product"]
model2.products[40].scores["adv_score"] = -10000
model2.products[30].scores["nterms"] = -10000
model2.products[15].scores["binomial_test"]["diabetes-"]["pvalue_corr"] = -1000
model2.products[20].scores["fisher_test"]["diabetes-"]["pvalue_corr"] = -1000

goterms_diff = model1.compare_to(model2, "goterms", compare_subfields=['name','description','weight','products'])
products_diff = model1.compare_to(model2, "products", compare_subfields=['id_synonyms','description','uniprot_id','ensg_id','enst_id','refseq_nt_id','mRNA','scores_adv-score','scores_nterms','scores_binomial-test','scores_fisher-test'])

## logger.info(f"Entire differences json: {goterms_diff}")
## logger.info(f"Entire differences json: {products_diff}")

i = 0
for diff in goterms_diff:
    print("\n")
    i+=1
    logger.info(f"[{i}]: difference at: {diff}, difference: {goterms_diff[diff]}")
    
for diff in products_diff:
    print("\n")
    i+=1
    logger.info(f"[{i}]: difference at: {diff}, difference: {products_diff[diff]}")

# compare goterms and products both in one function call
all_diff = model1.compare_to(model2, compare_field="", compare_subfields=[])
for diff in all_diff:
    print("\n")
    i+=1
    logger.info(f"[{i}]: difference at: {diff}, difference: {all_diff[diff]}")

