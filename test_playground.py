# Testing file

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

# *** test fetch info with non-ensembl values ***
# Cacher.init()
#model = ReverseLookup.load_model("diabetes_angio_4/model_async_scored.json")
#products_shortened = model.products[0:10]
#for product in products_shortened:
#    product.enst_id = ""
#    product.refseq_nt_id = ""
#    product.had_fetch_info_computed = ""
#    product.had_orthologs_computed = True
#    product.fetch_info()
#t = ""

#model = ReverseLookup.load_model("chronic_infl_cancer_1/data.json")
#goaf = GOAnnotiationsFile()
#fisher_test = fisher_exact_test(model, goaf)
#model.score_products([fisher_test])
#model.save_model("chronic_infl_cancer_1/data_v1.json")


def sorting_key(target_processes, product):
    """
    Sorting key used for the sorting of JSON data based on ascending pvalues.
    Fisher test MUST be calculated for this to work.
    """
    pvalue_sum = 0
    for process in target_processes:
        pvalue_process = product["scores"]["fisher_test"][f"{process['process']}{process['direction']}"].get('pvalue_corr',1)
        pvalue_sum += pvalue_process
    return pvalue_sum

model = ReverseLookup.load_model("chronic_infl_cancer_1/data_02-09-2023.json")
statistically_relevant_products = JsonUtil.load_json("chronic_infl_cancer_1/statistically_relevant_genes_02-09-2023.json")
statistically_relevant_products = statistically_relevant_products["chronic_inflammation+:cancer+"]

before_sorting = []
for element in statistically_relevant_products:
    before_sorting.append(element["genename"])

statistically_relevant_products_sorted = sorted(statistically_relevant_products, key=lambda product: sorting_key(model.target_processes, product))

after_sorting = []
for element in statistically_relevant_products_sorted:
    after_sorting.append(element["genename"])

logger.info(f"Sorted genes according to ascending pvalue.")
logger.info(f"Before sorting: {before_sorting}")
logger.info(f"After sorting: {after_sorting}")

JsonUtil.save_json(data_dictionary=statistically_relevant_products_sorted, filepath="chronic_infl_cancer_1/stat_analysis_sorted.json")







