# Testing file

from goreverselookuplib import ReverseLookup
from goreverselookuplib.GOTerm import GOTerm
from goreverselookuplib.AnnotationProcessor import GOAnnotiationsFile
from goreverselookuplib.Workflows import WorkflowOne, WorkflowTwo
# from goreverselookuplib.AnnotationProcessor import GOApi
import asyncio
import aiohttp
from goreverselookuplib.JsonUtil import SimpleNamespaceUtil, JsonToClass, JsonUtil
from goreverselookuplib.AnnotationProcessor import HumanOrthologFinder, UniProtAPI, EnsemblAPI, GOAnnotiationsFile
from goreverselookuplib.CacheUtils import ConnectionCacher, Cacher
from goreverselookuplib.FileUtil import FileUtil
from goreverselookuplib.Metrics import fisher_exact_test, adv_product_score, nterms, binomial_test
from goreverselookuplib.OboParser import OboParser

import logging
import requests

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
#goaf = GOAnnotiationsFile(go_categories=model.go_categories)
#fisher_test = fisher_exact_test(model, goaf)
#model.score_products([fisher_test])
#model.save_model("chronic_infl_cancer_1/data_v1.json")









