# This file showcases the entire program usage as per the latest version
# (with async queries).

from goreverselookuplib import ReverseLookup
from goreverselookuplib.AnnotationProcessor import GOAnnotiationsFile
from goreverselookuplib.Workflows import WorkflowOne, WorkflowTwo
from goreverselookuplib.AnnotationProcessor import GOApi
import asyncio
import aiohttp
from goreverselookuplib.JsonUtil import SimpleNamespaceUtil, JsonToClass
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

# Init Cacher for async queries
Cacher.init()

# Load model from input file
model = ReverseLookup.from_input_file("diabetes_angio_4/input.txt")
# or load pre-computed model: model = ReverseLookup.load_model("diabetes_angio_2/data_06-06-2023.json")

# Fetch all GO term names and descriptions
model.fetch_all_go_term_names_descriptions(run_async=True)
# Fetch all GO term products
model.fetch_all_go_term_products(web_download=True, run_async=True, recalculate=False, max_connections = 60, request_params={"rows":50000}, delay = 0.0)
# Create product instances from GO terms
model.create_products_from_goterms()
# Fetch human ortholog for products (either UniProtID, ENSG or genename)
model.fetch_ortholog_products(refetch=False, max_connections=15, req_delay=0.1, semaphore_connections=5)
model.prune_products()
# Fetch product information (from UniprotAPI or EnsemblAPI)
model.fetch_product_infos(refetch=False, run_async=True, max_connections=15, semaphore_connections=10, req_delay=0.1)
model.prune_products()
# Save model
model.save_model("diabetes_angio_4/model.json")

# Score products (only fisher_test is relevant here)
goaf = GOAnnotiationsFile()
# nterms_test = nterms(model)
# adv_product_score_test = adv_product_score(model)
fisher_test = fisher_exact_test(model, goaf)
# binom_test = binomial_test(model, goaf)
# model.score_products(score_classes=[nterms_test, adv_product_score_test, fisher_test, binom_test])
model.score_products(score_classes=[fisher_test])
model.save_model("diabetes_angio_4/model_scored.json")

# Perform statistical analysis of the scored products
model.perform_statistical_analysis(test_name="fisher_test", filepath="diabetes_angio_4/results.json")
model.save_model("diabetes_angio_4/model_analysed.json")


# Alternatively, use workflows:
use_workflows = False
if use_workflows:
    Cacher.init()
    workflow = WorkflowTwo(
        input_file_fpath="diabetes_angio_5/input.txt", 
        save_folder_dir="diabetes_angio_5",
        name="diabetes-angiogenesis"
    )
    # workflow = WorkflowTwo(
    #    input_file_fpath="diabetes_angio_5/data.json", 
    #    save_folder_dir="diabetes_angio_5",
    #    name="diabetes-angiogenesis"
    #    )
    workflow.run_workflow()