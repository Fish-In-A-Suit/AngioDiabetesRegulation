# Testing file

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

model = ReverseLookup.load_model("diabetes_angio_4/model_async_scored.json")
t = ""