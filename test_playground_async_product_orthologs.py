# This test file shows how to query product orthologs for a ReverseLookup model. You can also see the url response and ensembl/uniprot api function result caching by running in normal Run mode and either force stopping using CTRL-C or waiting for the program to finish.
# Make sure to run this file in normal Run (not in debug mode)
# Use CTRL-C if you wish to stop the program manually to test the atexit save of uniprot and ensembl cache.

from goreverselookuplib.Model import ReverseLookup
from goreverselookuplib.CacheUtils import Cacher

import logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s [%(levelname)s] %(funcName)s: %(message)s')
logger = logging.getLogger(__name__)

Cacher.init()
model_async = ReverseLookup.load_model("diabetes_angio_4/model_async_test.json") # or use a pre-computed async model

# if products for GO Terms have not been computed and saved yet, compute them:
# model_async.fetch_all_go_term_products(web_download=True, run_async=True, recalculate=False, delay=0.0, run_async_options="v3", request_params={"rows":50000}, max_connections=60)
# model_async.create_products_from_goterms()
# model_async.save_model("diabetes_angio_4/model_async_test.json")

# ortholog product fetch
model_async.fetch_ortholog_products(refetch=True, run_async=True, max_connections=15, req_delay=0.1, semaphore_connections=5) # semaphore_connections=10 works in 3min40s, semaphore_connections=15 results in 429:TooManyRequests errors