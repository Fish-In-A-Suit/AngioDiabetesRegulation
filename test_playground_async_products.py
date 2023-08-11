"""
This test file demonstrates the usage of async functions to obtain GO Term products and descriptions.
It also compares the newly obtained GO Term products against an older ReverseLookup model instance, which has been computed using
synchronous requests at a previous time. It analyses the differences between the newly obtained products and the products
obtained in the old ReverseLookup model instance. Finally, for all the changes at each GO Term (async obtained products for GO Terms which
are different than the synchronously obtained products), a new synchronous call is made.

In my observation, after comparing the new synchronously obtained products against the same async products,
no changes could be found. Thus it can be seen that Gene Ontology is constantly updating the genes (products)
associated with GO Terms.
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

reference_model = ReverseLookup.load_model("diabetes_angio_2/data_06-06-2023.json")
src_model = ReverseLookup.from_input_file("diabetes_angio_4/input.txt")
api = GOApi()

model_async = ReverseLookup.from_input_file("diabetes_angio_4/input.txt")
model_async.fetch_all_go_term_names_descriptions(run_async=True)
model_async.fetch_all_go_term_products(web_download=True, run_async=True, recalculate=False, delay=0.0, run_async_options="v3", 
                                       request_params={"rows":50000}, max_connections=60)
model_async.save_model("diabetes_angio_4/model_async_test.json")

# log differences
all_diff = model_async.compare_to(reference_model, compare_field="goterms", compare_subfields=["products"])
i=0
for diff in all_diff: 
    print("\n")
    i+=1
    logger.info(f"[{i}]: difference at: {diff}, difference: {all_diff[diff]}")

# log http errors
logger.info(f"Http errors when querying products: ")
j = 0
for goterm in model_async.goterms:
    if goterm.http_error_codes != {}:
        logger.info(f"    - {goterm.http_error_codes['products']}")
        j+=1
logger.info(f"In total, there were http errors for querying products of {j} goterms (total goterms = {len(model_async.goterms)})")

# FOR ALL ELEMENTS, WHICH ARE DIFFERENT THAN REFERENCE MODEL, RUN SYNCHRONOUS GET REQUESTS AND COMPARE
api = GOApi()
previous_goterms = {} # goterms obtained using async, which hold different products than reference model
diffs_after_sync = [] # differences after sync to async request products have been compared
for diff in all_diff:
    # note: diff = goterm id; all_diff[diff] = json string representing the mismatches
    goterm_id = diff
    # find goterm with goterm_id from model_async.goterms
    for goterm in model_async.goterms:
        if goterm.id == goterm_id:
            # append a separate copy to previous goterms, since products will change after synchronous fetch
            previous_goterms[goterm_id] = goterm.copy()
            # synchronously fetch products
            goterm.fetch_products(api)
            # compare new products of this goterm to the products inside previous_goterms
            products_sync_async_diff = goterm.compare_products(previous_goterms[goterm_id]) # src is goterm with sync request, ref is goterm with async request
            sync_products_not_in_async = products_sync_async_diff[0] # products obtained using synchronous http request, which weren't obtained async
            async_products_not_in_sync = products_sync_async_diff[1] # products obtained using asynchronous http request, which weren't obtained using sync
            if sync_products_not_in_async != [] or async_products_not_in_sync != []:
                diffs_after_sync.append(f"go_id: {goterm_id}, sync_not_async products = {sync_products_not_in_async}, async_not_sync products = {async_products_not_in_sync}")

logger.info(f"Differences after sync - async comparisons:")
i = 0
for diff in diffs_after_sync:
    logger.info(f"    - [{i}]: {diff}")
    i+=1