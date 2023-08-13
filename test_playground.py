from goreverselookuplib import ReverseLookup
from goreverselookuplib.AnnotationProcessor import GOAnnotiationsFile
from goreverselookuplib.Workflows import WorkflowOne
from goreverselookuplib.AnnotationProcessor import GOApi
import asyncio
import aiohttp
from goreverselookuplib.JsonUtil import SimpleNamespaceUtil, JsonToClass
from goreverselookuplib.AnnotationProcessor import HumanOrthologFinder, UniProtAPI, EnsemblAPI, GOAnnotiationsFile

import logging

from math import log10

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s [%(levelname)s] %(funcName)s: %(message)s')

# Create a file handler
file_handler = logging.FileHandler('./log_output/test_json_dump.log')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(logging.Formatter('%(asctime)s [%(filename)s:%(lineno)s] - [%(funcName)s] %(message)s'))

logger = logging.getLogger(__name__)

# before 06.06.2023
# model = ReverseLookup.load_model("diabetes_angio_2/data_web_full_scores.json")
# 06.06.2023: new code to perform statistical analysis on the model
# model = ReverseLookup.load_model("diabetes_angio_2/data_06-06-2023.json")
# model.perform_statistical_analysis(filepath="diabetes_angio_2/stat_analysis_06-06-2023.json")
# exit()

#test execution speed
#import timeit
#goaf = GOAnnotiationsFile()
#elapsed_time = timeit.timeit(lambda: goaf.get_products("GO:1903589"), number=100)
#print(elapsed_time)

#goaf = GOAnnotiationsFile()
#print(goaf.get_all_terms_for_product("ACAT2"))
#print(len(goaf.get_all_terms()))


""" # Ladi's code to perform statistical analysis for the products
if False:

    out_string = ""

    for p in model.products:
        entry = p.genename
        for pr in model.target_processes:
            for d in ['+','-']:
                entry = entry + "\t" + f"{-log10(p.scores['fisher_test'][pr['process']+d]['pvalue_corr']):.2e}"
        out_string = out_string + entry + "\n"

    with open('out.txt', 'w') as f:
        f.write(out_string)


    exit()

if True:
    model.pre

    for product in model.products:
        #processes = model.target_processes
        processes = [     
            {
                "process": "diabetes",
                "direction": "+"
            },
            {
                "process": "angio",
                "direction": "+"
            },
            {
                "process": "angio",
                "direction": "+"
            },
            #{
            #    "process": "obesity",
            #    "direction": "+"
            #}
            ]
        
    #    if (all(float(product.scores["fisher_test"][f"{process['process']}{process['direction']}"]["pvalue_corr"]) < 0.05 for process in processes)):
    #        print(f"Found significant product: {product.genename}")

        if "fisher_test" in product.scores:
            if (all(float(product.scores["fisher_test"][f"{process['process']}{process['direction']}"].get("pvalue_corr", 1)) < 0.05 for process in processes) and
                all(float(product.scores["fisher_test"][f"{process['process']}{'+' if process['direction'] == '-' else '-'}"].get("pvalue_corr", 1)) >= 0.05 for process in processes)):
                print(f"Found significant product: {product.genename}")
    exit()
model.prune_products()
"""

reference_model = ReverseLookup.load_model("diabetes_angio_2/data_06-06-2023.json")
src_model = ReverseLookup.from_input_file("diabetes_angio_4/input.txt")
api = GOApi()

### Async ortholog query testing
# load the existing model
# model_async = ReverseLookup.from_input_file("diabetes_angio_4/input.txt") # use this to work from the ground up
model_async = ReverseLookup.load_model("diabetes_angio_4/model_async_test.json") # or use a pre-computed async model

# if products for GO Terms have not been computed and saved yet, compute them:
# model_async.fetch_all_go_term_products(web_download=True, run_async=True, recalculate=False, delay=0.0, run_async_options="v3", request_params={"rows":50000}, max_connections=60)
# model_async.create_products_from_goterms()
# model_async.save_model("diabetes_angio_4/model_async_test.json")

# ortholog product fetch
# model_async.fetch_ortholog_products(refetch=True, run_async=True, max_connections=15, req_delay=0.5, semaphore_connections=5) # semaphore_connections=10 works in 3min40s, semaphore_connections=15 results in 429:TooManyRequests errors
# model_async.save_model("diabetes_angio_4/model_async_test_ortholog_query.json")
