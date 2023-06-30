from goreverselookuplib import ReverseLookup
from goreverselookuplib.AnnotationProcessor import GOAnnotiationsFile
from goreverselookuplib.Workflows import WorkflowOne
from goreverselookuplib.AnnotationProcessor import GOApi
import asyncio
import aiohttp
from goreverselookuplib.JsonUtil import SimpleNamespaceUtil, JsonToClass

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
"""
# This async approach works and is 1000x faster than synchronous
# TODO: IMPLEMENT THIS !!!!!
# goterms_selection = model.goterms[0:3]
goterms_selection = src_model.goterms
async def test_goterm_async():
    # Testing asnyc threading for server requests
    tasks = []
    for goterm in goterms_selection:
        # if goterm.name == None or goterm.description == None: # if goterm.name or description don't exist, then attempt fetch
        # goterm.fetch_name_description_async(api)
        task = asyncio.create_task(goterm.fetch_name_description_async(api))
        tasks.append(task)

    await asyncio.gather(*tasks)

asyncio.run(test_goterm_async())

# access updated goterms
for goterm in goterms_selection:
    print(goterm.name)

# goterms_diff = src_model.compare_to(reference_model, "goterms", compare_subfields=['name', 'description'])
#products_diff = src_model.compare_to(reference_model, "products")
#print(goterms_diff)
"""

### ReverseLookup model comparison testing
# load two models, change some data
#model1 = ReverseLookup.load_model("diabetes_angio_2/data_06-06-2023.json")
#model2 = ReverseLookup.load_model("diabetes_angio_2/data_06-06-2023.json")
#model2.goterms[0].name = "faulty-name"
#model2.goterms[3].description = "faulty-description"
#model2.goterms[32].products = ["faulty-product"]
#model2.products[40].scores["adv_score"] = -10000
#model2.products[30].scores["nterms"] = -10000
#model2.products[15].scores["binomial_test"]["diabetes-"]["pvalue_corr"] = -1000
#model2.products[20].scores["fisher_test"]["diabetes-"]["pvalue_corr"] = -1000

#goterms_diff = model1.compare_to(model2, "goterms", compare_subfields=['name','description','weight','products'])
#products_diff = model1.compare_to(model2, "products", compare_subfields=['id_synonyms','description','uniprot_id','ensg_id','enst_id','refseq_nt_id','mRNA','scores_adv-score','scores_nterms','scores_binomial-test','scores_fisher-test'])

## logger.info(f"Entire differences json: {goterms_diff}")
## logger.info(f"Entire differences json: {products_diff}")

#i = 0
#for diff in goterms_diff:
#    print("\n")
#    i+=1
#    logger.info(f"[{i}]: difference at: {diff}, difference: {goterms_diff[diff]}")
    
#for diff in products_diff:
#    print("\n")
#    i+=1
#    logger.info(f"[{i}]: difference at: {diff}, difference: {products_diff[diff]}")

# compare goterms and products both in one function call
#all_diff = model1.compare_to(model2, compare_field="", compare_subfields=[])
#for diff in all_diff:
#    print("\n")
#    i+=1
#    logger.info(f"[{i}]: difference at: {diff}, difference: {all_diff[diff]}")


# test async vs. reference model
model_async = ReverseLookup.from_input_file("diabetes_angio_4/input.txt")
#async def fetch_name_desc_test_async():
#    tasks = []
#    for goterm in model_async.goterms:
#        task = asyncio.create_task(goterm.fetch_name_description_async(api))
#        tasks.append(task)
#    await asyncio.gather(*tasks)
#asyncio.run(fetch_name_desc_test_async())
# above async code was placed into ReverseLookup and can now be called as:
# model_async.fetch_all_go_term_names_descriptions(run_async=True)

# model_reference = ReverseLookup.load_model("diabetes_angio_2/data_06-06-2023.json")

# all_diff = model_async.compare_to(model_reference, compare_field="goterms", compare_subfields=["name", "description"])
#i = 0
#for diff in all_diff:
#    print("\n")
#    i+=1
#    logger.info(f"[{i}]: difference at: {diff}, difference: {all_diff[diff]}")
#logger.info(f"total difference = {all_diff}")
### --- End of ReverseLookup model comparison testing

### Start of ReverseLookup fetch_all_goterm_products async testing
# model_async = ReverseLookup.from_input_file("diabetes_angio_4/input.txt")
# model_async.fetch_all_go_term_names_descriptions(run_async=True)
# model_reference = ReverseLookup.load_model("diabetes_angio_2/data_06-06-2023.json")
model_async.fetch_all_go_term_products(web_download=True, run_async=True, recalculate=False)