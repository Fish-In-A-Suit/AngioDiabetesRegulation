# This Python file enables the user to check differences between offline and online query algorithms for 
# product ortholog search. The final results are stored in root/debug_files/offline-online-ortholog-mismatches.txt
from goreverselookuplib import ReverseLookup
from goreverselookuplib.AnnotationProcessor import GOAnnotiationsFile
from goreverselookuplib.Workflows import WorkflowOne
from goreverselookuplib.AnnotationProcessor import GOApi
from goreverselookuplib.JsonUtil import SimpleNamespaceUtil, JsonToClass
from goreverselookuplib.AnnotationProcessor import HumanOrthologFinder, UniProtAPI, EnsemblAPI, GOAnnotiationsFile

import logging

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s [%(levelname)s] %(funcName)s: %(message)s')
logger = logging.getLogger(__name__)

human_ortholog_finder = HumanOrthologFinder()
uniprot_api = UniProtAPI()
ensembl_api = EnsemblAPI()
goaf = GOAnnotiationsFile()

model =  ReverseLookup.load_model("diabetes_angio_4/model_async_test.json")

#test_products = model_async.products[0:5]
test_products = model.products
for product in test_products:
    product.fetch_ortholog(human_ortholog_finder, uniprot_api, ensembl_api, goaf, prefer_goaf=True, _d_compare_goaf=True)

mismatches = 0
with open("debug_files/offline-online-ortholog-mismatches.txt", "w") as f:
    for product in test_products:
        if product._d_offline_online_ortholog_mismatch == True:
            mismatches += 1
            f.write(f"{product._d_offline_online_ortholog_mismatch_values}\n")
    
logger.info(f"number of online-offline mismatches: {mismatches}")