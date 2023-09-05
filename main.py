# The main class
from goreverselookuplib import ReverseLookup
from goreverselookuplib.AnnotationProcessor import GOAnnotiationsFile
from goreverselookuplib.Workflows import WorkflowOne, WorkflowTwo
from goreverselookuplib.AnnotationProcessor import GOApi
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

Cacher.init()
workflow = WorkflowTwo(input_file_fpath="chronic_infl_cancer_1/input_03-09-2023.txt", save_folder_dir="chronic_infl_cancer_1")
workflow.run_workflow()

#model = ReverseLookup.load_model("chronic_infl_cancer_1/data_31-08-2023.json")
#goaf = GOAnnotiationsFile()
#fisher_test = fisher_exact_test(model, goaf)
#model.score_products([fisher_test])
#model.perform_statistical_analysis(filepath="chronic_infl_cancer_1/statistically_relevant_genes_31-08-2023(1).json")
#model.save_model("chronic_infl_cancer_1/data_31-08-2023(1).json")


