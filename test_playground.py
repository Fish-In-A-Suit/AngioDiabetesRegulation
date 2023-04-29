from goreverselookuplib import ReverseLookup
from goreverselookuplib.AnnotationProcessor import GOAnnotiationsFile

import logging

from math import log10

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(funcName)s: %(message)s')

# Create a file handler
file_handler = logging.FileHandler('./log_output/test_json_dump.log')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(logging.Formatter('%(asctime)s [%(filename)s:%(lineno)s] - [%(funcName)s] %(message)s'))

# Add handlers to the root logger
#logging.getLogger().addHandler(file_handler)

logger = logging.getLogger(__name__)


model = ReverseLookup.load_model("diabetes_angio_2/data_web_full_scores.json")

#test execution speed
#import timeit
#goaf = GOAnnotiationsFile()
#elapsed_time = timeit.timeit(lambda: goaf.get_products("GO:1903589"), number=100)
#print(elapsed_time)

#goaf = GOAnnotiationsFile()
#print(goaf.get_all_terms_for_product("ACAT2"))
#print(len(goaf.get_all_terms()))

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

    for product in model.products:
        #processes = model.target_processes
        processes = [     
            {
                "process": "diabetes",
                "direction": "+"
            },
            #{
            #    "process": "angio",
            #    "direction": "+"
            #},
            {
                "process": "obesity",
                "direction": "+"
            }
            ]
        
    #    if (all(float(product.scores["fisher_test"][f"{process['process']}{process['direction']}"]["pvalue_corr"]) < 0.05 for process in processes)):
    #        print(f"Found significant product: {product.genename}")

        if "fisher_test" in product.scores:
            if (all(float(product.scores["fisher_test"][f"{process['process']}{process['direction']}"].get("pvalue_corr", 1)) < 0.05 for process in processes) and
                all(float(product.scores["fisher_test"][f"{process['process']}{'+' if process['direction'] == '-' else '-'}"].get("pvalue_corr", 1)) >= 0.05 for process in processes)):
                print(f"Found significant product: {product.genename}")
    exit()
model.prune_products()

model.save_model("diabetes_angio_1/data.json")