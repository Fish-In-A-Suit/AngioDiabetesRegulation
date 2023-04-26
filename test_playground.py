from goreverselookuplib import ReverseLookup
from goreverselookuplib.AnnotationProcessor import GOAnnotiationsFile

import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(funcName)s: %(message)s')

# Create a file handler
file_handler = logging.FileHandler('./log_output/test_json_dump.log')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(logging.Formatter('%(asctime)s [%(filename)s:%(lineno)s] - [%(funcName)s] %(message)s'))

# Add handlers to the root logger
#logging.getLogger().addHandler(file_handler)

logger = logging.getLogger(__name__)


model = ReverseLookup.load_model("diabetes_angio_2/data.json")

#test execution speed
#import timeit
#goaf = GOAnnotiationsFile()
#elapsed_time = timeit.timeit(lambda: goaf.get_products("GO:1903589"), number=100)
#print(elapsed_time)

#goaf = GOAnnotiationsFile()
#print(goaf.get_all_terms_for_product("ACAT2"))
#print(len(goaf.get_all_terms()))

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
        #{
        #    "process": "obesity",
        #    "direction": "+"
        #}
        ]
    
#    if (all(float(product.scores["fisher_test"][f"{process['process']}{process['direction']}"]["pvalue_corr"]) < 0.05 for process in processes)):
#        print(f"Found significant product: {product.genename}")


    if (all(float(product.scores["fisher_test"][f"{process['process']}{process['direction']}"]["pvalue_corr"]) < 0.05 for process in processes) and
        all(float(product.scores["fisher_test"][f"{process['process']}{'+' if process['direction'] == '-' else '-'}"]["pvalue_corr"]) >= 0.05 for process in processes)):
        print(f"Found significant product: {product.genename}")
exit()
model.prune_products()

model.save_model("diabetes_angio_1/data.json")