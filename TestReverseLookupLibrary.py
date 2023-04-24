from goreverselookuplib import ReverseLookup
from goreverselookuplib.AnnotationProcessor import GOAnnotiationsFile
from goreverselookuplib.Metrics import Metrics, adv_product_score, nterms, basic_mirna_score, binomial_test, fisher_exact_test
from goreverselookuplib.Report import ReportGenerator

import os

# Define model from input file
#model = ReverseLookup.from_input_file("diabetes_angio_1/input.txt")
model = ReverseLookup.load_model("diabetes_angio_1/data.json")

# Fetch all GO term names and descriptions
# model.fetch_all_go_term_names_descriptions()

# Fetch all GO term products
#model.fetch_all_go_term_products()

# Create products from GO terms
#model.create_products_from_goterms()

#odel.save_model("diabetes_angio_1/data.json")

# Fetch human ortholog for products (either UniProtID, ENSG or genename)
# model.fetch_ortholog_products()

# model.prune_products()

# model.save_model("diabetes_angio_1/data.json")

# Fetch product information (from UniprotAPI or EnsemblAPI)
# model.fetch_product_infos()

# Prune products
# model.prune_products()

# model.save_model("diabetes_angio_1/data.json")

# Score products
#adv_score = adv_product_score(model)
#nterms_score = nterms(model)
goaf = GOAnnotiationsFile()
binom_score = binomial_test(model, goaf)
fisher_score = fisher_exact_test(model, goaf)
model.score_products([fisher_score])

model.save_model("diabetes_angio_1/data.json")

# Optional: Fetch mRNA sequences
# model.fetch_mRNA_sequences()

# Predict miRNAs
# model.predict_miRNAs()

# Score miRNAs
#model.change_miRNA_overlap_treshold(0.6, True)
#basic_score = basic_mirna_score(model)
#model.score_miRNAs(basic_score)

# Generate report
#report = ReportGenerator(model, verbosity=3)
# report.general_report("diabetes_angio_1/general.txt") # this is bugged
#report.general_report("diabetes_angio_1/general.txt", product_score=binom_score)

# Save model
#model.save_model("diabetes_angio_1/data.json")