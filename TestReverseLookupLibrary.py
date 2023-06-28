from goreverselookuplib import ReverseLookup
from goreverselookuplib.AnnotationProcessor import GOAnnotiationsFile
from goreverselookuplib.Metrics import Metrics, adv_product_score, nterms, basic_mirna_score, binomial_test, fisher_exact_test
from goreverselookuplib.Report import ReportGenerator
from goreverselookuplib.Workflows import WorkflowOne

import os

"""
Example usage: First workflow (only adv_score and nterms score)
    # Pull GO data #
    model = ReverseLookup.load_model("diabetes_angio_2/data.json")
    model.fetch_all_go_term_names_descriptions()
    model.fetch_all_go_term_products(web_download=True)
    model.create_products_from_goterms()
    model.fetch_ortholog_products(refetch=False)
    model.fetch_product_infos(refetch=False)
    model.prune_products()
    model.save_model("diabetes_angio_2/data.json")

    # Score GO Term products #
    adv_score = adv_product_score(model)
    nterms_score = nterms(model)
    model.score_products([adv_score, nterms_score])
    model.save_model("diabetes_angio_2/data.json")

    # Score mRNA-miRNA #
    model.fetch_mRNA_sequences()
    model.predict_miRNAs()
    model.change_miRNA_overlap_treshold(0.6, True)
    basic_score = basic_mirna_score(model)
    model.score_miRNAs(basic_score)
    model.save_model("diabetes_angio_1/data.json")

    # Generate report #
    report = ReportGenerator(model, verbosity=3)
    report.general_report("diabetes_angio_1/general.txt", product_score=adv_score)
    
    
Example usage: Second workflow (adv_score, nterms_score, fisher_score, binomial_score)
# TODO: implement WorkflowTwo
    # Pull GO data #
    model = ReverseLookup.load_model("diabetes_angio_2/data.json")
    model.fetch_all_go_term_names_descriptions()
    model.fetch_all_go_term_products(web_download=True)
    model.create_products_from_goterms()
    model.fetch_ortholog_products(refetch=False)
    model.fetch_product_infos(refetch=False)
    model.prune_products()
    model.save_model("diabetes_angio_2/data.json")

    # Score GO Term products #
    adv_score = adv_product_score(model)
    nterms_score = nterms(model)
    goaf = GOAnnotiationsFile()
	binom_score = binomial_test(model, goaf)
	fisher_score = fisher_exact_test(model, goaf)
	model.score_products([adv_score, nterms_score, binom_score, fisher_score])
    model.save_model("diabetes_angio_2/data.json")

    # Score mRNA-miRNA #
    [TODO: implement mRNA-miRNA logic only for the valid GO Term products after scoring]
    // model.fetch_mRNA_sequences()
    // model.predict_miRNAs()
    // model.change_miRNA_overlap_treshold(0.6, True)
    // basic_score = basic_mirna_score(model)
    // model.score_miRNAs(basic_score)
    // model.save_model("diabetes_angio_1/data.json")

    # Generate report #
    [TODO: implement]
    // report = ReportGenerator(model, verbosity=3)
    // report.general_report("diabetes_angio_1/general.txt", product_score=adv_score)
"""

test_workflow = False

if test_workflow:
    workflow_one = WorkflowOne("diabetes_angio_3/data.json", "diabetes_angio_3", debug=True)
    workflow_one.run_workflow()
    exit()



# Define model from input file
# model = ReverseLookup.from_input_file("diabetes_angio_2/input_final.txt")
model = ReverseLookup.load_model("diabetes_angio_2/data.json")

#model.change_products_member_field("had_orthologs_computed", False)
#model.change_products_member_field("had_fetch_info_computed", False)
#model.save_model("diabetes_angio_2/data.json")
#exit()

# model._debug_shorten_GO_terms(5)

# Fetch all GO term names and descriptions
# model.fetch_all_go_term_names_descriptions()

# Fetch all GO term products
# model.fetch_all_go_term_products(web_download=True)

# Create products from GO terms
# model.create_products_from_goterms()

# model.save_model("diabetes_angio_2/data.json")

# Fetch human ortholog for products (either UniProtID, ENSG or genename)
# model.fetch_ortholog_products(refetch=False)

# model.prune_products()

# model.save_model("diabetes_angio_2/data.json")

# Fetch product information (from UniprotAPI or EnsemblAPI)
# model.fetch_product_infos(refetch=False)

# Prune products
# model.prune_products()

# model.save_model("diabetes_angio_2/data.json")

# Score products
adv_score = adv_product_score(model)
nterms_score = nterms(model)
# model.score_products([adv_score, nterms_score])
goaf = GOAnnotiationsFile()
binom_score = binomial_test(model, goaf)
fisher_score = fisher_exact_test(model, goaf)
model.score_products([adv_score, nterms_score, binom_score, fisher_score])
model.save_model("diabetes_angio_2/data_06-06-2023.json")
model.perform_statistical_analysis(filepath="diabetes_angio_2/stat_analysis")

#model.save_model("diabetes_angio_1/data.json")

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
#report.general_report("diabetes_angio_1/general.txt", product_score=fisher_score)

# Save model
#model.save_model("diabetes_angio_1/data.json")

# -- One-time call functions for any post-processing --
# model.change_products_member_field("had_orthologs_computed", True)
# model.change_products_member_field("had_fetch_info_computed", False)