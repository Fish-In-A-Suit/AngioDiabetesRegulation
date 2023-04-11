import ReverseLookup
from wakepy import keepawake

import logging
logger = logging.getLogger(__name__)

def test_build():
    model.fetch_all_go_term_names_descriptions()
    model.fetch_all_go_term_products()
    model.create_products_from_goterms()
    print(model.__dict__)

def test_save():
    #model.save_goterms_to_datafile("diabetes_angio_1/goterms.json")
    #model.save_products_to_datafile("diabetes_angio_1/products.json")
    model.save_model("diabetes_angio_1/data.json")

def test_load_old():
    model.load_go_term_datafile("diabetes_angio_1/goterms.json")
    model.create_products_from_goterms()
    model.load_products_datafile("diabetes_angio_1/products.json")

def test_uniprot():
    #model.fetch_UniprotID_products()
    model.fetch_Uniprot_infos()

def test_product_scoring():
    model.score_products()

def test_mRNA():
    model.fetch_mRNA_sequences()

def test_miRNA():
    model.predict_miRNAs()

def test_miRNA_scoring():
    model.score_miRNAs()

def test_report():
    report = ReverseLookup.ReportGenerator(model, verbosity=3)
    #report.generate_detailed_design_report("diabetes_angio_1/detaileddesign.txt")
    report.general_report("diabetes_angio_1/general.txt")

def test_prune():
    model.prune_products()

if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    logger.info("testing!")

    model = ReverseLookup.ReverseLookup.from_input_file("diabetes_angio_1/input.txt", mod_name="V1")
    # model.save_model("diabetes_angio_1/data.json")

    # Fetch all GO term names and descriptions
    #model.fetch_all_go_term_names_descriptions()

    # Fetch all GO term products
    #model.fetch_all_go_term_products()

    # Create products from GO terms
    #model.create_products_from_goterms()

    # Fetch UniProt ID products
    #model.fetch_UniprotID_products()

    # Prune products
    #model.prune_products()

    # Fetch UniProt information
    #model.fetch_Uniprot_infos()

    # Score products
    #model.score_products()

    # Optional: Fetch mRNA sequences
    #model.fetch_mRNA_sequences()

    # Predict miRNAs
    #model.predict_miRNAs()

    # Score miRNAs
    #model.score_miRNAs()

    # Generate report
    #report = ReverseLookup.ReportGenerator(model, verbosity=3)
    #report.general_report("diabetes_angio_1/general.txt")

    model.compute_all("diabetes_angio_1/report.txt")

    # Save model
    model.save_model("diabetes_angio_1/data.json")

    # future runs:
    model = ReverseLookup.ReverseLookup.load_model("diabetes_angio_1/data.json")



    #with keepawake(keep_screen_awake=False):
        #test_load_old()
        #test_load()
        #test_build()
        #test_prune()
        #test_uniprot()
        #test_product_scoring()
        #test_mRNA()
        #test_miRNA()
        #test_miRNA_scoring()
        # test_report()
        #test_save()