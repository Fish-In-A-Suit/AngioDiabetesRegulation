import ReverseLookup

import logging
logger = logging.getLogger(__name__)

def test_build():
    model.fetch_all_go_term_names_descriptions()
    model.fetch_all_go_term_products()
    model.create_products_from_goterms()
    print(model.__dict__)

def test_save():
    model.save_goterms_to_datafile("diabetes_angio_1/goterms.json")
    model.save_products_to_datafile("diabetes_angio_1/products.json")

def test_load():
    model.load_go_term_datafile("diabetes_angio_1/goterms.json")
    model.create_products_from_goterms()
    model.load_products_datafile("diabetes_angio_1/products.json")

def test_uniprotid():
    model.fetch_UniprotID_products()
    model.fetch_Uniprot_infos()

def test_scoring():
    model.score_products()

def test_report():
    report = ReverseLookup.ReportGenerator(model)
    #report.generate_detailed_design_report("diabetes_angio_1/detaileddesign.txt")
    report.generate_summary_report("diabetes_angio_1/summary.txt")

if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    logger.info("testing!")

    model = ReverseLookup.ReverseLookup.from_file("diabetes_angio_1/input.txt")

    test_load()
    #test_build()
    #test_save()
    test_uniprotid()
    #test_scoring()
    #test_report()
    test_save()