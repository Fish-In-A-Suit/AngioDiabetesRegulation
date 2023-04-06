import ReverseLookup

import logging
logger = logging.getLogger(__name__)

def test_build():
    model.fetch_all_go_term_names_descriptions()
    model.fetch_all_go_term_products()
    model.create_products_from_goterms()
    print(model.__dict__)

def test_save():
    model.create_products_from_goterms()
    print(model.__dict__)
    model.save_products_to_datafile("diabetes_angio_1/products.json")

def test_load():
    model.load_go_term_datafile("diabetes_angio_1/goterms.json")


if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    logger.info("testing!")

    model = ReverseLookup.ReverseLookup.from_file("diabetes_angio_1/input_2.txt")

    #test_build()
    test_load()
    test_save()