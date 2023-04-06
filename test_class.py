import logging
logger = logging.getLogger(__name__)

from full_analysis_classes import go_analysis_model

def main_test_add_key():
    append_dict = {"GO_id":"GO:1903589", "process":"angio", "direction":"+", "weight":1}
    append_list = ["GO:12345678", "diab", "+", 0]
    append_tuple = ("GO:123456789", "diab", "+", 0)
    append_no_GO = ["12345678", "diab", "+", 0]
    append_missing_tuple = ("diab", "+", 0)
    append_missing_dict = {"GO_id":"GO:1903589", "direction":"+", "weight":1}

    #model.add_goterm(append_list)
    #model.add_goterm(append_list)
    #model.add_goterm([append_tuple,append_dict])
    #model.add_goterm(append_missing_dict)
    print(model.go_terms_definitions)

def test_load_data():
    model.load_data("config")
    print(model.__dict__)



if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    logger.info("testing!")
    model = go_analysis_model(target_folder="diabetes_angio_1", config_filepath="diabetes_angio_1/input.txt")
    #main_test_add_key()
    test_load_data()