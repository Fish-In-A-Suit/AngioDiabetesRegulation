import dev_test_api_download
import util
import logging_config as cfg


def main():
#    dev_test_api_download.get_GO_genes_API("GO:1903670")
    terms_test = ['GO:1903670']
#    terms_angiogenesis_ids = util.get_array_terms("ANGIOGENESIS")
    dev_test_api_download.find_genes_related_to_GO_terms(terms_test)

if __name__ == '__main__':
    import logging.config
    logging.config.dictConfig(cfg.log_dict)
    main()