import dev_test_api_download
import logging_config as cfg


def main():
    dev_test_api_download.get_GO_genes_API("GO:1903670")

if __name__ == '__main__':
    import logging.config
    logging.config.dictConfig(cfg.log_dict)
    main()