log_dict = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'standard': {
            'format': '%(asctime)s [%(levelname)s] %(funcName)s: %(message)s'
        },
        'extended': {
            'format': '%(asctime)s [%(filename)s:%(lineno)s] - [%(funcName)s] %(message)s'
        }
    },
    'handlers': {
        'stream': {
            'level': 'INFO',
            'formatter': 'standard',
            'class': 'logging.StreamHandler',
        },
        'file': {
            'level': 'DEBUG',
            'filename': './log_output/test_json_dump.log',
            'class': 'logging.FileHandler',
            'formatter': 'extended',
            'mode': 'w',
        }
    },
    'loggers': {
        '': {
            'handlers': ['stream','file'],
            'level': 'DEBUG',
            'propagate': True
        },
    }
}