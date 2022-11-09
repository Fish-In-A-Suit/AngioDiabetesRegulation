log_dict = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'standard': {
            'format': '%(asctime)s [%(levelname)s] %(name)s: %(message)s'
        },
    },
    'handlers': {
        'stream': {
            'level': 'INFO',
            'formatter': 'standard',
            'class': 'logging.StreamHandler',
        },
        'file': {
            'level': 'INFO',
            'filename': './log_output/test_json_dump.log',
            'class': 'logging.FileHandler',
            'formatter': 'standard'
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