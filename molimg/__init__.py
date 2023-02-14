import logging.config

LOGGING_CONFIG = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "default": {
            "format": ("%(asctime)s:%(name)s:%(process)d:%(lineno)d " 
                       "%(levelname)s %(message)s"),
            "datefmt": "%Y-%m-%d %H:%M:%S",
        },
        "simple": {
            "format": "%(message)s",
        },
        "json": {
            "class": "pythonjsonlogger.jsonlogger.JsonFormatter",
            "format": """
                message: %(message)s
                name: %(name)s
                filename: %(filename)s
                funcName: %(funcName)s
                levelname: %(levelname)s
                lineno: %(lineno)d
                asctime: %(asctime)s
            """,
            "datefmt": "%Y-%m-%d %H:%M:%S",
        }
    },
    "handlers": {
        "verbose_output": {
            "formatter": "simple",
            "level": "DEBUG",
            "class": "logging.StreamHandler",
            "stream": "ext://sys.stdout",
        },
        "json": {
            "formatter": "json",
            "level": "INFO",
            "class": "logging.StreamHandler",
            "stream": "ext://sys.stdout",
        }
    },
    "loggers": {
        "molimg": {
            "level": "INFO",
            "handlers": [
                "json"
            ],
        },
        "root": {
            "level": "DEBUG",
            "handlers": [
                "json"
            ],
        },
    },
}
logging.config.dictConfig(LOGGING_CONFIG)
logger = logging.getLogger(__name__)
