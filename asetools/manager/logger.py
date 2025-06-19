# Logger

import logging
import datetime

def configure_logger(logger_name: str = 'logger'):
    """
    Configure the logger to log to both a file and stdout.
    The log file will be named with the current date and time.
    """
    log_filename = f"{logger_name}_{datetime.datetime.now():%Y%m%d_%H%M}.log"
    logging.basicConfig(
        level    = logging.INFO,
        format   = "%(asctime)s [%(levelname)s] %(message)s",
        handlers = [
            logging.FileHandler(log_filename),
            logging.StreamHandler()      # also print to stdout
        ]
    )
    return logging.getLogger(logger_name)