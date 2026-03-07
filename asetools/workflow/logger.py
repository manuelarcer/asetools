# Logger

import logging
import datetime

def configure_logging(
    *,
    project_logger: str = "asetools.manager",
    file_prefix:    str = "run",
    level:          int = logging.INFO,
):
    logfile = f"{file_prefix}_{datetime.datetime.now():%Y%m%d_%H%M}.log"

    logging.basicConfig(
        level    = level,
        format   = "%(asctime)s [%(name)s %(levelname)s] %(message)s",
        handlers = [
            logging.FileHandler(logfile),
            logging.StreamHandler()
        ]
    )

    # Optionally enforce a minimum level on your project’s root logger
    proj = logging.getLogger(project_logger)
    proj.setLevel(level)