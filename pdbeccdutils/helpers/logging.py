import logging


def set_up_logger(logger_name):
    """Set up logger of the application

    Args:
        logger_name (str): Name of the logger to be used

    Returns:
        logging.logger: Application log.
    """
    frm = "[%(asctime)-15s]  %(message)s"
    logging.basicConfig(
        level=logging.DEBUG, format=frm, datefmt="%a, %d %b %Y %H:%M:%S"
    )

    logging.getLogger().disabled = True
    logging.getLogger("neo4j").disabled = True

    return logging.getLogger(logger_name)


logger = set_up_logger("pdbeccdutils")
