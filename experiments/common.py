import logging

def set_verbose():
    logging.getLogger("").setLevel(logging.DEBUG)
    # Enhance the format
    fmt = logging.Formatter(
        "%(levelname)-8s | %(asctime)s | %(name)-20s | %(message)s")
    logging.getLogger("").handlers[0].setFormatter(fmt)
