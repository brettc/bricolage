from python_log_indenter import IndentedLoggerAdapter
import logging
import os
import inspect

# TODO: Testing this
# import coloredlogs
# coloredlogs.install(level='DEBUG',
#                     fmt="%(levelname)-8s | %(asctime)s | %(message)s")
# coloredlogs.install(level='DEBUG')


def get_logger(fname=None):
    """Pass in a file name or we'll try to detect it."""

    # Magically get the filename from the calling function
    if fname is None:
        caller_frame = inspect.stack()[1][0]
        fname = caller_frame.f_globals["__file__"]

    # Strip the beginning and the extension
    head_tail = os.path.split(fname)
    base_ext = os.path.splitext(head_tail[1])

    log_name = base_ext[0]

    # Trim it max 10 characters
    log_name = log_name[:10]

    # Now wrap it and return it
    return IndentedLoggerAdapter(logging.getLogger(log_name))


def set_logging(verbose=False):
    logging.basicConfig(
        format="%(levelname)-8s | %(asctime)s | %(message)s", level=logging.INFO
    )
    # fmt = logging.Formatter("%(levelname)-8s | %(asctime)s | %(message)s")
    # print logging.getLogger("sqlalchemy").handlers
    #
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

        # We need to do this separately. This is equivalent to setting "echo"
        # on the database creation.
        logging.getLogger("sqlalchemy.engine").setLevel(logging.INFO)
