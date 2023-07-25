"""
loggin utilities
"""
import typing
from contextlib import contextmanager
import logging
import time


FORMAT = "%(asctime)s - %(pathname)s:%(lineno)d - %(name)s - %(levelname)s: %(message)s"


def progressbar(
    iterable: typing.Iterable, logger: typing.Optional[logging.Logger] = None
):
    """
    Args:
        iterable
        logger: optional

    Examples:
        >>> import logging
        >>> import rna
        >>> import sys
        >>> sys.modules['tqdm'] = None
        >>> log = logging.getLogger(__name__)

        >>> a = range(3)
        >>> for value in rna.log.progressbar(a, logger=log):
        ...     _ = value * 3

    """
    if logger is None:
        logger = logging
    try:
        # tqdm not required for the module to work
        from tqdm import tqdm as progressor  # pylint: disable=import-outside-toplevel

        tqdm_exists = True
    except ImportError:

        def progressor(iterable):
            """
            dummy function. Does nothing
            """
            return iterable

        tqdm_exists = False
    try:
        n_total = len(iterable)
    except TypeError:
        n_total = None

    for i in progressor(iterable):
        if not tqdm_exists:
            if n_total is None:
                logger.info("Progress: item %d", i)
            else:
                logger.info("Progress: {%d} / {%d}", i, n_total)
        yield i


@contextmanager
def timeit(
    process_name: str,
    *args,
    logger: typing.Optional[logging.Logger] = None,
    precision=1,
) -> None:
    """
    Context manager for autmated timeing

    Args:
        process_name (str): message to customize the log message
        logger
        precision (int): show until 10^-<precision> digits

    Examples:
        >>> import rna
        >>> with rna.log.timeit("Logger Name"):
        ...     # Some very long calculation
        ...     pass
    """
    if logger is None:
        logger = logging
    start_time = time.time()
    logger.log(logging.INFO, "Timing Process: " + process_name + " ...", *args)

    yield

    logger.log(
        logging.INFO,
        "\t\t\t\t\t\t... Process Duration: %1." + str(int(precision)) + "f s",
        time.time() - start_time,
    )


class ColorFormatter(logging.Formatter):
    """
    Formatter with colors corresponding to severity level.

    This can be quickly used like so:
    """

    GREY = "\x1b[38;5;249m"
    WHITE = "\x1b[38;21m"
    YELLOW = "\x1b[33;21m"
    RED = "\x1b[31;21m"
    BOLD_RED = "\x1b[31;1m"
    RESET = "\x1b[0m"

    COLORS = {
        logging.DEBUG: GREY,
        logging.INFO: WHITE,
        logging.WARNING: YELLOW,
        logging.ERROR: RED,
        logging.CRITICAL: BOLD_RED,
    }

    def format(self, record):
        return self.COLORS.get(record.levelno) + super().format(record) + self.RESET


class TerminalHandler(logging.StreamHandler):
    """
    Stream handler with ColorFormatter set already.

    Args:
        fmt: if str, passed to ColorFormatter, else logging.Formatter expected.


    Examples:

        >>> import logging
        >>> import rna.log
        >>> logger = logging.getLogger()

        Add a stream handler to your logger (colored by default).
        >>> logger.addHandler(rna.log.TerminalHandler())

        Now error logs will be red
        >>> logger.error("error logs")

    """

    def __init__(
        self,
        *args,
        fmt: typing.Optional[typing.Union[str, logging.Formatter]] = None,
        **kwargs,
    ):
        if fmt is None:
            fmt = FORMAT
        if isinstance(fmt, str):
            formatter = ColorFormatter(fmt)
        else:
            formatter = fmt

        level = kwargs.pop("level", None)
        super().__init__(*args, **kwargs)
        self.setFormatter(formatter)
        if level is not None:
            # same signature as the logging.Handler (StreamHandler has different signature)
            self.setLevel(level)
