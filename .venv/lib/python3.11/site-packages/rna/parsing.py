"""
Collection of standard parent parsers.
"""
import argparse
import logging
import rna.log


# pylint:disable=protected-access
class LoggingAction(argparse._StoreAction):
    """
    Call basic config with namespace.level
    """

    def __call__(self, parser, namespace, values, option_string=None):
        super().__call__(parser, namespace, values, option_string=option_string)
        # logging.basicConfig(level=namespace.level)
        root = logging.getLogger()
        root.setLevel(namespace.log_level)

        formatter = rna.log.ColorFormatter(namespace.log_format or rna.log.FORMAT)

        terminal_handler = logging.StreamHandler()
        terminal_handler.setFormatter(formatter)
        root.addHandler(terminal_handler)

        if namespace.log_file:
            file_handler = logging.FileHandler(filename=namespace.log_file, mode="w")
            file_handler.setFormatter(formatter)
            root.addHandler(file_handler)


class LogParser(argparse.ArgumentParser):
    """
    ArgumentParser concerning logging setup.
    Actions are already setting handlers and formatters.

    Examples:
        >>> import argparse
        >>> import rna.parsing
        >>> parser = argparse.ArgumentParser(parents=[rna.parsing.LogParser()])
        >>> args, _ = parser.parse_known_args("--log_level 42".split())
        >>> args.log_level
        42
    """

    def __init__(self, add_help=False, **kwargs):
        super().__init__(add_help=add_help, **kwargs)
        self.add_argument(
            "--log_level",
            type=int,
            default=logging.INFO,
            help="Level of root logger.",
        )
        self.add_argument(
            "--log_file",
            type=str,
            help="FileHandler to be added to root logger",
        )
        self.add_argument(
            "--log_format",
            type=str,
            help="Logging format. Default: rna.log.FORMAT",
            action=LoggingAction,
        )
