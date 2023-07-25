"""
Main entry point for the project when run from the command line
i.e. via $ python -m rna
"""

import sys
import argparse
import rna


class SomeAction(argparse.Action):
    """Some actions."""

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        print(
            f"Example action invoked by manage in namespace: {namespace} with values {values}"
            " and option string {option_string}."
        )
        setattr(namespace, self.dest, values)

    def showcase_dummy(self):
        """
        You can define a method to expose functionality of the class
        """
        print(self)


def manage(args_):
    """Example function."""
    print("Managing!")
    print(args_.x * args_.y)


def parse_args(args_):
    """Parse args."""
    # create the top-level parser
    parser = argparse.ArgumentParser(prog="rna app")
    parser.add_argument(
        "--version",
        action="version",
        version="v" + rna.__version__,
        help="Show program's version number and exit",
    )
    parser = argparse.ArgumentParser(prog="rna app")

    # subparsers
    subparsers = parser.add_subparsers(help="sub-command help")

    # create the parser for the "test" command
    example_sub_parser = subparsers.add_parser("manage", help="manage something")
    example_sub_parser.add_argument("-x", type=int, default=1)
    example_sub_parser.add_argument("-y", type=float, default=42.0)
    example_sub_parser.set_defaults(func=manage)

    # If no arguments were used, print base-level help with possible commands.
    if len(args_) == 0:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    args_ = parser.parse_args(args_)
    # let argparse do the job of calling the appropriate function after
    # argument parsing is complete
    return args_.func(args_)


if __name__ == "__main__":
    _ = parse_args(sys.argv[1:])
