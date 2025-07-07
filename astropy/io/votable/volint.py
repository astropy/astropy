# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Script support for validating a VO file.
"""

import argparse


def main(args=None):
    from . import table

    parser = argparse.ArgumentParser(
        description="Check a VOTable file for compliance to the VOTable specification"
    )
    # TODO: pass color and suggest_on_error as kwargs when PYTHON_LT_14 is dropped
    parser.color = True
    parser.suggest_on_error = True

    parser.add_argument("filename", nargs=1, help="Path to VOTable file to check")
    args = parser.parse_args(args)

    table.validate(args.filename[0])
