# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Script support for validating the WCS keywords in a FITS file.
"""


def main(args=None):
    import argparse

    from . import wcs

    parser = argparse.ArgumentParser(
        description=(
            "Check the WCS keywords in a FITS file for compliance against the standards"
        )
    )
    # TODO: pass color and suggest_on_error as kwargs when PYTHON_LT_14 is dropped
    parser.color = True
    parser.suggest_on_error = True

    parser.add_argument("filename", nargs=1, help="Path to FITS file to check")
    args = parser.parse_args(args)

    print(wcs.validate(args.filename[0]))
