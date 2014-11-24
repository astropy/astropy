# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Script support for validating a VO file.
"""

from ... import __version__

def main(args=None):
    from . import table
    from astropy.utils.compat import argparse

    parser = argparse.ArgumentParser(
        description=("Check a VOTable file for compliance to the "
                     "VOTable specification"
                     "This script is part of the Astropy package " + __version__ +  "."
                     "More documentation can be found here: "
                     "http://astropy.readthedocs.org/en/latest/io/votable/index.htmll"))
    parser.add_argument(
        'filename', nargs=1, help='Path to VOTable file to check')
    args = parser.parse_args(args)

    table.validate(args.filename[0])
