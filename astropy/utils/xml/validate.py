# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Functions to do XML schema and DTD validation.  At the moment, this
makes a subprocess call to xmllint.  This could use a Python-based
library at some point in the future, if something appropriate could be
found.
"""


import os
import subprocess


def validate_schema(filename, schema_file):
    """
    Validates an XML file against a schema or DTD.

    Parameters
    ----------
    filename : str
        The path to the XML file to validate

    schema_file : str
        The path to the XML schema or DTD

    Returns
    -------
    returncode, stdout, stderr : int, str, str
        Returns the returncode from xmllint and the stdout and stderr
        as strings
    """
    base, ext = os.path.splitext(schema_file)
    if ext == ".xsd":
        schema_part = "--schema"
    elif ext == ".dtd":
        schema_part = "--dtdvalid"
    else:
        raise TypeError("schema_file must be a path to an XML Schema or DTD")

    p = subprocess.Popen(
        ["xmllint", "--noout", "--nonet", schema_part, schema_file, filename],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdout, stderr = p.communicate()

    if p.returncode == 127:
        raise OSError("xmllint not found, so can not validate schema")
    elif p.returncode < 0:
        from astropy.utils.misc import signal_number_to_name

        raise OSError(
            f"xmllint was terminated by signal '{signal_number_to_name(-p.returncode)}'"
        )

    return p.returncode, stdout, stderr
