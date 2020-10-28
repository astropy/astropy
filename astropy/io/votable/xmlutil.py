# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Various XML-related utilities
"""


# ASTROPY
from astropy.logger import log
from astropy.utils import data
from astropy.utils.xml import check as xml_check
from astropy.utils.xml import validate

# LOCAL
from .exceptions import (warn_or_raise, vo_warn, W02, W03, W04, W05)


__all__ = [
    'check_id', 'fix_id', 'check_token', 'check_mime_content_type',
    'check_anyuri', 'validate_schema'
    ]


def check_id(ID, name='ID', config=None, pos=None):
    """
    Raises a `~astropy.io.votable.exceptions.VOTableSpecError` if *ID*
    is not a valid XML ID_.

    *name* is the name of the attribute being checked (used only for
    error messages).
    """
    if (ID is not None and not xml_check.check_id(ID)):
        warn_or_raise(W02, W02, (name, ID), config, pos)
        return False
    return True


def fix_id(ID, config=None, pos=None):
    """
    Given an arbitrary string, create one that can be used as an xml id.

    This is rather simplistic at the moment, since it just replaces
    non-valid characters with underscores.
    """
    if ID is None:
        return None
    corrected = xml_check.fix_id(ID)
    if corrected != ID:
        vo_warn(W03, (ID, corrected), config, pos)
    return corrected


_token_regex = r"(?![\r\l\t ])[^\r\l\t]*(?![\r\l\t ])"


def check_token(token, attr_name, config=None, pos=None):
    """
    Raises a `ValueError` if *token* is not a valid XML token.

    As defined by XML Schema Part 2.
    """
    if (token is not None and not xml_check.check_token(token)):
        return False
    return True


def check_mime_content_type(content_type, config=None, pos=None):
    """
    Raises a `~astropy.io.votable.exceptions.VOTableSpecError` if
    *content_type* is not a valid MIME content type.

    As defined by RFC 2045 (syntactically, at least).
    """
    if (content_type is not None and
        not xml_check.check_mime_content_type(content_type)):
        warn_or_raise(W04, W04, content_type, config, pos)
        return False
    return True


def check_anyuri(uri, config=None, pos=None):
    """
    Raises a `~astropy.io.votable.exceptions.VOTableSpecError` if
    *uri* is not a valid URI.

    As defined in RFC 2396.
    """
    if (uri is not None and not xml_check.check_anyuri(uri)):
        warn_or_raise(W05, W05, uri, config, pos)
        return False
    return True


def validate_schema(filename, version='1.1'):
    """
    Validates the given file against the appropriate VOTable schema.

    Parameters
    ----------
    filename : str
        The path to the XML file to validate

    version : str, optional
        The VOTABLE version to check, which must be a string \"1.0\",
        \"1.1\", \"1.2\" or \"1.3\".  If it is not one of these,
        version \"1.1\" is assumed.

        For version \"1.0\", it is checked against a DTD, since that
        version did not have an XML Schema.

    Returns
    -------
    returncode, stdout, stderr : int, str, str
        Returns the returncode from xmllint and the stdout and stderr
        as strings
    """
    if version not in ('1.0', '1.1', '1.2', '1.3'):
        log.info(f'{filename} has version {version}, using schema 1.1')
        version = '1.1'

    if version in ('1.1', '1.2', '1.3'):
        schema_path = data.get_pkg_data_filename(
            f'data/VOTable.v{version}.xsd')
    else:
        schema_path = data.get_pkg_data_filename(
            'data/VOTable.dtd')

    return validate.validate_schema(filename, schema_path)
