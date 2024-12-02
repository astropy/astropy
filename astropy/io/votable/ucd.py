# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This file contains routines to verify the correctness of UCD strings.
"""

# STDLIB
import re

# LOCAL
from astropy.utils import data

__all__ = ["check_ucd", "parse_ucd"]


class UCDWords:
    """
    Manages a list of acceptable UCD words.

    Works by reading in a data file exactly as provided by IVOA.  This
    file resides in data/ucd1p-words.txt.
    """

    def __init__(self):
        self._primary = set()
        self._secondary = set()
        self._descriptions = {}
        self._capitalization = {}

        with data.get_pkg_data_fileobj("data/ucd1p-words.txt", encoding="ascii") as fd:
            for line in fd.readlines():
                if line.startswith("#"):
                    continue

                type, name, descr = (x.strip() for x in line.split("|"))
                name_lower = name.lower()
                if type in "QPEVC":
                    self._primary.add(name_lower)
                if type in "QSEVC":
                    self._secondary.add(name_lower)
                self._descriptions[name_lower] = descr
                self._capitalization[name_lower] = name

    def is_primary(self, name):
        """
        Returns True if *name* is a valid primary name.
        """
        return name.lower() in self._primary

    def is_secondary(self, name):
        """
        Returns True if *name* is a valid secondary name.
        """
        return name.lower() in self._secondary

    def get_description(self, name):
        """
        Returns the official English description of the given UCD
        *name*.
        """
        return self._descriptions[name.lower()]

    def normalize_capitalization(self, name):
        """
        Returns the standard capitalization form of the given name.
        """
        return self._capitalization[name.lower()]


_ucd_singleton = None


def parse_ucd(ucd, check_controlled_vocabulary=False, has_colon=False):
    """
    Parse the UCD into its component parts.

    Parameters
    ----------
    ucd : str
        The UCD string

    check_controlled_vocabulary : bool, optional
        If `True`, then each word in the UCD will be verified against
        the UCD1+ controlled vocabulary, (as required by the VOTable
        specification version 1.2), otherwise not.

    has_colon : bool, optional
        If `True`, the UCD may contain a colon (as defined in earlier
        versions of the standard).

    Returns
    -------
    parts : list
        The result is a list of tuples of the form:

            (*namespace*, *word*)

        If no namespace was explicitly specified, *namespace* will be
        returned as ``'ivoa'`` (i.e., the default namespace).

    Raises
    ------
    ValueError
        if *ucd* is invalid
    """
    global _ucd_singleton
    if _ucd_singleton is None:
        _ucd_singleton = UCDWords()

    if has_colon:
        m = re.search(r"[^A-Za-z0-9_.:;\-]", ucd)
    else:
        m = re.search(r"[^A-Za-z0-9_.;\-]", ucd)
    if m is not None:
        raise ValueError(f"UCD has invalid character '{m.group(0)}' in '{ucd}'")

    word_component_re = r"[A-Za-z0-9][A-Za-z0-9\-_]*"
    word_re = rf"{word_component_re}(\.{word_component_re})*"

    parts = ucd.split(";")
    words = []
    for i, word in enumerate(parts):
        colon_count = word.count(":")
        if colon_count == 1:
            ns, word = word.split(":", 1)
            if not re.match(word_component_re, ns):
                raise ValueError(f"Invalid namespace '{ns}'")
            ns = ns.lower()
        elif colon_count > 1:
            raise ValueError(f"Too many colons in '{word}'")
        else:
            ns = "ivoa"

        if not re.match(word_re, word):
            raise ValueError(f"Invalid word '{word}'")

        if ns == "ivoa" and check_controlled_vocabulary:
            if i == 0:
                if not _ucd_singleton.is_primary(word):
                    if _ucd_singleton.is_secondary(word):
                        raise ValueError(
                            f"Secondary word '{word}' is not valid as a primary word"
                        )
                    else:
                        raise ValueError(f"Unknown word '{word}'")
            else:
                if not _ucd_singleton.is_secondary(word):
                    if _ucd_singleton.is_primary(word):
                        raise ValueError(
                            f"Primary word '{word}' is not valid as a secondary word"
                        )
                    else:
                        raise ValueError(f"Unknown word '{word}'")

        try:
            normalized_word = _ucd_singleton.normalize_capitalization(word)
        except KeyError:
            normalized_word = word
        words.append((ns, normalized_word))

    return words


def check_ucd(ucd, check_controlled_vocabulary=False, has_colon=False):
    """
    Returns False if *ucd* is not a valid `unified content descriptor`_.

    Parameters
    ----------
    ucd : str
        The UCD string

    check_controlled_vocabulary : bool, optional
        If `True`, then each word in the UCD will be verified against
        the UCD1+ controlled vocabulary, (as required by the VOTable
        specification version 1.2), otherwise not.

    has_colon : bool, optional
        If `True`, the UCD may contain a colon (as defined in earlier
        versions of the standard).

    Returns
    -------
    valid : bool
    """
    if ucd is None:
        return True

    try:
        parse_ucd(
            ucd,
            check_controlled_vocabulary=check_controlled_vocabulary,
            has_colon=has_colon,
        )
    except ValueError:
        return False
    return True
