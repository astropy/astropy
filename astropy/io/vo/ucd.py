"""
This file contains routines to verify the correctness of UCD strings.
"""

from __future__ import with_statement, absolute_import

#STDLIB
import os
import re

class UCDWords:
    """
    A class to manage the list of acceptable UCD words.  Works by
    reading in a data file exactly as provided by IVOA.  This file
    resides in data/ucd1p-words.txt.
    """
    def __init__(self):
        self._primary = set()
        self._secondary = set()
        self._descriptions = {}
        self._capitalization = {}

        ucd_words_filepath = os.path.join(
            os.path.dirname(__file__),
            "data", "ucd1p-words.txt")
        with open(ucd_words_filepath, 'r') as fd:
            for line in fd.readlines():
                type, name, descr = [x.strip() for x in line.split('|')]
                name_lower = name.lower()
                if type in 'QPEV':
                    self._primary.add(name_lower)
                if type in 'QSEV':
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
    Parse the UCD into its component parts.  The result is a list of
    tuples of the form:

      (*namespace*, *word*)

    If no namespace was explicitly specified, *namespace* will be
    returned as ``'ivoa'`` (i.e., the default namespace).

    If *check_controlled_vocabulary* is ``True``, then each word in
    the UCD will be verified against the UCD1+ controlled vocabulary,
    (as required by the VOTable specification version 1.2), otherwise
    not.

    Will raise `ValueError` if *ucd* is invalid.
    """
    global _ucd_singleton
    if _ucd_singleton is None:
        _ucd_singleton = UCDWords()

    if ucd is None:
        return True

    if has_colon:
        m = re.search('[^A-Za-z0-9_.:;\-]', ucd)
    else:
        m = re.search('[^A-Za-z0-9_.;\-]', ucd)
    if m is not None:
        raise ValueError("UCD has invalid character '%s' in '%s'" %
                         (m.group(0), ucd))

    word_component_re = '[A-Za-z0-9][A-Za-z0-9\-_]*'
    word_re = '%s(\.%s)*' % (word_component_re, word_component_re)

    parts = ucd.split(';')
    words = []
    for i, word in enumerate(parts):
        colon_count = word.count(':')
        if colon_count == 1:
            ns, word = word.split(':', 1)
            if not re.match(word_component_re, ns):
                raise ValueError("Invalid namespace '%s'" % ns)
            ns = ns.lower()
        elif colon_count > 1:
            raise ValueError("Too many colons in '%s'" % word)
        else:
            ns = 'ivoa'

        if not re.match(word_re, word):
            raise ValueError("Invalid word '%s'" % word)

        if ns == 'ivoa' and check_controlled_vocabulary:
            if i == 0:
                if not _ucd_singleton.is_primary(word):
                    if _ucd_singleton.is_secondary(word):
                        raise ValueError(
                            "Secondary word '%s' is not valid as a primary word" %
                            word)
                    else:
                        raise ValueError("Unknown word '%s'" % word)
            else:
                if not _ucd_singleton.is_secondary(word):
                    if _ucd_singleton.is_primary(word):
                        raise ValueError(
                            "Primary word '%s' is not valid as a secondary word" %
                            word)
                    else:
                        raise ValueError("Unknown word '%s'" % word)

        try:
            normalized_word = _ucd_singleton.normalize_capitalization(word)
        except KeyError:
            normalized_word = word
        words.append((ns, normalized_word))

    return words


def check_ucd(ucd, check_controlled_vocabulary=False):
    """
    Returns False if *ucd* is not a valid `unified content
    descriptor`_.

    If *check_controlled_vocabulary* is ``True``, then each word in
    the UCD will be verified against the UCD1+ controlled vocabulary,
    (as required by the VOTable specification version 1.2), otherwise
    not.
    """
    try:
        parse_ucd(ucd, check_controlled_vocabulary=check_controlled_vocabulary)
    except ValueError as e:
        return False
    return True

