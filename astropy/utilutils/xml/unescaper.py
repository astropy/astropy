# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""URL unescaper functions."""

# STDLIB
from xml.sax import saxutils


__all__ = ['unescape_all']

# This is DIY
_bytes_entities = {b'&amp;': b'&', b'&lt;': b'<', b'&gt;': b'>',
                   b'&amp;&amp;': b'&', b'&&': b'&', b'%2F': b'/'}
_bytes_keys = [b'&amp;&amp;', b'&&', b'&amp;', b'&lt;', b'&gt;', b'%2F']

# This is used by saxutils
_str_entities = {'&amp;&amp;': '&', '&&': '&', '%2F': '/'}
_str_keys = ['&amp;&amp;', '&&', '&amp;', '&lt;', '&gt;', '%2F']


def unescape_all(url):
    """Recursively unescape a given URL.

    .. note:: '&amp;&amp;' becomes a single '&'.

    Parameters
    ----------
    url : str or bytes
        URL to unescape.

    Returns
    -------
    clean_url : str or bytes
        Unescaped URL.

    """
    if isinstance(url, bytes):
        func2use = _unescape_bytes
        keys2use = _bytes_keys
    else:
        func2use = _unescape_str
        keys2use = _str_keys
    clean_url = func2use(url)
    not_done = [clean_url.count(key) > 0 for key in keys2use]
    if True in not_done:
        return unescape_all(clean_url)
    else:
        return clean_url


def _unescape_str(url):
    return saxutils.unescape(url, _str_entities)


def _unescape_bytes(url):
    clean_url = url
    for key in _bytes_keys:
        clean_url = clean_url.replace(key, _bytes_entities[key])
    return clean_url
