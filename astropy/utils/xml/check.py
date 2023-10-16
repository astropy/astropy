# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A collection of functions for checking various XML-related strings for
standards compliance.
"""


import re
import urllib.parse


def check_id(ID):
    """
    Returns `True` if *ID* is a valid XML ID.
    """
    return re.match(r"^[A-Za-z_][A-Za-z0-9_\.\-]*$", ID) is not None


def fix_id(ID):
    """
    Given an arbitrary string, create one that can be used as an xml
    id.  This is rather simplistic at the moment, since it just
    replaces non-valid characters with underscores.
    """
    if re.match(r"^[A-Za-z_][A-Za-z0-9_\.\-]*$", ID):
        return ID
    if len(ID):
        corrected = ID
        if not len(corrected) or re.match("^[^A-Za-z_]$", corrected[0]):
            corrected = "_" + corrected
        corrected = re.sub(r"[^A-Za-z_]", "_", corrected[0]) + re.sub(
            r"[^A-Za-z0-9_\.\-]", "_", corrected[1:]
        )
        return corrected
    return ""


_token_regex = r"(?![\r\l\t ])[^\r\l\t]*(?![\r\l\t ])"


def check_token(token):
    """
    Returns `True` if *token* is a valid XML token, as defined by XML
    Schema Part 2.
    """
    return (
        token == ""
        or re.match(r"[^\r\n\t ]?([^\r\n\t ]| [^\r\n\t ])*[^\r\n\t ]?$", token)
        is not None
    )


def check_mime_content_type(content_type):
    """
    Returns `True` if *content_type* is a valid MIME content type
    (syntactically at least), as defined by RFC 2045.
    """
    ctrls = "".join(chr(x) for x in range(0x20))
    token_regex = f'[^()<>@,;:\\"/[\\]?= {ctrls}\x7f]+'
    return (
        re.match(rf"(?P<type>{token_regex})/(?P<subtype>{token_regex})$", content_type)
        is not None
    )


def check_anyuri(uri):
    """
    Returns `True` if *uri* is a valid URI as defined in RFC 2396.
    """
    if (
        re.match(
            (
                r"(([a-zA-Z][0-9a-zA-Z+\-\.]*:)?/{0,2}[0-9a-zA-Z;"
                r"/?:@&=+$\.\-_!~*'()%]+)?(#[0-9a-zA-Z;/?:@&=+$\.\-_!~*'()%]+)?"
            ),
            uri,
        )
        is None
    ):
        return False
    try:
        urllib.parse.urlparse(uri)
    except Exception:
        return False
    return True
