#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""webquery: Get the output of a specified URL.

This program gets the output of a specified URL.  The user
specifies the host and URL as well as any parameters that might
be used by the URL.  The output of URL is returned either
to standard output or to a specified file.  Five special keyword
parameters may be specified: url, host, method, port and file.
All other parameters are passed to the URL.

Webquery does no reformatting or error checking.  Typically
the returned data will be an HTML formatted document.  If an
error is encountered, the returned data may be an HTML formatted
error message.

Usage::

    % webquery.py url=URL host=HOST method=METHOD port=PORT file=FILE
        [key=VALUE key=VALUE]

The default URL is the null string. The default host is the local
host.  The method and port keywords may be used
to override the defaults of POST and 80 respectively.
The file keyword specifies an output file to which the response is
to be sent.  If not specified, output is sent to standard output.

Additional keywords are appended to the URL request as
keyword values for a forms request.  Webquery will appropriately
escape special characters.

Has the same command-line interface as the Perl script webquery.pl,
but is vastly simpler since it uses much higher-level Python modules
as an interface.

:Authors: Perry Greenfield, Michael Droettboom

:Organization: Space Telescope Science Institute

See Also
--------
ftp://legacy.gsfc.nasa.gov/heasarc/software/web_batch/webquery.pl

"""
# STDLIB
import socket
import sys
import urllib
import urllib2

__all__ = ['webquery_open', 'webquery', 'webget_open', 'webget']

URLLIB2_HAS_TIMEOUT = (sys.hexversion >= 0x02060000)


class WebQueryError(Exception):  # pragma: no cover
    pass


def webquery_open(args=(), **kw):
    """
    Return a read-only file descriptor to read the results of a web
    query.

    Keywords for query may be specified as a sequence of pairs in args
    or as keywords.  Special keywords that define the URL include:

        * host (default 'localhost')
        * url (default null)
        * method (default 'POST')
        * port (default 80)
        * timeout (default None)

    Additional keywords are passed as parameters to the query.

    If a parameter keyword has a list as its value, the parameter is
    included multiple times in the query, once for each argument.

    """
    args = list(args)
    for key, value in kw.iteritems():
        args.append((key, value))
    port = 80
    method = "POST"
    url = ""
    host = urllib.localhost()
    query = []
    timeout = None
    for key, value in args:
        if key == "port":
            port = int(value)
        elif key == "method":
            method = value.upper()
        elif key == "url":
            url = value
        elif key == "host":
            host = value
        elif key == "timeout":
            timeout = value
        elif value is None:
            query.append(urllib.quote(key))
        elif isinstance(value, list):
            qkey = urllib.quote(key)
            for v in value:
                query.append('{}={}'.format(qkey, urllib.quote_plus(str(v))))
        else:
            query.append('{}={}'.format(
                urllib.quote(key), urllib.quote_plus(str(value))))
    query = '&'.join(query)

    if url[:1] == "/":
        # don't add an extra slash (purely for aesthetic purposes)
        url = "http://{}:{}{}".format(host, port, url)
    else:
        url = "http://{}:{}/{}".format(host, port, url)

    if not query:
        query = None
    elif method == "GET":
        url = "{}?{}".format(url, query)
        query = None

    if URLLIB2_HAS_TIMEOUT:
        return urllib2.urlopen(url, query, timeout)
    else:
        # This is the old way to set a socket timeout prior to Python
        # 2.6.  NOTE THIS IS NOT THREADSAFE
        old_timeout = socket.getdefaulttimeout()
        socket.setdefaulttimeout(timeout)
        try:
            req = urllib2.urlopen(url, query)
        finally:
            socket.setdefaulttimeout(old_timeout)
        return req


def webquery(args=(), **kw):
    """
    Write output of a specified URL to stdout or file.

    Keywords for query may be specified as a sequence of pairs in args
    or as keywords.  Special keywords that define the URL include:

        * host (default 'localhost')
        * url (default null)
        * method (default 'POST')
        * port (default 80)
        * timeout (default None)
        * file (output filename or file handle; default sys.stdout)

    Additional keywords are passed as parameters to the query.

    If a parameter keyword has a list as its value, the parameter is
    included multiple times in the query, once for each argument.

    """
    args = list(args)
    for key, value in kw.iteritems():
        args.append((key, value))
    outfile = sys.stdout
    outargs = []
    for key, value in args:
        if key == "file":
            outfile = value
        else:
            outargs.append((key, value))

    close_outfile = False
    if isinstance(outfile, basestring):
        outfile = open(outfile, "w")
        close_outfile = True

    inurl = webquery_open(outargs)

    try:
        s = inurl.read(102400)
        while s:
            outfile.write(s)
            s = inurl.read(102400)
    finally:
        inurl.close()
        if close_outfile:
            outfile.close()


def webget_open(url, timeout=None, method='GET', **keywords):
    """
    Simplified version of webquery that presumes the url is well formed with
    some optional GET key/value pairs.

    The url can be comprised of the host, path, and fixed key=value
    pairs. Keyword/value pairs are simply appended (thus the url
    should explicitly have '?' in it already. If there are already
    parameters part of the url, the url string should end with '&' if
    more are expected).

    Returns a read-only file-like object to stream over the net.  This
    object should be closed when no longer needed.

    """
    if len(keywords) and not (url.endswith('?') or url.endswith('&')):
        raise WebQueryError("url should already end with '?' or '&'")

    query = []
    for key, value in keywords.iteritems():
        value = keywords[key]
        if value is None:
            query.append(urllib.quote(key))
        elif isinstance(value, list):
            qkey = urllib.quote(key)
            for v in value:
                query.append('{}={}'.format(qkey, urllib.quote_plus(str(v))))
        else:
            query.append('{}={}'.format(
                urllib.quote(key), urllib.quote_plus(str(value))))
    query = '&'.join(query)

    if method == 'GET':
        url += query
        query = None
    elif method == 'POST':
        pass
    else:
        raise WebQueryError("method kwarg must be 'GET' or 'POST'")

    if URLLIB2_HAS_TIMEOUT:
        return urllib2.urlopen(url, query, timeout)
    else:
        # This is the old way to set a socket timeout prior to Python
        # 2.6.  NOTE THIS IS NOT THREADSAFE
        old_timeout = socket.getdefaulttimeout()
        socket.setdefaulttimeout(timeout)
        try:
            req = urllib2.urlopen(url, query)
        finally:
            socket.setdefaulttimeout(old_timeout)
        return req


def webget(url, file=None, timeout=None, method='GET', **keywords):
    """
    Simplified version of webquery that presumes the url is well formed with
    some optional GET key/value pairs.

    The url can be comprised of the host, path, and fixed key=value
    pairs. Keyword/value pairs are simply appended (thus the url
    should explicitly have '?' in it already. If there are already
    parameters part of the url, the url string should end with '&' if
    more are expected).  Returns the info obtained as part of the
    retrieval.

    """
    inurl = webget_open(url, timeout=timeout, method=method, **keywords)

    close_outfile = False
    if isinstance(file, basestring):
        file = open(file, "w")
        close_outfile = True

    try:
        s = inurl.read(102400)
        while s:
            file.write(s)
            s = inurl.read(102400)
    finally:
        inurl.close()
        if close_outfile:
            file.close()

    return inurl.info()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise WebQueryError(__doc__)
    else:
        arglist = []
        for arg in sys.argv[1:]:
            f = arg.split('=', 1)
            if len(f) == 1:
                arglist.append((arg, None))
            else:
                arglist.append(f)
        webquery(arglist)
