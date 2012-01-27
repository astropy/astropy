# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Utilities for managing the Astropy project - e.g. github tricks, release
scripts, etc.
"""

__all__ = ['issue_to_pr']


def issue_to_pr(issuenum, srcbranch, repo='astropy', targetuser='astropy',
                targetbranch='master', username=None, pw=None,
                baseurl='https://api.github.com'):
    """
    Attaches code to an issue, converting a regular issue into a pull request.

    Parameters
    ----------
    issuenum : int
        The issue number (in `targetuser`/`repo`) onto which the code should be
        attached.
    srcbranch : str
        The branch (in `username`/`repo`) the from which to attach the code.
        After this is complete, further updates to this branch will be passed
        on to the pull request.
    repo : str
        The name of the repository (the same-name repo should be present for
        both `targetuser` and `username`).
    targetuser : str
        The name of the user/organization that has the issue.
    targetbranch : str
        The name of the branch (in `targetuser`/`repo`) that the pull request
        should merge into.
    username : str or None
        The name of the user sourcing the pull request.  If None, the caller
        will be prompted for their name at the terminal.
    pw : str or None
        The password of the user sourcing the pull request.  If None, the
       caller will be prompted for their name at the terminal (text will be
       hidden).
   baseurl : str
       The URL to use to access the github site (including protocol).

    .. warning::
        Be cautious supplying `pw` as a string - if you do this in an ipython
        session, for example, it will be logged in the input history, revealing
        your password in clear text.  When in doubt, leave it as `None`,
        as this will securely prompt you for your password.

    Returns
    -------
    response : str
        The json-decoded response from the github server.
    error message : str, optional
        If present, indicates github responded with an HTTP error.

    """

    import urllib
    import urllib2
    import getpass
    import json

    if username is None:
        username = raw_input('Username: ')
    if pw is None:
        pw = getpass.getpass()

    data = {'issue': str(issuenum),
            'head': username + ':' + srcbranch,
            'base': targetbranch}

    datajson = json.dumps(data)

    suburl = 'repos/{user}/{repo}/pulls'.format(user=targetuser, repo=repo)
    url = urllib.basejoin(baseurl, suburl)

    req = urllib2.Request(url)
    req.add_data(datajson)
    _add_basic_auth_header(req, username, pw)

    try:
        response = urllib2.urlopen(req)
    except urllib2.HTTPError, e:
        print 'HTTP Error', e
        res = e.fp.read()
        return json.loads(res), str(e)
    res = response.read()
    return json.loads(res)


def _add_basic_auth_header(req, username, pw):
    from base64 import b64encode

    upwstr = username + ':' + pw
    upwstrenc = b64encode(upwstr.encode('utf-8')).strip().decode('utf-8')

    req.add_header('Authorization', 'Basic ' + upwstrenc)
