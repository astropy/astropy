import datetime
import getpass
import logging
import os
import re
import sys
import xmlrpclib

from ConfigParser import ConfigParser

try:
    from docutils.core import publish_parts
except ImportError:
    print >> sys.stderr, \
           'docutils is required to convert the PyFITS changelog to HTML ' \
           'for updating the PyFITS homepage\n\n' \
           'Try `pip install docutils` or `easy_install docutils`.'
    sys.exit(1)

from zest.releaser.choose import version_control
from zest.releaser.utils import get_last_tag, ask


log = None


PYFITS_HOMEPAGE_BASE_URL = \
    'http://www.stsci.edu/resources/software_hardware/pyfits'
# These are the pages to run find/replace of the version number on
PYFITS_HOMEPAGE_SUBPAGES = ['localProductDescription', 'Download']

# The website will only have final release up on it, so we can use a simplified
# version regexp
VERSION_RE = re.compile(r'(?P<MAJOR>\d+)\.(?P<MINOR>\d+)(?:\.(?P<MICRO>\d+))?')

# This is the format used to search for/replace the previous version
# This is based on simply a manual analysis of where the PyFITS version number
# appears on the website; note that the version alone shouldn't be used since
# other version strings (i.e. Python versions) can appear on the site
# NOTE: The MICRO version number is appended to the format later, since it is
# optional if the micro format is 0
SEARCH_VERSION_RE_FORMAT = (r'(?P<prefix>v|V|[vV]ersion\s+|pyfits-)'
                             '%(major)s\.%(minor)s(?:\s*\((?P<date>.+)\))?')

DATE_FORMAT = '%B %d %Y'


class ReleaseManager(object):
    def __init__(self):
        self.vcs = version_control()
        self.history_lines = []
        self.previous_version = ''

    def prereleaser_before(self, data):
        """Set tag_svn_revision to False."""

        if data['name'] != 'pyfits':
            return

        global log
        log = logging.getLogger('prerelease')

        # TODO: This bit should belong to a generic release hook somewhere in
        # the stsci namespace, or even should maybe be suggested as a built in
        # action of zest.releaser
        def callback(section, option, value, lineno, line):
            if section == 'egg_info' and option == 'tag_svn_revision':
                log.info('Disabling tag_svn_revision in setup.cfg')
                return 'tag_svn_revision = False\n'
            else:
                return line

        config_parser('setup.cfg', callback)

    def prereleaser_after(self, data):
        """Before preforming the release, get the previously released version
        from the latest tag in version control.
        """

        if data['name'] != 'pyfits':
            return

        self.previous_version = get_last_tag(self.vcs)
        self.history_lines = data['history_lines']

    def postreleaser_before(self, data):
        """Restore tag_svn_revision"""

        if data['name'] != 'pyfits':
            return

        global log
        log = logging.getLogger('postrelease')

        def callback(section, option, value, lineno, line):
            if section == 'egg_info' and option == 'tag_svn_revision':
                return 'tag_svn_revision = True'
            else:
                return line

        config_parser('setup.cfg', callback)

    def postreleaser_after(self, data):
        """Used to update the PyFITS website.

        TODO: If at any point we get a Windows build machine that we can remote
        into, use this as a point to create Windows builds as well.
        """

        if data['name'] != 'pyfits':
            return

        if not ask('Update PyFITS homepage'):
            return


        previous_version = raw_input(
            'Enter previous version [%s]: ' % self.previous_version).strip()
        if not previous_version:
            previous_version = self.previous_version

        new_version = raw_input(
            'Enter new version [%s]: ' % data['new_version']).strip()
        if not new_version:
            new_version = data['new_version']

        username = raw_input(
                'Enter your Zope username [%s]: ' % getpass.getuser()).strip()
        if not username:
            username = getpass.getuser()

        password = getpass.getpass(
            'Enter your Zope password (password will not be displayed): ')

        match = VERSION_RE.match(previous_version)
        if not match:
            log.error('Previous version (%s) is invalid: Version must be in '
                      'the MAJOR.MINOR[.MICRO] format.' % previous_version)
            sys.exit()

        micro = int(match.group('MICRO')) if match.group('MICRO') else 0

        previous_version = (int(match.group('MAJOR')),
                            int(match.group('MINOR')), micro)

        match = VERSION_RE.match(new_version)
        if not match:
            log.error('New version (%s) is invalid: Version must be in '
                      'the MAJOR.MINOR[.MICRO] format.' % new_version)
            sys.exit()

        micro = int(match.group('MICRO')) if match.group('MICRO') else 0

        new_version = (int(match.group('MAJOR')), int(match.group('MINOR')),
                       micro)

        # This is the regular expression to search for version replacement
        format_args = {'major': str(previous_version[0]),
                       'minor': str(previous_version[1])}
        if previous_version[2] != 0:
            # Append the micro version after the minor version if nonzero
            format_args['minor'] += r'\.%d' % previous_version[2]
        search_version_re = re.compile(SEARCH_VERSION_RE_FORMAT % format_args)

        new_version_str = '.'.join((str(new_version[0]), str(new_version[1])))
        if new_version[2] != 0:
            new_version_str += '.%d' % new_version[2]

        def version_replace(match):
            repl = match.group('prefix') + new_version_str
            if match.group('date'):
                today = datetime.datetime.today().strftime(DATE_FORMAT)
                repl += ' (%s)' % today
            return repl

        # Go ahead and do the find/replace on supported subpages
        for page in PYFITS_HOMEPAGE_SUBPAGES:
            try:
                url = os.path.join(PYFITS_HOMEPAGE_BASE_URL, page)
                proxy = _ZopeProxy(url, username, password)
                content = proxy.retrieve()
                content = search_version_re.sub(version_replace, content)
                proxy.update(content)
            except Exception, e:
                continue

        # Update the release notes
        parts = publish_parts('\n'.join(self.history_lines),
                              writer_name='html')
        # Get just the body of the HTML and convert headers to <h3> tags
        # instead of <h1> (there might be a 'better' way to do this, but
        # this is a simple enough case to suffice for our purposes
        content = parts['html_body']

        # A quickie regexp--no good for general use, but should work fine
        # in this case; this will prevent replacement of the <h1> tag in
        # the title, but will take care of all the others
        content = re.sub(r'<h1>([^<]+)</h1>', r'<h3>\1</h3>', content)

        try:
            url = os.path.join(PYFITS_HOMEPAGE_BASE_URL, 'release')
            proxy = _ZopeProxy(url, username, password)
            # And upload...
            proxy.update(content)
        except Exception, e:
            pass


# TODO: This is also a handy utility that could probaby be used elsewhere
def config_parser(filename, callback):
    """This is a very simplified config file parser that can update a config
    file in-place so that order and comments are preserved.

    The callback should be a function that takes a section, option, value,
    line number, and raw line as its input (the current config section the
    parser is in, the current option, its value, the line number of the config
    file, and the actual line string the parser is on).

    The callback function is called for each line of the file.  For multi-line
    option values the function is still called once for each line of the
    option, so the callback needs to know how to handle these if it desires to.
    A line in which section, option, and value are None is either a comment or
    a blank line.

    The return value of the callback function should be the raw line to output
    or an iterable of lines to output.  In most cases the callback function
    will just return the same line that was passed in.
    """

    config = open(filename).readlines()
    new_config = []
    current_section = None
    current_option = None
    updated = False

    for lineno, line in enumerate(config):
        match = ConfigParser.SECTCRE.match(line)
        if match:
            current_section = match.group('header')
            current_option = None
            section, option, value = current_section, None, None
        else:
            if re.match(r'^\s*#', line) or not line.strip():
                section, option, value = None, None, None
            elif re.match(r'^\s+', line):
                # A new line in the current option
                section, option, value = (current_section, current_option,
                                          line.strip())
            else:
                option, value = (item.strip() for item in line.split('=', 1))
                section = current_section
                current_option = option
        lines = callback(section, option, value, lineno, line)
        if lines != line:
            updated = True
        if isinstance(lines, basestring):
            new_config.append(lines)
        else:
            new_config.extend(lines)

    if updated:
        open(filename, 'w').writelines(new_config)


class _ZopeProxy(object):
    """This is a simple class for handling retriving and updating of pages on a
    Zope2 site.  This only handles updates to static content, and not
    directories or anything like that.
    """

    def __init__(self, url, username=None, password=None):
        if username and password:
            protocol, rest = url.split('://', 1)
            self.url = '%s://%s:%s@%s' % (protocol, username, password, rest)
            self.masked_url = '%s//%s:%s@%s' % (protocol, username, '*' * 8,
                                                rest)
        else:
            self.url = self.masked_url = url

        self.proxy = None

    def connect(self):
        if self.proxy is not None:
            return
        try:
            self.proxy = xmlrpclib.ServerProxy(self.url)
        except Exception, e:
            # TODO: Catch bad authentication and let the user enter a new
            # username/password
            if log:
                log.error('Failed to connect to %s: %s' %
                          (self.masked_url, str(e)))
            raise

    def retrieve(self):
        """Retrieves the static page contents at the proxy's URL."""

        self.connect()
        if log:
            log.info('Retrieving %s...' % self.masked_url)
        try:
            return self.proxy.document_src()
        except Exception, e:
            if log:
                log.error('Failed to download content at %s: %s' %
                          (self.masked_url, str(e)))
            raise

    def update(self, content):
        """Updates the static page content at the proxy's URL."""

        self.connect()
        if log:
             log.info('Updating %s...' % self.masked_url)
        try:
            self.proxy.manage_upload(content)
        except Exception, e:
            if log:
                log.error('Failed to update content at %s: %s' %
                          (self.masked_url, str(e)))
            raise


releaser = ReleaseManager()
