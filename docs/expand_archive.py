# This Sphinx extension expands all the tar files from the archive into the
# build directory.

import os
import glob
from tarfile import TarFile

from sphinx.util.console import bold
from sphinx.util import logging
logger = logging.getLogger(__name__)


def setup(app):
    from sphinx.application import Sphinx
    if not isinstance(app, Sphinx):
        return
    app.connect('build-finished', expand_archive)


META_NOINDEX = '<meta name="robots" content="noindex, nofollow">'

LATEST_WARNING = """
<div class="admonition warning">
  <p class="first admonition-title">Note</p>
  <p class="last">
    This is an old version of the documentation. See
    <a href="/en/stable/index.html">http://docs.astropy.org/en/stable</a>
    for the latest version.
  </p>
</div>
"""


def expand_archive(app, exc):

    # Expand all tgz files directly into the build output

    archive_dir = os.path.join(app.builder.srcdir, 'archive')

    logger.info(bold('scanning {0} for archives...'.format(archive_dir)))

    for filename in sorted(glob.glob(os.path.join(archive_dir, '*.tgz'))):
        logger.info('   extracting {0}'.format(filename))
        tar_file = TarFile.open(filename)
        tar_file.extractall(app.builder.outdir)

    # Go through all html files in the built output and add the meta tag
    # to prevent search engines from crawling the pages. In addition, we also
    # rename some of the CSS classes otherwise the RTD javascript inserts
    # a message about this being an old version of the docs with an incorrect
    # link. We rename the classes (body -> old_body and document -> old_document)
    # to prevent this (since the RTD js uses these tags) and then insert the
    # warning in manually afterwards with the correct link.

    logger.info(bold('adding {0} tag to pages...').format(META_NOINDEX))

    for filename in glob.glob(os.path.join(app.builder.outdir, 'v*', '**', '*.html'), recursive=True):

        with open(filename, 'r') as f:
            content = f.read()

        # Insert meta noindex tag

        if META_NOINDEX not in content:
            if '<head>' in content:
                content = content.replace('<head>', '<head>{0}'.format(META_NOINDEX))
            else:
                raise Exception("Could not determine start of <head> section in {0}".format(filename))

        # Rename CSS classes to prevent automated warning

        content = content.replace('class="document"', 'class="old_document"')
        content = content.replace('class="body"', 'class="old_body"')

        # Insert manual warning about this being an old version of the docs

        pos = content.index('<div class="old_body"')
        pos = content.index('>', pos)

        content = content[:pos + 1] + LATEST_WARNING + content[pos + 1:]

        # Write file back out

        with open(filename, 'w') as f:
            f.write(content)

    # Finally we update the actual CSS files to change the body and document
    # class names.

    logger.info(bold('updating CSS files...'))

    for filename in glob.glob(os.path.join(app.builder.outdir, 'v*', '**',
                                           'bootstrap-astropy.css'), recursive=True):

        with open(filename, 'r') as f:
            content = f.read()

        # Note that the spaces after the strings are required to avoid matching
        # e.g. div.documentwrapper.
        content = content.replace('div.body ', 'div.old_body ')
        content = content.replace('div.document ', 'div.old_document ')

        with open(filename, 'w') as f:
            f.write(content)
