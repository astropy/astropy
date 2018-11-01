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


def expand_archive(app, exc):

    archive_dir = os.path.join(app.builder.srcdir, 'archive')

    logger.info(bold('scanning {0} for archives...'.format(archive_dir)))

    for filename in sorted(glob.glob(os.path.join(archive_dir, '*.tgz'))):
        logger.info('   extracting {0}'.format(filename))
        tar_file = TarFile.open(filename)
        tar_file.extractall(app.builder.outdir)

    logger.info(bold('...done'))
