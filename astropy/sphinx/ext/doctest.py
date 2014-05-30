# **Please Note**: ``astropy.sphinx`` exists only for backward-compatibility
# purposes - it has now been moved to the separate astropy-helpers package,
# located at https://github.com/astropy/astropy-helpers. Any new development or
# bug fixes should be done there.
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This is a set of three directives that allow us to insert metadata
about doctests into the .rst files so the testing framework knows
which tests to skip.

This is quite different from the doctest extension in Sphinx itself,
which actually does something.  For astropy, all of the testing is
centrally managed from py.test and Sphinx is not used for running
tests.
"""
from docutils.nodes import literal_block
from sphinx.util.compat import Directive


class DoctestSkipDirective(Directive):
    has_content = True

    def run(self):
        code = '\n'.join(self.content)
        return [literal_block(code, code)]


class DoctestRequiresDirective(DoctestSkipDirective):
    # This is silly, but we really support an unbounded number of
    # optional arguments
    optional_arguments = 64


def setup(app):
    app.add_directive('doctest-requires', DoctestRequiresDirective)
    app.add_directive('doctest-skip', DoctestSkipDirective)
    app.add_directive('doctest-skip-all', DoctestSkipDirective)
