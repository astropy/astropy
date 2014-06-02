# **Please Note**: ``astropy.sphinx`` exists only for backward-compatibility
# purposes - it has now been moved to the separate astropy-helpers package,
# located at https://github.com/astropy/astropy-helpers. Any new development or
# bug fixes should be done there.
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from ....tests.helper import pytest
pytest.importorskip('sphinx')  # skips these tests if sphinx not present


class FakeEnv(object):
    """
    Mocks up a sphinx env setting construct for automodapi tests
    """
    def __init__(self, **kwargs):
        for k, v in kwargs.iteritems():
            setattr(self, k, v)


class FakeBuilder(object):
    """
    Mocks up a sphinx builder setting construct for automodapi tests
    """
    def __init__(self, **kwargs):
        self.env = FakeEnv(**kwargs)


class FakeApp(object):
    """
    Mocks up a `sphinx.application.Application` object for automodapi tests
    """
    def __init__(self, srcdir, automodapipresent=True):
        self.builder = FakeBuilder(srcdir=srcdir)
        self.info = []
        self.warnings = []
        self._extensions = []
        if automodapipresent:
            self._extensions.append('astropy.sphinx.ext.automodapi')

    def info(self, msg, loc):
        self.info.append((msg, loc))

    def warn(self, msg, loc):
        self.warnings.append((msg, loc))


ams_to_asmry_str = """
Before

.. automodsumm:: astropy.sphinx.ext.automodsumm
    :p:

And After
"""

ams_to_asmry_expected = """.. currentmodule:: astropy.sphinx.ext.automodsumm

.. autosummary::
    :p:

    Automoddiagram
    Automodsumm
    automodsumm_to_autosummary_lines
    generate_automodsumm_docs
    process_automodsumm_generation
    setup"""


def test_ams_to_asmry(tmpdir):
    from ..automodsumm import automodsumm_to_autosummary_lines

    fi = tmpdir.join('automodsumm.rst')
    fi.write(ams_to_asmry_str)

    fakeapp = FakeApp(srcdir='')
    resultlines = automodsumm_to_autosummary_lines(str(fi), fakeapp)

    assert '\n'.join(resultlines) == ams_to_asmry_expected
