# Licensed under a 3-clause BSD style license - see LICENSE.rst
from pytest import importorskip
importorskip('sphinx')  # skips these tests if sphinx not present

class FakeConfig(object):
    """
    Mocks up a sphinx configuration setting construct for automodapi tests
    """
    def __init__(self, **kwargs):
        for k, v in kwargs.iteritems():
            setattr(self, k, v)


class FakeApp(object):
    """
    Mocks up a `sphinx.application.Application` object for automodapi tests
    """
    def __init__(self, **configs):
        self.config = FakeConfig(**configs)
        self.info = []
        self.warnings = []

    def info(self, msg, loc):
        self.info.append((msg, loc))

    def warn(self, msg, loc):
        self.warnings.append((msg, loc))


am_replacer_str = """
This comes before

.. automodapi:: astropy.sphinx.ext.tests.test_automodapi
{options}

This comes after
"""

am_replacer_basic_expected = """
This comes before

astropy.sphinx.ext.tests.test_automodapi Module
-----------------------------------------------

.. automodule:: astropy.sphinx.ext.tests.test_automodapi

Functions
^^^^^^^^^

.. automodsumm:: astropy.sphinx.ext.tests.test_automodapi
    :functions-only:
    :toctree: _generated/

Classes
^^^^^^^

.. automodsumm:: astropy.sphinx.ext.tests.test_automodapi
    :classes-only:
    :toctree: _generated/

Class Inheritance Diagram
^^^^^^^^^^^^^^^^^^^^^^^^^

.. automod-diagram:: astropy.sphinx.ext.tests.test_automodapi
    :private-bases:

This comes after
"""


def test_am_replacer_basic():
    """
    Tests replacing an ".. automodapi::" with the automodapi no-option
    template
    """
    from ..automodapi import automodapi_replace

    fakeapp = FakeApp(automodapi_toctreedirnm='_generated/')
    result = automodapi_replace(am_replacer_str.format(options=''), fakeapp)

    assert result == am_replacer_basic_expected

am_replacer_noinh_expected = """
This comes before

astropy.sphinx.ext.tests.test_automodapi Module
-----------------------------------------------

.. automodule:: astropy.sphinx.ext.tests.test_automodapi

Functions
^^^^^^^^^

.. automodsumm:: astropy.sphinx.ext.tests.test_automodapi
    :functions-only:
    :toctree: _generated/

Classes
^^^^^^^

.. automodsumm:: astropy.sphinx.ext.tests.test_automodapi
    :classes-only:
    :toctree: _generated/


This comes after
"""


def test_am_replacer_noinh():
    """
    Tests replacing an ".. automodapi::" with no-inheritance-diagram
    option
    """
    from ..automodapi import automodapi_replace

    fakeapp = FakeApp(automodapi_toctreedirnm='_generated/')
    ops = ['', ':no-inheritance-diagram:']
    ostr = '\n    '.join(ops)
    result = automodapi_replace(am_replacer_str.format(options=ostr), fakeapp)

    assert result == am_replacer_noinh_expected

am_replacer_titleandhdrs_expected = """
This comes before

astropy.sphinx.ext.tests.test_automodapi Module
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

.. automodule:: astropy.sphinx.ext.tests.test_automodapi

Functions
*********

.. automodsumm:: astropy.sphinx.ext.tests.test_automodapi
    :functions-only:
    :toctree: _generated/

Classes
*******

.. automodsumm:: astropy.sphinx.ext.tests.test_automodapi
    :classes-only:
    :toctree: _generated/

Class Inheritance Diagram
*************************

.. automod-diagram:: astropy.sphinx.ext.tests.test_automodapi
    :private-bases:


This comes after
"""


def test_am_replacer_titleandhdrs():
    """
    Tests replacing an ".. automodapi::" entry with title-setting and header
    character options.
    """
    from ..automodapi import automodapi_replace

    fakeapp = FakeApp(automodapi_toctreedirnm='_generated/')
    ops = ['', ':title: A new title', ':headings: &*']
    ostr = '\n    '.join(ops)
    result = automodapi_replace(am_replacer_str.format(options=ostr), fakeapp)

    assert result == am_replacer_titleandhdrs_expected


am_replacer_ss_str = """
This comes before

.. automodapi:: astropy.sphinx.ext
    :subsections: automodsumm,automodapi,tests,falsemod
    :no-main-section:

This comes after
"""

am_replacer_subsections_expected = """
This comes before

astropy.sphinx.ext.automodsumm Module
-------------------------------------

.. automodule:: astropy.sphinx.ext.automodsumm

Functions
^^^^^^^^^

.. automodsumm:: astropy.sphinx.ext.automodsumm
    :functions-only:
    :toctree: _generated/

Classes
^^^^^^^

.. automodsumm:: astropy.sphinx.ext.automodsumm
    :classes-only:
    :toctree: _generated/

Class Inheritance Diagram
^^^^^^^^^^^^^^^^^^^^^^^^^

.. automod-diagram:: astropy.sphinx.ext.automodsumm
    :private-bases:

astropy.sphinx.ext.automodapi Module
------------------------------------

.. automodule:: astropy.sphinx.ext.automodapi

Functions
^^^^^^^^^

.. automodsumm:: astropy.sphinx.ext.automodapi
    :functions-only:
    :toctree: _generated/

astropy.sphinx.ext.tests Module
-------------------------------

.. automodule:: astropy.sphinx.ext.tests


This comes after
"""


def test_am_replacer_subsections():
    """
    Tests replacing an ".. automodapi::" with "subsections" and
    "no-main-section" options.
    """
    from ..automodapi import automodapi_replace

    fakeapp = FakeApp(automodapi_toctreedirnm='_generated/')
    result = automodapi_replace(am_replacer_ss_str, fakeapp)

    expected_warnings = [('Attempted to add documentation section for '
      'astropy.sphinx.ext.falsemod, which is not importable. Skipping.', None)]

    assert fakeapp.warnings == expected_warnings
    assert result == am_replacer_subsections_expected
