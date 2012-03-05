# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This sphinx extension adds a tools to simplify generating the API
documentation for Astropy packages and affiliated packages.

Specifically, this extension ads an 'automodapi' directive.  The
directive takes a single argument that must be module or package. It
will produce a Documentation section named "Reference/API" that includes
the docstring for the package, an 'automodsumm' directive, and an
'automod-diagram' if there are any classes in the module.

"""

# Implementation note:
# The 'automodapi' directive is not actually implemented as a docutils
# directive. Instead, this extension searches for the 'automodapi' text in
# all sphinx documents, and replaces it where necessary from a template built
# into this extension. This is necessary because automodsumm (and autosummary)
# use the "builder-inited" event, which comes before the directives are
# actually built.

import re

from .automodsumm import find_mod_objs

automod_templ = """
Reference/API
-------------

{modname} Package
{modcrts}^^^^^^^^

.. automodule:: {modname}

Classes and Functions
^^^^^^^^^^^^^^^^^^^^^

.. automodsumm:: {modname}
    {toctree}

"""
toctreestr = ':toctree: _generated/'
automod_inh_templ = """
Class Inheritance Diagram
^^^^^^^^^^^^^^^^^^^^^^^^^

.. automod-diagram:: {modname}
    :private-bases:

"""

_automodapirex = re.compile(
    r'^(\s*\.\.\s+automodapi::\s*)([A-Za-z0-9_.]+)\s*$', flags=re.MULTILINE)


def automodapi_replace(sourcestr, dotoctree=True):
    """
    replaces sourcestr's entries of automodapi with automodsumm entries
    """
    from inspect import isclass

    spl = _automodapirex.split(sourcestr)
    if len(spl) > 1:  # automodsumm is in this document
        for grp in reversed(range(len(spl) // 3)):
            modnm = spl[grp * 3 + 2]
            del spl[grp * 3 + 2]

            templstr = automod_templ.format(
                modname=modnm, modcrts='^' * len(modnm),
                toctree=toctreestr if dotoctree else '')

            # add in the inheritance diagram if any classes are in the module
            if any([isclass(obj) for obj in find_mod_objs(modnm, False)]):
                templinhstr = automod_inh_templ.format(modname=modnm)
                spl[grp * 3 + 1] = templstr + templinhstr
            else:
                spl[grp * 3 + 1] = templstr
    return ''.join(spl)


def process_automodapi(app, docname, source):
    source[0] = automodapi_replace(source[0], app.config.automodsumm_generate)


def setup(app):
    # need automodsumm for automodapi
    app.setup_extension('astropy.sphinx.ext.automodsumm')

    app.connect('source-read', process_automodapi)
