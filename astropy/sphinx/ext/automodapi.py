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
toctreedirnm = '_generated/'
automod_inh_templ = """
Class Inheritance Diagram
^^^^^^^^^^^^^^^^^^^^^^^^^

.. automod-diagram:: {modname}
    :private-bases:

"""

_automodapirex = re.compile(r'^(?:\s*\.\.\s+automodapi::\s*)([A-Za-z0-9_.]+)'
                            r'\s*$((?:\n\s+:[a-zA-Z_\-]+:.*$)*)',
                            flags=re.MULTILINE)
#the last group of the above regex is intended to go into finall with the below
_automodapiargsrex = re.compile(r':([a-zA-Z_\-]+):(.*)$', flags=re.MULTILINE)

def automodapi_replace(sourcestr, dotoctree=True):
    """
    replaces sourcestr's entries of automodapi with automodsumm entries
    """
    from inspect import isclass

    spl = _automodapirex.split(sourcestr)
    if len(spl) > 1:  # automodsumm is in this document

        if dotoctree:
            toctreestr = ':toctree: '
            if isinstance(dotoctree, basestring):
                toctreestr += '../' * dotoctree.count('/') + toctreedirnm
            else:
                toctreestr += toctreedirnm
        else:
            toctreestr = ''

        newstrs = [spl[0]]
        for grp in range(len(spl) // 3):
            modnm = spl[grp * 3 + 1]
            #each of modops is a optionname, arguments tuple
            modops = _automodapiargsrex.findall(spl[grp * 3 + 2])

            #do something with modops

            newstrs.append(automod_templ.format(
                modname=modnm, modcrts='^' * len(modnm),
                toctree=toctreestr))

            # add in the inheritance diagram if any classes are in the module
            if any([isclass(obj) for obj in find_mod_objs(modnm, False)]):
                templinhstr = automod_inh_templ.format(modname=modnm)
                newstrs[-1] += templinhstr
            newstrs.append(spl[grp * 3 + 3])
        return ''.join(newstrs)
    else:
        return sourcestr


def process_automodapi(app, docname, source):
    dotoctree = docname if app.config.automodsumm_generate else False
    source[0] = automodapi_replace(source[0], dotoctree)


def setup(app):
    # need automodsumm for automodapi
    app.setup_extension('astropy.sphinx.ext.automodsumm')

    app.connect('source-read', process_automodapi)
