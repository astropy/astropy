# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This sphinx extension adds a tools to simplify generating the API
documentationfor Astropy packages and affiliated packages.

======================
`automodapi` directive
======================
This directive takes a single argument that must be module or package.
Itwill produce a Documentation section named "Reference/API" that
includes the docstring for the package, an `automodsumm` directive, and
an `automod-diagram` if there are any classes in the module.

It accepts the following options:

    * ``:no-inheritance-diagram:``
        If present, the inheritance diagram will not be shown even if
        themodule/packagehashas classes.

    * ``:show-subsections:``
        If present and the specified package has modules or subpackages,
        the summary table will be split into multiple tables for each
        subpackage.This does not nest - only one layer of
        subpackages/moduleswill be shown.

    * ``:title: [str]``
        Specifies the top-level title for the section. Defaults to
        "Reference/API".

"""

# Implementation note:
# The 'automodapi' directive is not actually implemented as a docutils
# directive. Instead, this extension searches for the 'automodapi' text in
# all sphinx documents, and replaces it where necessary from a template built
# into this extension. This is necessary because automodsumm (and autosummary)
# use the "builder-inited" event, which comes before the directives are
# actually built.

import re

toctreedirnm = '_generated/'

automod_templ_header = """
{title}
{titledashes}
"""

automod_templ_docs = """
{modname} {pkgormod}
{modcrts}^{pkgormodcrts}

.. automodule:: {modname}

{classesandfunctions}
{classesandfunctionsudrsc}

.. automodsumm:: {modname}
    {toctree}
"""

automod_inh_templ = """
Class Inheritance Diagram
_________________________

.. automod-diagram:: {modname}
    :private-bases:
"""

_automodapirex = re.compile(r'^(?:\s*\.\.\s+automodapi::\s*)([A-Za-z0-9_.]+)'
                            r'\s*$((?:\n\s+:[a-zA-Z_\-]+:.*$)*)',
                            flags=re.MULTILINE)
#the last group of the above regex is intended to go into finall with the below
_automodapiargsrex = re.compile(r':([a-zA-Z_\-]+):(.*)$', flags=re.MULTILINE)


def automodapi_replace(sourcestr, dotoctree=True, docname=None, app=None):
    """
    replaces sourcestr's entries of automodapi with automodsumm entries
    if docname is None, _generated may not be in the right place
    if app is None, warnings will pass silently
    """
    from inspect import ismodule
    from .automodsumm import find_mod_objs

    spl = _automodapirex.split(sourcestr)
    if len(spl) > 1:  # automodsumm is in this document

        if dotoctree:
            toctreestr = ':toctree: '
            if docname is not None:
                toctreestr += '../' * docname.count('/') + toctreedirnm
            else:
                toctreestr += toctreedirnm
        else:
            toctreestr = ''

        newstrs = [spl[0]]
        for grp in range(len(spl) // 3):
            basemodnm = spl[grp * 3 + 1]

            #findall yields an optionname, arguments tuple
            modops = dict(_automodapiargsrex.findall(spl[grp * 3 + 2]))

            inhdiag = 'no-inheritance-diagram' not in modops
            modops.pop('no-inheritance-diagram', None)
            subsecs = 'show-subsections' in modops
            modops.pop('show-subsections', None)
            sectitle = modops.pop('sectitle', 'Reference/API')

            #tell sphinx that the remaining args are invalid.
            if len(modops) > 0 and app is not None:
                opsstrs = ','.join(modops.keys())
                msg = 'Found additional options ' + opsstrs + ' in automodapi.'
                if docname is None:
                    location = None
                else:
                    location = (docname, spl[0].count('\n'))
                app.warn(msg, location)

            #now actually populate the templates
            newstrs.append(automod_templ_header.format(title=sectitle,
                titledashes='-' * len(sectitle)))

            # construct the list of modules to document based on the
            # show-subsections argument
            modnames = [basemodnm]
            if subsecs:
                raise NotImplementedError('subsecs not yet supported')
                for obj in find_mod_objs(basemodnm, False):
                    if ismodule(obj):
                        modnames.append(obj.__name__)

            for modnm in modnames:
                ispkg, hascls, hasfunc = _mod_info(modnm)

                if hascls and not hasfunc:
                    clsfuncstr = 'Classes'
                if not hascls and hasfunc:
                    clsfuncstr = 'Functions'
                else:
                    clsfuncstr = 'Classes and Functions'

                newstrs.append(automod_templ_docs.format(modname=modnm,
                               modcrts='^' * len(modnm),
                               pkgormod='Package' if ispkg else 'Module',
                               pkgormodcrts='^' * (7 if ispkg else 6),
                               classesandfunctions=clsfuncstr,
                               classesandfunctionsudrsc='_' * len(clsfuncstr),
                               toctree=toctreestr))

                if inhdiag and hascls:
                    # add inheritance diagram if any classes are in the module
                    newstrs.append(automod_inh_templ.format(modname=modnm))

            newstrs.append(spl[grp * 3 + 3])
        return ''.join(newstrs)
    else:
        return sourcestr


def _mod_info(modname):
    """
    Determines if a module is a module or a package and whether or not it has
    classes or functions.
    """
    import sys

    from os.path import split
    from inspect import isclass, isfunction
    from .automodsumm import find_mod_objs

    hascls = hasfunc = False
    for obj in find_mod_objs(modname, False):
        hascls = hascls or isclass(obj)
        hasfunc = hasfunc or isfunction(obj)
        if hascls and hasfunc:
            break

    #find_mod_objs has already imported modname
    pkg = sys.modules[modname]
    ispkg = '__init__.' in split(pkg.__name__)[1]

    return ispkg, hascls, hasfunc


def process_automodapi(app, docname, source):
    source[0] = automodapi_replace(source[0],
                                   app.config.automodsumm_generate,
                                   docname)


def setup(app):
    # need automodsumm for automodapi
    app.setup_extension('astropy.sphinx.ext.automodsumm')

    app.connect('source-read', process_automodapi)
