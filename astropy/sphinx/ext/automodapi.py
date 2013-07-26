# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This sphinx extension adds a tools to simplify generating the API
documentationfor Astropy packages and affiliated packages.

======================
`automodapi` directive
======================
This directive takes a single argument that must be module or package.
It will produce a block of documentation that includes the docstring for
the package, an `automodsumm` directive, and an `automod-diagram` if
there are any classes in the module.

It accepts the following options:

    * ``:no-inheritance-diagram:``
        If present, the inheritance diagram will not be shown even if
        the module/package has classes.

    * ``:skip: str``
        This option results in the
        specified object being skipped, that is the object will *not* be
        included in the generated documentation. This option may appear
        any number of times to skip multiple objects.

    * ``:no-main-docstr:``
        If present, the docstring for the module/package will not be generated.
        The function and class tables will still be used, however.

    * ``:headings: str``
        Specifies the characters (in one string) used as the heading
        levels used for the generated section. This must have at least 2
        characters (any after 2 will be ignored). This also *must* match
        the rest of the documentation on this page for sphinx to be
        happy. Defaults to "-^", which matches the convention used for
        Python's documentation, assuming the automodapi call is inside a
        top-level section (which usually uses '=').

This extension also adds a sphinx configuration option
`automodapi_toctreedirnm`. It must be a string that specifies the name of
the directory the automodsumm generated documentation ends up in. This
directory path should be relative to the documentation root (e.g., same
place as ``index.rst``). It defaults to '_generated'

"""

# Implementation note:
# The 'automodapi' directive is not actually implemented as a docutils
# directive. Instead, this extension searches for the 'automodapi' text in
# all sphinx documents, and replaces it where necessary from a template built
# into this extension. This is necessary because automodsumm (and autosummary)
# use the "builder-inited" event, which comes before the directives are
# actually built.

import re

automod_templ_modheader = """
{modname} {pkgormod}
{modhds}{pkgormodhds}

{automoduleline}
"""

automod_templ_classes = """
Classes
{clshds}

.. automodsumm:: {modname}
    :classes-only:
    {toctree}
    {skips}
"""

automod_templ_funcs = """
Functions
{funchds}

.. automodsumm:: {modname}
    :functions-only:
    {toctree}
    {skips}
"""

automod_templ_inh = """
Class Inheritance Diagram
{clsinhsechds}

.. automod-diagram:: {modname}
    :private-bases:
"""

_automodapirex = re.compile(r'^(?:\s*\.\.\s+automodapi::\s*)([A-Za-z0-9_.]+)'
                            r'\s*$((?:\n\s+:[a-zA-Z_\-]+:.*$)*)',
                            flags=re.MULTILINE)
#the last group of the above regex is intended to go into finall with the below
_automodapiargsrex = re.compile(r':([a-zA-Z_\-]+):(.*)$', flags=re.MULTILINE)


def automodapi_replace(sourcestr, app, dotoctree=True, docname=None,
                       warnings=True):
    """
    Replaces `sourcestr`'s entries of ".. automdapi::" with the
    automodapi template form based on provided options.

    This is used with the sphinx event 'source-read' to replace
    `automodapi` entries before sphinx actually processes them, as
    automodsumm needs the code to be present to generate stub
    documentation.

    Parameters
    ----------
    sourcestr : str
        The string with sphinx source to be checked for automodapi
        replacement.
    app : `sphinx.application.Application`
        The sphinx application.
    dotoctree : bool
        If True, a ":toctree:" option will be added in the "..
        automodsumm::" sections of the template, pointing to the
        appropriate "generated" directory based on the Astropy convention
        (e.g. in ``docs/_generated``)
    docname : str
        The name of the file for this `sourcestr` (if known - if not, it
        can be None). If not provided and `dotoctree` is True, the
        generated files may end up in the wrong place.
    warnings : bool
        If False, all warnings that would normally be issued are
        silenced.

    Returns
    -------
    newstr :str
        The string with automodapi entries replaced with the correct
        sphinx markup.
    """
    from os import sep

    spl = _automodapirex.split(sourcestr)
    if len(spl) > 1:  # automodsumm is in this document

        if dotoctree:
            toctreestr = ':toctree: '
            dirnm = app.config.automodapi_toctreedirnm
            if not dirnm.endswith(sep):
                dirnm += sep
            if docname is not None:
                toctreestr += '../' * docname.count('/') + dirnm
            else:
                toctreestr += dirnm
        else:
            toctreestr = ''

        newstrs = [spl[0]]
        for grp in range(len(spl) // 3):
            modnm = spl[grp * 3 + 1]

            #find where this is in the document for warnings
            if docname is None:
                location = None
            else:
                location = (docname, spl[0].count('\n'))

            #initialize default options
            toskip = []
            inhdiag = maindocstr = True
            hds = '-^'

            #look for actual options
            unknownops = []
            for opname, args in _automodapiargsrex.findall(spl[grp * 3 + 2]):
                if opname == 'skip':
                    toskip.append(args.strip())
                elif opname == 'no-inheritance-diagram':
                    inhdiag = False
                elif opname == 'no-main-docstr':
                    maindocstr = False
                elif opname == 'headings':
                    hds = args
                else:
                    unknownops.append(opname)

            # get the two heading chars
            if len(hds) < 2:
                msg = 'Not enough headings (got {0}, need 2), using default -^'
                if warnings:
                    app.warn(msg.format(len(hds)), location)
                hds = '-^'
            h1, h2 = hds.lstrip()[:2]

            #tell sphinx that the remaining args are invalid.
            if len(unknownops) > 0 and app is not None:
                opsstrs = ','.join(unknownops)
                msg = 'Found additional options ' + opsstrs + ' in automodapi.'
                if warnings:
                    app.warn(msg, location)

            ispkg, hascls, hasfuncs = _mod_info(modnm, toskip)

            #add automodule directive only if no-main-docstr isn't present
            if maindocstr:
                automodline = '.. automodule:: {modname}'.format(modname=modnm)
            else:
                automodline = ''

            newstrs.append(automod_templ_modheader.format(modname=modnm,
                modhds=h1 * len(modnm),
                pkgormod='Package' if ispkg else 'Module',
                pkgormodhds=h1 * (8 if ispkg else 7),
                automoduleline=automodline))

            if hasfuncs:
                newstrs.append(automod_templ_funcs.format(modname=modnm,
                    funchds=h2 * 9,
                    toctree=toctreestr,
                    skips=':skip: ' + ','.join(toskip) if toskip else ''))

            if hascls:
                newstrs.append(automod_templ_classes.format(modname=modnm,
                    clshds=h2 * 7,
                    toctree=toctreestr,
                    skips=':skip: ' + ','.join(toskip) if toskip else ''))

            if inhdiag and hascls:
                # add inheritance diagram if any classes are in the module
                newstrs.append(automod_templ_inh.format(
                    modname=modnm, clsinhsechds=h2 * 25))

            newstrs.append(spl[grp * 3 + 3])
        return ''.join(newstrs)
    else:
        return sourcestr


def _mod_info(modname, toskip=[]):
    """
    Determines if a module is a module or a package and whether or not
    it has classes or functions.
    """
    import sys

    from os.path import split
    from inspect import isclass, isfunction
    from ...utils.misc import find_mod_objs

    hascls = hasfunc = False

    for localnm, fqnm, obj in zip(*find_mod_objs(modname, onlylocals=True)):
        if localnm not in toskip:
            hascls = hascls or isclass(obj)
            hasfunc = hasfunc or isfunction(obj)
            if hascls and hasfunc:
                break

    #find_mod_objs has already imported modname
    pkg = sys.modules[modname]
    ispkg = '__init__.' in split(pkg.__name__)[1]

    return ispkg, hascls, hasfunc


def process_automodapi(app, docname, source):
    source[0] = automodapi_replace(source[0], app, True, docname)


def setup(app):
    # need automodsumm for automodapi
    app.setup_extension('astropy.sphinx.ext.automodsumm')

    app.connect('source-read', process_automodapi)

    app.add_config_value('automodapi_toctreedirnm', '_generated', True)
