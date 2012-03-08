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

    * ``:subsections: mod1[,mod2,subpkg3]``
        If present, this generates separate documentation sections for the
        requested submodules or subpackages.

    * ``:no-main-section:``
        If present, the documentation and summary table for the main module or
        package will not be generated (this would generally only be used with
        ``:subsections:`` to document a set of subsections only.)

    * ``:title: [str]``
        Specifies the top-level title for the section. Defaults to
        "Reference/API".

    * ``:headings: [str]``
        Specifies the characters (all in one string) to use for the heading
        levels.  This *must* have at least 3 characters (any after 3 will be
        ignored).  Defaults to "-^_".  Note that this must match the rest of
        the documentation page.


This extension also adds a sphinx configuration option
`automodapi_toctreedirnm`. It must be a string that specifies the name of the
directory the automodsumm generated documentation ends up in.  This directory
path should be relative to the documentation root (e.g., same place as
``index.rst``).  It defaults to '_generated'

"""

# Implementation note:
# The 'automodapi' directive is not actually implemented as a docutils
# directive. Instead, this extension searches for the 'automodapi' text in
# all sphinx documents, and replaces it where necessary from a template built
# into this extension. This is necessary because automodsumm (and autosummary)
# use the "builder-inited" event, which comes before the directives are
# actually built.

import re

automod_templ_header = """
{title}
{titlehd}
"""

automod_templ_docs = """
{modname} {pkgormod}
{modhdr}{pkgormodhds}

.. automodule:: {modname}

{classesandfunctions}
{classesandfunctionsudrsc}

.. automodsumm:: {modname}
    {toctree}
"""

automod_inh_templ = """
Class Inheritance Diagram
{clsinhsechdr}

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
    Replaces `sourcestr`'s entries of ".. automdapi::" with the automodapi
    template form based on provided options.

    This is used with the sphinx event 'source-read' to replace `automodapi`
    entries before sphinx actually processes them, as automodsumm needs the
    code to be present to generate stub documentation.

    Parameters
    ----------
    sourcestr : str
        The string with sphinx source to be checked for automodapi replacement.
    app : `sphinx.application.Application`
        The sphinx application.
    dotoctree : bool
        If True, a ":toctree:" option will be added in the ".. automodsumm::"
        sections of the template, pointing to the appropriate "generated"
        directory based on the Astropy convention (e.g. in ``docs/_generated``)
    docname : str
        The name of the file for this `sourcestr` (if known - if not, it can be
        None). If not provided and `dotoctree` is True, the generated files may
        end up in the wrong place.
    warnings : bool
        If False, all warnings that would normally be issued are silenced.

    Returns
    -------
    newstr :str
        The string with automodapi entries replaced with the correct sphinx
        markup.
    """
    from os import sep
    from inspect import ismodule

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
            basemodnm = spl[grp * 3 + 1]

            #find where this is in the document for warnings
            if docname is None:
                location = None
            else:
                location = (docname, spl[0].count('\n'))

            #findall yields an optionname, arguments tuple
            modops = dict(_automodapiargsrex.findall(spl[grp * 3 + 2]))

            inhdiag = 'no-inheritance-diagram' not in modops
            modops.pop('no-inheritance-diagram', None)
            subsecs = modops.pop('show-subsections', None)
            nomain = 'no-main-section' in modops
            modops.pop('no-main-section', None)
            sectitle = modops.pop('sectitle', 'Reference/API')
            hds = modops.pop('headings', '-^_')

            if len(hds) < 3:
                msg = 'not enough headings (got {0}, need 3), using default -^_'
                if warnings:
                    app.warn(msg.format(len(hds)), location)
                hds = '-^_'
            h1, h2, h3 = hds[:3]

            #tell sphinx that the remaining args are invalid.
            if len(modops) > 0 and app is not None:
                opsstrs = ','.join(modops.keys())
                msg = 'Found additional options ' + opsstrs + ' in automodapi.'

                if warnings:
                    app.warn(msg, location)

            #now actually populate the templates
            newstrs.append(automod_templ_header.format(title=sectitle,
                titlehd=h1 * len(sectitle)))

            # construct the list of modules to document based on the
            # show-subsections argument
            modnames = [] if nomain else [basemodnm]
            if subsecs is not None:
                for ss in subsecs.replace(' ', '').split(','):
                    submodnm = basemodnm + '.' + ss
                    try:
                        mod = __import__(submodnm)
                        if ismodule(mod):
                            modnames.append(mod.__name__)
                        else:
                            msg = 'Attempted to add documentation section for '
                            '{0}, which is neither module nor package. '
                            'Skipping.'
                            if warnings:
                                app.warn(msg.format(submodnm), location)
                    except ImportError:
                        msg = 'Attempted to add documentation section for '
                        '{0}, which is not importable. Skipping.'
                        if warnings:
                            app.warn(msg.format(submodnm), location)

            for modnm in modnames:
                ispkg, hascls, hasfunc = _mod_info(modnm)

                if hascls and not hasfunc:
                    clsfuncstr = 'Classes'
                if not hascls and hasfunc:
                    clsfuncstr = 'Functions'
                else:
                    clsfuncstr = 'Classes and Functions'

                newstrs.append(automod_templ_docs.format(modname=modnm,
                               modhds=h2 * len(modnm),
                               pkgormod='Package' if ispkg else 'Module',
                               pkgormodhds=h2 * (8 if ispkg else 7),
                               classesandfunctions=clsfuncstr,
                               classesandfunctionsudrsc=h3 * len(clsfuncstr),
                               toctree=toctreestr))

                if inhdiag and hascls:
                    # add inheritance diagram if any classes are in the module
                    newstrs.append(automod_inh_templ.format(
                        modname=modnm, clsinhsechd=h3 * 25))

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
    source[0] = automodapi_replace(source[0], app, True, docname)


def setup(app):
    # need automodsumm for automodapi
    app.setup_extension('astropy.sphinx.ext.automodsumm')

    app.connect('source-read', process_automodapi)

    app.add_config_value('automodapi_toctreedirnm', '_generated', True)
