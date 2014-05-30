# **Please Note**: ``astropy.sphinx`` exists only for backward-compatibility
# purposes - it has now been moved to the separate astropy-helpers package,
# located at https://github.com/astropy/astropy-helpers. Any new development or
# bug fixes should be done there.
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This sphinx extension adds two directives for summarizing the public
members of a module or package.

These directives are primarily for use with the `automodapi`_ extension,
but can be used independently.

.. _automodsumm:

=======================
automodsumm directive
=======================

This directive will produce an "autosummary"-style table for public
attributes of a specified module. See the `sphinx.ext.autosummary`_ extension
for details on this process. The main difference from the `autosummary`_
directive is that `autosummary`_ requires manually inputting all attributes
that appear in the table, while this captures the entries automatically.

This directive requires a single argument that must be a module or
package.

It also accepts any options supported by the `autosummary`_ directive-
see `sphinx.ext.autosummary`_ for details. It also accepts two additional
options:

    * ``:classes-only:``
        If present, the autosummary table will only contain entries for
        classes. This cannot be used at the same time with
        ``:functions-only:`` .

    * ``:functions-only:``
        If present, the autosummary table will only contain entries for
        functions. This cannot be used at the same time with
        ``:classes-only:`` .

    * ``:skip: obj1, [obj2, obj3, ...]``
        If present, specifies that the listed objects should be skipped
        and not have their documentation generated, nor be includded in
        the summary table.

    * ``:allowed-package-names: pkgormod1, [pkgormod2, pkgormod3, ...]``
        Specifies the packages that functions/classes documented here are
        allowed to be from, as comma-separated list of package names. If not
        given, only objects that are actually in a subpackage of the package
        currently being documented are included.

This extension also adds one sphinx configuration option:

* ``automodsumm_writereprocessed``
    Should be a bool, and if True, will cause `automodsumm`_ to write files
    with any ``automodsumm`` sections replaced with the content Sphinx
    processes after ``automodsumm`` has run.  The output files are not
    actually used by sphinx, so this option is only for figuring out the
    cause of sphinx warnings or other debugging.  Defaults to `False`.

.. _sphinx.ext.autosummary: http://sphinx-doc.org/latest/ext/autosummary.html
.. _autosummary: http://sphinx-doc.org/latest/ext/autosummary.html#directive-autosummary

.. _automod-diagram:

===========================
automod-diagram directive
===========================

This directive will produce an inheritance diagram like that of the
`sphinx.ext.inheritance_diagram`_ extension.

This directive requires a single argument that must be a module or
package. It accepts no options.

.. note::
    Like 'inheritance-diagram', 'automod-diagram' requires
    `graphviz <http://www.graphviz.org/>`_ to generate the inheritance diagram.

.. _sphinx.ext.inheritance_diagram: http://sphinx-doc.org/latest/ext/inheritance.html
"""

import inspect
import os
import re

from sphinx.ext.autosummary import Autosummary
from sphinx.ext.inheritance_diagram import InheritanceDiagram
from docutils.parsers.rst.directives import flag

from ...utils.misc import find_mod_objs
from .astropyautosummary import AstropyAutosummary


def _str_list_converter(argument):
    """
    A directive option conversion function that converts the option into a list
    of strings. Used for 'skip' option.
    """
    if argument is None:
        return []
    else:
        return [s.strip() for s in argument.split(',')]


class Automodsumm(AstropyAutosummary):
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = False
    has_content = False
    option_spec = dict(Autosummary.option_spec)
    option_spec['functions-only'] = flag
    option_spec['classes-only'] = flag
    option_spec['skip'] = _str_list_converter
    option_spec['allowed-package-names'] = _str_list_converter

    def run(self):
        env = self.state.document.settings.env
        modname = self.arguments[0]

        self.warnings = []
        nodelist = []

        try:
            localnames, fqns, objs = find_mod_objs(modname)
        except ImportError:
            self.warnings = []
            self.warn("Couldn't import module " + modname)
            return self.warnings

        try:
            # set self.content to trick the Autosummary internals.
            # Be sure to respect functions-only and classes-only.
            funconly = 'functions-only' in self.options
            clsonly = 'classes-only' in self.options

            skipnames = []
            if 'skip' in self.options:
                option_skipnames = set(self.options['skip'])
                for lnm in localnames:
                    if lnm in option_skipnames:
                        option_skipnames.remove(lnm)
                        skipnames.append(lnm)
                if len(option_skipnames) > 0:
                    self.warn('Tried to skip objects {objs} in module {mod}, '
                              'but they were not present.  Ignoring.'.format(
                              objs=option_skipnames, mod=modname))

            if funconly and not clsonly:
                cont = []
                for nm, obj in zip(localnames, objs):
                    if nm not in skipnames and inspect.isfunction(obj):
                        cont.append(nm)
            elif clsonly:
                cont = []
                for nm, obj in zip(localnames, objs):
                    if nm not in skipnames and inspect.isclass(obj):
                        cont.append(nm)
            else:
                if clsonly and funconly:
                    self.warning('functions-only and classes-only both '
                                 'defined. Skipping.')
                cont = [nm for nm in localnames if nm not in skipnames]

            self.content = cont

            #for some reason, even though ``currentmodule`` is substituted in, sphinx
            #doesn't necessarily recognize this fact.  So we just force it
            #internally, and that seems to fix things
            env.temp_data['py:module'] = modname

            #can't use super because Sphinx/docutils has trouble
            #return super(Autosummary,self).run()
            nodelist.extend(Autosummary.run(self))
            return self.warnings + nodelist
        finally:  # has_content = False for the Automodsumm
            self.content = []


#<-------------------automod-diagram stuff------------------------------------>
class Automoddiagram(InheritanceDiagram):

    option_spec = dict(InheritanceDiagram.option_spec)
    option_spec['allowed-package-names'] = _str_list_converter

    def run(self):
        try:
            ols = self.options.get('allowed-package-names', [])
            ols = True if len(ols) == 0 else ols  # if none are given, assume only local

            nms, objs = find_mod_objs(self.arguments[0], onlylocals=ols)[1:]
        except ImportError:
            self.warnings = []
            self.warn("Couldn't import module " + self.arguments[0])
            return self.warnings

        clsnms = []
        for n, o in zip(nms, objs):

            if inspect.isclass(o):
                clsnms.append(n)

        oldargs = self.arguments
        try:
            if len(clsnms) > 0:
                self.arguments = [u' '.join(clsnms)]
            return InheritanceDiagram.run(self)
        finally:
            self.arguments = oldargs


#<---------------------automodsumm generation stuff--------------------------->
def process_automodsumm_generation(app):
    env = app.builder.env
    ext = app.config.source_suffix

    filestosearch = [x + ext for x in env.found_docs
                     if os.path.isfile(env.doc2path(x))]\

    liness = []
    for sfn in filestosearch:
        lines = automodsumm_to_autosummary_lines(sfn, app)
        liness.append(lines)
        if app.config.automodsumm_writereprocessed:
            if lines:  # empty list means no automodsumm entry is in the file
                outfn = os.path.join(app.srcdir, sfn) + '.automodsumm'
                with open(outfn, 'w') as f:
                    for l in lines:
                        f.write(l)
                        f.write('\n')

    for sfn, lines in zip(filestosearch, liness):
        if len(lines) > 0:
            generate_automodsumm_docs(lines, sfn, builder=app.builder,
                                      warn=app.warn, info=app.info,
                                      suffix=app.config.source_suffix,
                                      base_path=app.srcdir)

#_automodsummrex = re.compile(r'^(\s*)\.\. automodsumm::\s*([A-Za-z0-9_.]+)\s*'
#                             r'\n\1(\s*)(\S|$)', re.MULTILINE)
_lineendrex = r'(?:\n|$)'
_hdrex = r'^\n?(\s*)\.\. automodsumm::\s*(\S+)\s*' + _lineendrex
_oprex1 = r'(?:\1(\s+)\S.*' + _lineendrex + ')'
_oprex2 = r'(?:\1\4\S.*' + _lineendrex + ')'
_automodsummrex = re.compile(_hdrex + '(' + _oprex1 + '?' + _oprex2 + '*)',
                             re.MULTILINE)


def automodsumm_to_autosummary_lines(fn, app):
    """
    Generates lines from a file with an "automodsumm" entry suitable for
    feeding into "autosummary".

    Searches the provided file for `automodsumm` directives and returns
    a list of lines specifying the `autosummary` commands for the modules
    requested. This does *not* return the whole file contents - just an
    autosummary section in place of any :automodsumm: entries. Note that
    any options given for `automodsumm` are also included in the
    generated `autosummary` section.

    Parameters
    ----------
    fn : str
        The name of the file to search for `automodsumm` entries.
    app : sphinx.application.Application
        The sphinx Application object

    Return
    ------
    lines : list of str
        Lines for all `automodsumm` entries with the entries replaced by
        `autosummary` and the module's members added.


    """
    fullfn = os.path.join(app.builder.env.srcdir, fn)

    with open(fullfn) as fr:
        if 'astropy.sphinx.ext.automodapi' in app._extensions:
            from astropy.sphinx.ext.automodapi import automodapi_replace
            # Must do the automodapi on the source to get the automodsumm
            # that might be in there
            filestr = automodapi_replace(fr.read(), app, True, fn, False)
        else:
            filestr = fr.read()

    spl = _automodsummrex.split(filestr)
    #0th entry is the stuff before the first automodsumm line
    indent1s = spl[1::5]
    mods = spl[2::5]
    opssecs = spl[3::5]
    indent2s = spl[4::5]
    remainders = spl[5::5]

    # only grab automodsumm sections and convert them to autosummary with the
    # entries for all the public objects
    newlines = []

    #loop over all automodsumms in this document
    for i, (i1, i2, modnm, ops, rem) in enumerate(zip(indent1s, indent2s, mods,
                                                    opssecs, remainders)):
        allindent = i1 + ('' if i2 is None else i2)

        #filter out functions-only and classes-only options if present
        oplines = ops.split('\n')
        toskip = []
        allowedpkgnms = []
        funcsonly = clssonly = False
        for i, ln in reversed(list(enumerate(oplines))):
            if ':functions-only:' in ln:
                funcsonly = True
                del oplines[i]
            if ':classes-only:' in ln:
                clssonly = True
                del oplines[i]
            if ':skip:' in ln:
                toskip.extend(_str_list_converter(ln.replace(':skip:', '')))
                del oplines[i]
            if ':allowed-package-names:' in ln:
                allowedpkgnms.extend(_str_list_converter(ln.replace(':allowed-package-names:', '')))
                del oplines[i]
        if funcsonly and clssonly:
            msg = ('Defined both functions-only and classes-only options. '
                   'Skipping this directive.')
            lnnum = sum([spl[j].count('\n') for j in range(i * 5 + 1)])
            app.warn('[automodsumm]' + msg, (fn, lnnum))
            continue

        # Use the currentmodule directive so we can just put the local names
        # in the autosummary table.  Note that this doesn't always seem to
        # actually "take" in Sphinx's eyes, so in `Automodsumm.run`, we have to
        # force it internally, as well.
        newlines.extend([i1 + '.. currentmodule:: ' + modnm,
                         '',
                         '.. autosummary::'])
        newlines.extend(oplines)

        ols = True if len(allowedpkgnms) == 0 else allowedpkgnms
        for nm, fqn, obj in zip(*find_mod_objs(modnm, onlylocals=ols)):
            if nm in toskip:
                continue
            if funcsonly and not inspect.isfunction(obj):
                continue
            if clssonly and not inspect.isclass(obj):
                continue
            newlines.append(allindent + nm)

    return newlines


def generate_automodsumm_docs(lines, srcfn, suffix='.rst', warn=None,
                              info=None, base_path=None, builder=None,
                              template_dir=None):
    """
    This function is adapted from
    `sphinx.ext.autosummary.generate.generate_autosummmary_docs` to
    generate source for the automodsumm directives that should be
    autosummarized. Unlike generate_autosummary_docs, this function is
    called one file at a time.
    """

    from sphinx.jinja2glue import BuiltinTemplateLoader
    from sphinx.ext.autosummary import import_by_name, get_documenter
    from sphinx.ext.autosummary.generate import (find_autosummary_in_lines,
                                                 _simple_info, _simple_warn)
    from sphinx.util.osutil import ensuredir
    from sphinx.util.inspect import safe_getattr
    from jinja2 import FileSystemLoader, TemplateNotFound
    from jinja2.sandbox import SandboxedEnvironment

    if info is None:
        info = _simple_info
    if warn is None:
        warn = _simple_warn

    #info('[automodsumm] generating automodsumm for: ' + srcfn)

    # Create our own templating environment - here we use Astropy's
    # templates rather than the default autosummary templates, in order to
    # allow docstrings to be shown for methods.
    template_dirs = [os.path.join(os.path.dirname(__file__), 'templates'),
                     os.path.join(base_path, '_templates')]
    if builder is not None:
        # allow the user to override the templates
        template_loader = BuiltinTemplateLoader()
        template_loader.init(builder, dirs=template_dirs)
    else:
        if template_dir:
            template_dirs.insert(0, template_dir)
        template_loader = FileSystemLoader(template_dirs)
    template_env = SandboxedEnvironment(loader=template_loader)

    # read
    #items = find_autosummary_in_files(sources)
    items = find_autosummary_in_lines(lines, filename=srcfn)
    if len(items) > 0:
        msg = '[automodsumm] {1}: found {0} automodsumm entries to generate'
        info(msg.format(len(items), srcfn))

#    gennms = [item[0] for item in items]
#    if len(gennms) > 20:
#        gennms = gennms[:10] + ['...'] + gennms[-10:]
#    info('[automodsumm] generating autosummary for: ' + ', '.join(gennms))

    # remove possible duplicates
    items = dict([(item, True) for item in items]).keys()

    # keep track of new files
    new_files = []

    # write
    for name, path, template_name in sorted(items):
        if path is None:
            # The corresponding autosummary:: directive did not have
            # a :toctree: option
            continue

        path = os.path.abspath(path)
        ensuredir(path)

        try:
            name, obj, parent = import_by_name(name)
        except ImportError, e:
            warn('[automodsumm] failed to import %r: %s' % (name, e))
            continue

        fn = os.path.join(path, name + suffix)

        # skip it if it exists
        if os.path.isfile(fn):
            continue

        new_files.append(fn)

        f = open(fn, 'w')

        try:
            doc = get_documenter(obj, parent)

            if template_name is not None:
                template = template_env.get_template(template_name)
            else:
                tmplstr = 'autosummary/%s.rst'
                try:
                    template = template_env.get_template(tmplstr % doc.objtype)
                except TemplateNotFound:
                    template = template_env.get_template(tmplstr % 'base')

            def get_members_mod(obj, typ, include_public=[]):
                """
                typ = None -> all
                """
                items = []
                for name in dir(obj):
                    try:
                        documenter = get_documenter(safe_getattr(obj, name),
                                                    obj)
                    except AttributeError:
                        continue
                    if typ is None or documenter.objtype == typ:
                        items.append(name)
                public = [x for x in items
                          if x in include_public or not x.startswith('_')]
                return public, items

            def get_members_class(obj, typ, include_public=[],
                                  include_base=False):
                """
                typ = None -> all
                include_base -> include attrs that are from a base class
                """
                items = []

                # using dir gets all of the attributes, including the elements
                # from the base class, otherwise use __slots__ or __dict__
                if include_base:
                    names = dir(obj)
                else:
                    if hasattr(obj, '__slots__'):
                        names = tuple(getattr(obj, '__slots__'))
                    else:
                        names = getattr(obj, '__dict__').keys()

                for name in names:
                    try:
                        documenter = get_documenter(safe_getattr(obj, name),
                                                    obj)
                    except AttributeError:
                        continue
                    if typ is None or documenter.objtype == typ:
                        items.append(name)
                public = [x for x in items
                          if x in include_public or not x.startswith('_')]
                return public, items

            ns = {}

            if doc.objtype == 'module':
                ns['members'] = get_members_mod(obj, None)
                ns['functions'], ns['all_functions'] = \
                                   get_members_mod(obj, 'function')
                ns['classes'], ns['all_classes'] = \
                                 get_members_mod(obj, 'class')
                ns['exceptions'], ns['all_exceptions'] = \
                                   get_members_mod(obj, 'exception')
            elif doc.objtype == 'class':
                api_class_methods = ['__init__', '__call__']
                ns['members'] = get_members_class(obj, None)
                ns['methods'], ns['all_methods'] = \
                                 get_members_class(obj, 'method', api_class_methods)
                ns['attributes'], ns['all_attributes'] = \
                                 get_members_class(obj, 'attribute')
                ns['methods'].sort()
                ns['attributes'].sort()

            parts = name.split('.')
            if doc.objtype in ('method', 'attribute'):
                mod_name = '.'.join(parts[:-2])
                cls_name = parts[-2]
                obj_name = '.'.join(parts[-2:])
                ns['class'] = cls_name
            else:
                mod_name, obj_name = '.'.join(parts[:-1]), parts[-1]

            ns['fullname'] = name
            ns['module'] = mod_name
            ns['objname'] = obj_name
            ns['name'] = parts[-1]

            ns['objtype'] = doc.objtype
            ns['underline'] = len(name) * '='

            # We now check whether a file for reference footnotes exists for
            # the module being documented. We first check if the
            # current module is a file or a directory, as this will give a
            # different path for the reference file. For example, if
            # documenting astropy.wcs then the reference file is at
            # ../wcs/references.txt, while if we are documenting
            # astropy.config.logging_helper (which is at
            # astropy/config/logging_helper.py) then the reference file is set
            # to ../config/references.txt
            if '.' in mod_name:
                mod_name_dir = mod_name.replace('.', '/').split('/', 1)[1]
            else:
                mod_name_dir = mod_name
            if not os.path.isdir(os.path.join(base_path, mod_name_dir)) \
               and os.path.isdir(os.path.join(base_path, mod_name_dir.rsplit('/', 1)[0])):
                mod_name_dir = mod_name_dir.rsplit('/', 1)[0]

            # We then have to check whether it exists, and if so, we pass it
            # to the template.
            if os.path.exists(os.path.join(base_path, mod_name_dir, 'references.txt')):
                # An important subtlety here is that the path we pass in has
                # to be relative to the file being generated, so we have to
                # figure out the right number of '..'s
                ndirsback = path.replace(base_path, '').count('/')
                ref_file_rel_segments = ['..'] * ndirsback
                ref_file_rel_segments.append(mod_name_dir)
                ref_file_rel_segments.append('references.txt')
                ns['referencefile'] = os.path.join(*ref_file_rel_segments)

            rendered = template.render(**ns)
            f.write(rendered)
        finally:
            f.close()


def setup(app):
    # need our autosummary
    app.setup_extension('astropy.sphinx.ext.astropyautosummary')
    # need inheritance-diagram for automod-diagram
    app.setup_extension('sphinx.ext.inheritance_diagram')

    app.add_directive('automod-diagram', Automoddiagram)
    app.add_directive('automodsumm', Automodsumm)
    app.connect('builder-inited', process_automodsumm_generation)

    app.add_config_value('automodsumm_writereprocessed', False, True)
