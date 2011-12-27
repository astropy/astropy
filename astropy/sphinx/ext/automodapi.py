# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This sphinx extensions adds tools to simplify generating the API documentation
for Astropy packages.

It adds three components:

* An 'automodsumm' directive that requires a single argument that must be a
  module or package.  It will then produce an "autosummary" table for public
  attributes of that module, as produced by the `sphinx.ext.autosummary`
  extension.  The main difference from the autosummary directive is that
  autosummary requires manually inputting all attributes that appear in the
  table, while this captures the entries automatically.  Otherwise, it has the
  same options as the 'autosummary' directive.

* An 'automod-diagram' directive that requires a single argument that must be a
  module or package.  The will produce an inheritance diagram like that of the
  `sphinx.ext.inheritance_diagram` extension.

.. note::
  Like 'inheritance-diagram', 'automod-diagram' requires graphviz to generate
  the inheritance diagram.

* An 'automodapi' directive that requires a single argument that must be a
  module or package. This will produce an 'Reference/API' section following the
  standard Astropy scheme, including the docstring for the package, an
  'automodsumm' directive, and an 'automod-diagram' if there are any classes in
  the module.

.. note::
    The 'automodapi' directive is, strictly speaking, not a formal docutils
    directive.  Instead, this extension searches for the 'automodapi' text in
    all sphinx documents, and replaces it where necessary by the appropriate
    documentation sections.  This probably doesn't matter unless you want to
    extend this extension somehow, but be aware that this could cause oddities
    if you ever plan to change its behavior.

To control autogeneration of stub documentation for automodsumm directives,
this extension adds a 'automodsumm_generate' configuration option. If it is
True, 'automodsumm' (and 'automodapi') directives where the ':toctree:' option
is set cause stub documentation to be produced for the items in the
summary table. For 'automodapi', the generated files are placed in a directory
'_generated', while for 'automodsumm', the 'toctree' option can be used to
specify where the documentation goes. See the `sphinx.ext.autosummary`
documentation for more details on controlling the output of this generation.

"""
import re

from sphinx.ext.autosummary import Autosummary
from sphinx.ext.inheritance_diagram import InheritanceDiagram


class Automodsumm(Autosummary):
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = False
    has_content = False
    #option_spec = dict(Autosummary.option_spec)

    def run(self):
        try:
            modnms = find_mod_objs(self.arguments[0], True)
        except ImportError:
            self.warnings = []
            self.warn("Couldn't import module " + self.arguments[0])
            return self.warnings

        try:
            #set self.content to trick the Autosummary internals
            self.content = ['~' + objname for objname in modnms]
            #can't use super because Sphinx/docutils has trouble
            #return super(Autosummary,self).run()
            return Autosummary.run(self)
        finally:  # has_content = False for the Automodsumm
            self.content = []


def find_mod_objs(modname, names=False):
    """ Find the all public attributes of a module or package.

    Parameters
    ----------
    modname : str
        The name of the module to search.
    names : bool
        If True, the attribute's names will be returned, otherwise the objects
        themselves.

    Returns
    -------
    objs : list
        A list with the fully-qualified names of the module's public attributes
        If `names` is True, or the objects themselves if it is False.

    """
    import sys

    __import__(modname)
    pkg = sys.modules[modname]
    pkgd = dict([(k, v) for k, v in pkg.__dict__.iteritems() if k[0] != '_'])

    #full-qualified names
    fullnames = []
    for n in pkgd:
        obj = pkgd[n]
        if hasattr(obj, '__module__'):
            fullnames.append(obj.__module__ + '.' + obj.__name__)
        else:
            fullnames.append(None)

    #sort on fqn
    sortnms = sorted(zip(fullnames, pkgd))
    fullnames = [nmtuple[0] for nmtuple in sortnms if nmtuple[0] is not None]
    attrnames = [nmtuple[1] for nmtuple in sortnms if nmtuple[0] is not None]

    if names:
        return fullnames
    else:
        return [pkgd[n] for n in attrnames]


#<-------------------automod-diagram stuff------------------------------------->
class Automoddiagram(InheritanceDiagram):
    def run(self):
        from inspect import isclass

        try:
            nms = find_mod_objs(self.arguments[0], True)
            objs = find_mod_objs(self.arguments[0], False)
        except ImportError:
            self.warnings = []
            self.warn("Couldn't import module " + self.arguments[0])
            return self.warnings

        clsnms = []
        for n, o in zip(nms, objs):

            if isclass(o):
                clsnms.append(n)

        oldargs = self.arguments
        try:
            if len(clsnms) > 0:
                self.arguments = [u' '.join(clsnms)]
            return InheritanceDiagram.run(self)
        finally:
            self.arguments = oldargs

#<------------------automodapi "directive" stuff------------------------------->
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


#<---------------------automodsumm generation stuff--------------------------->


def process_automodsumm_generation(app):
    import os

    if app.config.automodsumm_generate:
        env = app.builder.env
        ext = app.config.source_suffix

        filestosearch = [x + ext for x in env.found_docs
                         if os.path.isfile(env.doc2path(x))]\

        liness = []
        for sfn in filestosearch:
            fullfn = os.path.join(env.srcdir, sfn)
            liness.append(automodsumm_to_autosummary_lines(fullfn, app))

        for sfn, lines in zip(filestosearch, liness):
            if len(lines) > 0:
                generate_automodsumm_docs(lines, sfn, builder=app.builder,
                                          warn=app.warn, info=app.info,
                                          suffix=app.config.source_suffix,
                                          base_path=app.srcdir)

_automodsummrex = re.compile(
    r'^(\s*)\.\.\s+(automodsumm)::\s*([A-Za-z0-9_.]+)\s*$')


def automodsumm_to_autosummary_lines(fn, app):
    """
    Searches the provided file for automodsumm directives and returns a list of
    lines where they've been replaced by appropriate autosummary directives
    """
    config = app.config

    lines = []
    singleindent = ' ' * 4

    with open(fn) as fr:
        # Must do the automodapi on the sourcex to get the automodsumm
        # that might be in there
        filestr = automodapi_replace(fr.read(), config.automodsumm_generate)

    indentrex = None
    modnm = None
    indent = ''

    def finish_content(lines, modnm, indent):
        """Adds the content items for autosummary to the lines"""
        lines.append(indent + singleindent)
        for objnm in find_mod_objs(modnm, names=True):
            lines.append(indent + singleindent + '~' + objnm)
        lines.append(indent + singleindent)

    # iterate over lines in the source file, searching for automodsumm
    for l in filestr.split('\n'):
        if indentrex is not None:  # inside automodsumm
            if not indentrex.match(l):
                finish_content(lines, modnm, indent)
                indentrex = None
                modnm = None
                indent = ''
            else:
                lines.append(l)
        else:
            m = _automodsummrex.match(l)
            if m:
                indent = m.group(1)
                modnm = m.group(3)

                amsstart, amsend = m.regs[2]
                mnmstart, mnmend = m.regs[3]
                lines.append((l[:amsstart] + 'autosummary' +
                              l[amsend:mnmstart] + l[mnmend:]))

                #now make the regex for checking if still in automodsumm
                rexstr = '({0})|({1})'.format(r'\s*$',
                                              r'^' + indent + r'\s+')
                indentrex = re.compile(rexstr)
    if indentrex is not None:  # the file ends while still in automodsumm
        finish_content(lines, modnm, indent)

    return lines


def generate_automodsumm_docs(lines, srcfn, suffix='.rst', warn=None,
                              info=None, base_path=None, builder=None,
                              template_dir=None):
    """
    This function is adapted from
    `sphinx.ext.autosummary.generate.generate_autosummmary_docs` to generate
    source for the automodsumm directives that should be autosummarized.
    Unlike generate_autosummary_docs, this function is called one file at a
    time.
    """
    import os

    from sphinx import package_dir
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

    info('[automodapi] generating automodsumm for: ' + srcfn)

    # create our own templating environment
    template_dirs = [os.path.join(package_dir, 'ext',
                                  'autosummary', 'templates')]
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

#    gennms = [item[0] for item in items]
#    if len(gennms) > 20:
#        gennms = gennms[:10] + ['...'] + gennms[-10:]
#    info('[automodapi] generating autosummary for: ' + ', '.join(gennms))

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
            warn('[automodapi] failed to import %r: %s' % (name, e))
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
                try:
                    template = template_env.get_template('autosummary/%s.rst'
                                                         % doc.objtype)
                except TemplateNotFound:
                    template = template_env.get_template('autosummary/base.rst')

            def get_members(obj, typ, include_public=[]):
                items = []
                for name in dir(obj):
                    try:
                        documenter = get_documenter(safe_getattr(obj, name),
                                                    obj)
                    except AttributeError:
                        continue
                    if documenter.objtype == typ:
                        items.append(name)
                public = [x for x in items
                          if x in include_public or not x.startswith('_')]
                return public, items

            ns = {}

            if doc.objtype == 'module':
                ns['members'] = dir(obj)
                ns['functions'], ns['all_functions'] = \
                                   get_members(obj, 'function')
                ns['classes'], ns['all_classes'] = \
                                 get_members(obj, 'class')
                ns['exceptions'], ns['all_exceptions'] = \
                                   get_members(obj, 'exception')
            elif doc.objtype == 'class':
                ns['members'] = dir(obj)
                ns['methods'], ns['all_methods'] = \
                                 get_members(obj, 'method', ['__init__'])
                ns['attributes'], ns['all_attributes'] = \
                                 get_members(obj, 'attribute')

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

            rendered = template.render(**ns)
            f.write(rendered)
        finally:
            f.close()


def setup(app):
    # need autosummary for automodapi
    app.setup_extension('sphinx.ext.autosummary')
    app.add_directive('automodsumm', Automodsumm)

    #need inheritance-diagram for automod-diagram
    app.setup_extension('sphinx.ext.inheritance_diagram')
    app.add_directive('automod-diagram', Automoddiagram)

    app.connect('source-read', process_automodapi)
    app.connect('builder-inited', process_automodsumm_generation)

    app.add_config_value('automodsumm_generate', False, True)
