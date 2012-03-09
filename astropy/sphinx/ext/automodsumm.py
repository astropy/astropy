# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This sphinx extension adds two directives for summarizing the public
members of a module or package.

These directives are primarily for use with the `automodapi` extension,
but can be used independently.

=======================
`automodsumm` directive
=======================

This directive will produce an "autosummary"-style table for public
attributesof a specified module. See the `sphinx.ext.autosummary`
extensionfor details on this process. The main difference from the
`autosummary`directive is that `autosummary` requires manually inputting
all attributes that appear in the table, while this captures the entries
automatically.

This directive requires a single argument that must be a module or
package. It also accepts the same arguments as the `autosummary`
directive- see `sphinx.ext.autosummary` for details.


===========================
`automod-diagram` directive
===========================

This directive will produce an inheritance diagram like that of the
`sphinx.ext.inheritance_diagram`extension.

This directive requires a single argument that must be a module or
package.It accepts no options.

.. note::
    Like 'inheritance-diagram', 'automod-diagram' requires graphviz to
    generate the inheritance diagram.

"""
import re

from sphinx.ext.autosummary import Autosummary
from sphinx.ext.inheritance_diagram import InheritanceDiagram

from ...utils.misc import find_mod_objs


class Automodsumm(Autosummary):
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = False
    has_content = False
    #option_spec = dict(Autosummary.option_spec)

    def run(self):
        try:
            modnms = find_mod_objs(self.arguments[0])[1]
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


#<-------------------automod-diagram stuff------------------------------------>
class Automoddiagram(InheritanceDiagram):
    def run(self):
        from inspect import isclass

        try:
            nms, objs = find_mod_objs(self.arguments[0], onlylocals=True)[1:]
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


#<---------------------automodsumm generation stuff--------------------------->


def process_automodsumm_generation(app):
    import os

    env = app.builder.env
    ext = app.config.source_suffix

    filestosearch = [x + ext for x in env.found_docs
                     if os.path.isfile(env.doc2path(x))]\

    liness = [automodsumm_to_autosummary_lines(sfn, app) for
              sfn in filestosearch]

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
    import os

    fullfn = os.path.join(app.builder.env.srcdir, fn)

    lines = []
    singleindent = ' ' * 4

    with open(fullfn) as fr:
        if 'astropy.sphinx.ext.automodapi' in app._extensions:
            from astropy.sphinx.ext.automodapi import automodapi_replace
            # Must do the automodapi on the source to get the automodsumm
            # that might be in there
            filestr = automodapi_replace(fr.read(), app, True, fn, False)
        else:
            filestr = fr.read()

    indentrex = None
    modnm = None
    indent = ''

    def finish_content(lines, modnm, indent):
        """Adds the content items for autosummary to the lines"""
        lines.append(indent + singleindent)
        for objnm in find_mod_objs(modnm, onlylocals=True)[1]:
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
                tmplstr = 'autosummary/%s.rst'
                try:
                    template = template_env.get_template(tmplstr % doc.objtype)
                except TemplateNotFound:
                    template = template_env.get_template(tmplstr % 'base')

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
