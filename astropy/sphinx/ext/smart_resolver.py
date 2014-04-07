# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The classes in the astropy docs are documented by their API location,
which is not necessarily where they are defined in the source.  This
causes a problem when certain automated features of the doc build,
such as the inheritance diagrams or the `Bases` list of a class
reference a class by its canonical location rather than its "user"
location.

In the `autodoc-process-docstring` event, a mapping from the actual
name to the API name is maintained.  Later, in the `missing-reference`
enent, unresolved references are looked up in this dictionary and
corrected if possible.
"""

from docutils.nodes import literal


def process_docstring(app, what, name, obj, options, lines):
    if what in ('class', 'exception'):
        env = app.env
        if not hasattr(env, 'class_name_mapping'):
            env.class_name_mapping = {}
        mapping = env.class_name_mapping
        mapping[obj.__module__ + '.' + obj.__name__] = name


def missing_reference_handler(app, env, node, contnode):
    if not hasattr(env, 'class_name_mapping'):
        env.class_name_mapping = {}
    mapping = env.class_name_mapping
    reftype = node['reftype']
    reftarget = node['reftarget']
    if reftype in ('obj', 'class', 'exc', 'meth'):
        reftarget = node['reftarget']
        suffix = ''
        if reftarget not in mapping:
            if reftype in ('obj', 'meth') and '.' in reftarget:
                front, suffix = reftarget.rsplit('.', 1)
                if front in mapping:
                    reftarget = front
                    suffix = '.' + suffix

        if reftarget in mapping:
            newtarget = mapping[reftarget] + suffix
            if not node['refexplicit'] and not '~' in node.rawsource:
                contnode = literal(text=newtarget)
            newnode = env.domains['py'].resolve_xref(
                env, node['refdoc'], app.builder, 'class', newtarget,
                node, contnode)
            newnode['reftitle'] = reftarget
            return newnode


def setup(app):
    app.connect('autodoc-process-docstring', process_docstring)

    app.connect('missing-reference', missing_reference_handler)
