# -*- coding: utf-8 -*-
"""
This extension makes it easy to edit documentation on github.

It adds links associated with each docstring that go to the
corresponding view source page on Github.  From there, the user can
push the "Edit" button, edit the docstring, and submit a pull request.

It has the following configuration parameters:

edit_on_github_project:
    The name of the github project, in the form
    "username/projectname".

edit_on_github_branch:
    The name of the branch to edit.  If this is a released version,
    this should be a git tag referring to that version.  For a
    dev version, it often makes sense for it to be "master".  It
    may also be a git hash.

edit_on_github_source_root:
    The location within the source tree of the root of the
    Python package.  Defaults to "lib".

edit_on_github_doc_root:
    The location within the source tree of the root of the
    documentation source.  Defaults to "doc", but it may make sense to
    set it to "doc/source" if the project uses a separate source
    directory.

edit_on_github_docstring_message:
    The phrase displayed in the links to edit a docstring.

edit_on_github_page_message:
    The phrase displayed in the links to edit a RST page.
"""
import inspect
import os
import sys

from docutils import nodes

from sphinx import addnodes


def import_object(modname, name):
    """
    Import the object given by *modname* and *name* and return it.
    If not found, or the import fails, returns None.
    """
    try:
        __import__(modname)
        mod = sys.modules[modname]
        obj = mod
        for part in name.split('.'):
            obj = getattr(obj, part)
        return obj
    except:
        return None


def doctree_read(app, doctree):
    if app.config.edit_on_github_project == 'REQUIRED':
        raise ValueError(
            "The edit_on_github_project configuration variable must be "
            "provided in the conf.py")

    source_root = app.config.edit_on_github_source_root
    if source_root != '' and not source_root.endswith('/'):
        source_root += '/'
    doc_root = app.config.edit_on_github_doc_root
    if doc_root != '' and not doc_root.endswith('/'):
        doc_root += '/'
    url = 'http://github.com/%s/tree/%s/' % (
        app.config.edit_on_github_project,
        app.config.edit_on_github_branch)

    docstring_message = app.config.edit_on_github_docstring_message
    page_message = app.config.edit_on_github_page_message

    doc_path = os.path.relpath(doctree.get('source'), app.builder.srcdir)
    path = url + doc_root + doc_path
    section = nodes.section()
    para = nodes.paragraph()
    section += para
    onlynode = addnodes.only(expr='html')
    para += onlynode
    onlynode += nodes.reference('', refuri=path)
    onlynode[0] += nodes.inline(
        '', page_message, classes=['edit-on-github'])
    doctree += section

    for objnode in doctree.traverse(addnodes.desc):
        if objnode.get('domain') != 'py':
            continue
        names = set()
        for signode in objnode:
            if not isinstance(signode, addnodes.desc_signature):
                continue
            modname = signode.get('module')
            if not modname:
                continue
            fullname = signode.get('fullname')
            if fullname in names:
                # only one link per name, please
                continue
            names.add(fullname)
            obj = import_object(modname, fullname)
            anchor = None
            if obj is not None:
                try:
                    lines, lineno = inspect.getsourcelines(obj)
                except:
                    pass
                else:
                    anchor = '#L%d' % lineno
            if anchor:
                path = '%s%s%s.py%s' % (
                    url, source_root, modname.replace('.', '/'), anchor)
                onlynode = addnodes.only(expr='html')
                onlynode += nodes.reference('', refuri=path)
                onlynode[0] += nodes.inline(
                    '', '', nodes.raw('', '&nbsp;', format='html'),
                    nodes.Text(docstring_message),
                    classes=['edit-on-github', 'viewcode-link'])
                signode += onlynode


def setup(app):
    app.add_config_value('edit_on_github_project', 'REQUIRED', True)
    app.add_config_value('edit_on_github_branch', 'master', True)
    app.add_config_value('edit_on_github_source_root', 'lib', True)
    app.add_config_value('edit_on_github_doc_root', 'doc', True)
    app.add_config_value('edit_on_github_docstring_message',
                         '[edit on github]', True)
    app.add_config_value('edit_on_github_page_message',
                         '[edit this page on github]', True)

    app.connect('doctree-read', doctree_read)
