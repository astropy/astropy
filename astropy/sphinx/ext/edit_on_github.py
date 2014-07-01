# **Please Note**: ``astropy.sphinx`` exists only for backward-compatibility
# purposes - it has now been moved to the separate astropy-helpers package,
# located at https://github.com/astropy/astropy-helpers. Any new development or
# bug fixes should be done there.
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This extension makes it easy to edit documentation on github.

It adds links associated with each docstring that go to the
corresponding view source page on Github.  From there, the user can
push the "Edit" button, edit the docstring, and submit a pull request.

It has the following configuration options (to be set in the project's
``conf.py``):

* ``edit_on_github_project``
    The name of the github project, in the form
    "username/projectname".

* ``edit_on_github_branch``
    The name of the branch to edit.  If this is a released version,
    this should be a git tag referring to that version.  For a
    dev version, it often makes sense for it to be "master".  It
    may also be a git hash.

* ``edit_on_github_source_root``
    The location within the source tree of the root of the
    Python package.  Defaults to "lib".

* ``edit_on_github_doc_root``
    The location within the source tree of the root of the
    documentation source.  Defaults to "doc", but it may make sense to
    set it to "doc/source" if the project uses a separate source
    directory.

* ``edit_on_github_docstring_message``
    The phrase displayed in the links to edit a docstring.  Defaults
    to "[edit on github]".

* ``edit_on_github_page_message``
    The phrase displayed in the links to edit a RST page.  Defaults
    to "[edit this page on github]".

* ``edit_on_github_help_message``
    The phrase displayed as a tooltip on the edit links.  Defaults to
    "Push the Edit button on the next page"

* ``edit_on_github_skip_regex``
    When the path to the .rst file matches this regular expression,
    no "edit this page on github" link will be added.  Defaults to
    ``"_.*"``.
"""
import inspect
import os
import re
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


def get_url_base(app):
    return  'http://github.com/%s/tree/%s/' % (
        app.config.edit_on_github_project,
        app.config.edit_on_github_branch)


def doctree_read(app, doctree):
    # Get the configuration parameters
    if app.config.edit_on_github_project == 'REQUIRED':
        raise ValueError(
            "The edit_on_github_project configuration variable must be "
            "provided in the conf.py")

    source_root = app.config.edit_on_github_source_root
    url = get_url_base(app)

    docstring_message = app.config.edit_on_github_docstring_message

    # Handle the docstring-editing links
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
                real_modname = inspect.getmodule(obj).__name__
                path = '%s%s%s.py%s' % (
                    url, source_root, real_modname.replace('.', '/'), anchor)
                onlynode = addnodes.only(expr='html')
                onlynode += nodes.reference(
                    reftitle=app.config.edit_on_github_help_message,
                    refuri=path)
                onlynode[0] += nodes.inline(
                    '', '', nodes.raw('', '&nbsp;', format='html'),
                    nodes.Text(docstring_message),
                    classes=['edit-on-github', 'viewcode-link'])
                signode += onlynode


def html_page_context(app, pagename, templatename, context, doctree):
    if (templatename == 'page.html' and
        not re.match(app.config.edit_on_github_skip_regex, pagename)):

        doc_root = app.config.edit_on_github_doc_root
        if doc_root != '' and not doc_root.endswith('/'):
            doc_root += '/'
        doc_path = os.path.relpath(doctree.get('source'), app.builder.srcdir)
        url = get_url_base(app)

        page_message = app.config.edit_on_github_page_message

        context['edit_on_github'] = url + doc_root + doc_path
        context['edit_on_github_page_message'] = (
            app.config.edit_on_github_page_message)


def setup(app):
    app.add_config_value('edit_on_github_project', 'REQUIRED', True)
    app.add_config_value('edit_on_github_branch', 'master', True)
    app.add_config_value('edit_on_github_source_root', 'lib', True)
    app.add_config_value('edit_on_github_doc_root', 'doc', True)
    app.add_config_value('edit_on_github_docstring_message',
                         '[edit on github]', True)
    app.add_config_value('edit_on_github_page_message',
                         'Edit This Page on Github', True)
    app.add_config_value('edit_on_github_help_message',
                         'Push the Edit button on the next page', True)
    app.add_config_value('edit_on_github_skip_regex',
                         '_.*', True)

    app.connect('doctree-read', doctree_read)
    app.connect('html-page-context', html_page_context)
