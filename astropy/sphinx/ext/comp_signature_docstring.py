# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Extension for Sphinx to check the completness of docstrings.

This function compares the signature and the parameters documented in the 
doctring for every method and every function in the sphinx documentation.

It it limited to functions and methods that actually appear in the documentation.
Private members that are not be referenced in the sphinx API docs, will not be 
checked. Usually (unless sphinx is instructed otherwise), this also means that
special functions (such as ``__getattr__``) will not be included in this check.

For every method and function checked the parameters in the calling signature
and in the docstring are compared. 
The docstring is parsed for every field list whose title contains the string
`'Parameters'` (such as "Parameters" or "Other Parameters" which are common
with numpydoc.

In order to pass the test, the following conditions are tested::

- The signature matches the documented arguments in full.
- If the signature contains `*args` and/or `**kwargs`, all named signature
  arguments must match the documented arguments in full. Further entries in the
  docstring are allowed.
- Private arguments (starting with `_`) in the signature are ignored and should
  not appear in the documentation.

Example
-------
sudo python setup.py build_sphinx -b signaturedocstring
"""
from __future__ import print_function

import sys

import docutils.nodes
from sphinx.builders import Builder
from sphinx import addnodes

import inspect



def setup(app):
    # Inside the test suite, "app" will be a module - do I need this?
    if inspect.ismodule(app):
        return
    app.info('Initializing Documentation checker')
    app.add_builder(MatchSignatureDocstringBuilder)
    app.info('If the build succeeded, all parameters in all calling signature are documented.')
    return


class MatchSignatureDocstringBuilder(Builder):
    """Checks if documented paremeters match calling signatures.
    """
    name = 'signaturedocstring'

    def get_outdated_docs(self):
        return 'all documents'

    def prepare_writing(self, docnames):
        return

    def get_target_uri(self, docname, typ=None):
        return ''


    def write_doc(self, docname, doctree):
        print('')
        for node in doctree.traverse(addnodes.desc):
            if node['desctype'] in ['function', 'method','class']:
                docstring = []
                signature = []
                for subnode in node.children[0].traverse(addnodes.desc_parameterlist):
                    signature = [str(t.astext().split('=')[0]) for t 
                                 in subnode.traverse(docutils.nodes.Text)]
                for c in node.children[1].children:
                    # If node is a class, if contains info on all
                    # methods, but at a lower level of hirachy.
                    if isinstance(c, docutils.nodes.field_list):
                        for fields in c.traverse(docutils.nodes.field):
                            if 'Parameters'in fields.children[0].astext():
                                docstring.extend([str(f.astext()) for f 
                                                  in fields.traverse(docutils.nodes.strong)])
                # Remove *args, **kwargs and private keywords from the signature
                signature_required = [s for s in signature if (s[0] not in '*_')]
                n_pars = len(signature)
                n_pars_named = len(signature_required)

                if not(
                    # no *args or **kwargs or other special values
                    ((n_pars_named == n_pars) and (signature == docstring))
                    or
                    # *args or **kwargs
                    ((n_pars_named < n_pars) and
                     (signature_required == docstring[:n_pars_named]))
                    ):
                    # Where are we?
                    if node.children[0]['ids']:
                        name = node.children[0]['ids'][0]
                    else:
                        name = node.children[0]['module'][0] + '.' + node.children[0]['fullname'][0]
                    print('{0}: Parameters in signature and docstring do not match'.format(name))
                    print('    Parameters in signature:', signature)
                    print('    Parameters in docstring:', docstring)
                    # This parameter list is inconsistent with the documentation.
                    # Set status code for error.
                    self.app.statuscode = 1
        return


