# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Extension for Sphinx to check the completness of docstrings.

This function compares the signature and the parameters documented in the 
doctring for every method and every function in the sphinx documentation.

It it limited to functions and methods that actually appear in the documentation.
Private members that are not be referenced in the sphinx API docs, will not be 
checked. Usually (unless sphinx is instructed otherwise), this also means that
special functions (such as ``__getattr__``) will not be included in this check.

For every methods and function checked the 
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
                # print ''
                # print node
                docstring = []
                signature = []
                for subnode in node.children[0].traverse(addnodes.desc_parameterlist):
                    signature = [str(t.astext().split('=')[0]) for t 
                                 in subnode.traverse(docutils.nodes.Text)]
                for fields in node.children[1].traverse(docutils.nodes.field):
                    if fields.children[0].astext() == 'Parameters':
                        docstring = [str(f.astext()) for f 
                                 in fields.traverse(docutils.nodes.strong)]
                    # All methods are listed within a class and then again as
                    # methods. Thus, classes have several field lists, but
                    # only the first one is for __init__ and needs to be checked here.
                    break
                # Remove *args and **kwargs from the signature
                signature_required = [s for s in signature if (s[0]!='*')]
                n_pars = len(signature)
                n_pars_named = len(signature_required)
                # split docstring where several params are described on one line
                temp = []
                for x in docstring:
                    temp.extend(x.split(','))
                docstring = [x.strip() for x in temp]

                # print(signature)
                # print(signature_required)
                # print(docstring)
                # print(n_pars, n_pars_named)
                # print(((n_pars_named == n_pars) and (signature == docstring)))
                # print(((n_pars_named < n_pars) and
                #      (signature_required == docstring[:n_pars_named])))
                # print(not(
                #     # no *args or **kwargs
                #     ((n_pars_named == n_pars) and (signature == docstring))
                #     or
                #     # *args or **kswargs
                #     ((n_pars_named < n_pars) and
                #      (signature_required == docstring[:n_pars_named]))
                #     ))
                if not(
                    # no *args or **kwargs
                    ((n_pars_named == n_pars) and (signature == docstring))
                    or
                    # *args or **kswargs
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

