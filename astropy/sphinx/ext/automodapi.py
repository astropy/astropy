"""
This sphinx extensions adds tools to simplify generating the API documentation
for Astropy packages.

It adds three main pieces:

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


"""
import re

from sphinx.ext.autosummary import Autosummary
from sphinx.ext.inheritance_diagram import InheritanceDiagram


class Automodsumm(Autosummary):
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = False
    has_content = False
    option_spec = dict(Autosummary.option_spec)

    def __init__(self, name, arguments, options, content, lineno,
                 content_offset, block_text, state, state_machine):
        if 'toctree' in options and options['toctree'] == '':
            options['toctree'] = '_generated/'
        super(Autosummary, self).__init__(
              name, arguments, options, content, lineno, content_offset,
              block_text, state, state_machine)

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
    fullnames = [nmtuple[0] for nmtuple in sortnms]
    attrnames = [nmtuple[1] for nmtuple in sortnms]

    if names:
        return [n for n in sorted(fullnames) if n is not None]
    else:
        return [pkgd[n] for n in sorted(attrnames) if n is not None]


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
"""
automod_inh_templ = """
Class Inheritance Diagram
^^^^^^^^^^^^^^^^^^^^^^^^^

.. automod-diagram:: {modname}
    :private-bases:

"""

_automodapirex = re.compile(
    r'^(\s*\.\.\s+automodapi::\s*)([A-Za-z0-9_.]+)\s*$', flags=re.MULTILINE)


def process_automodapi(app, docname, source):
    from inspect import isclass

    spl = _automodapirex.split(source[0])
    if len(spl) > 1:  # automodsumm is in this document
        for grp in reversed(range(len(spl) // 3)):
            modnm = spl[grp * 3 + 2]
            del spl[grp * 3 + 2]

            templstr = automod_templ.format(modname=modnm,
                                            modcrts='^' * len(modnm))

            # add in the inheritance diagram if any classes are in the module
            if any([isclass(obj) for obj in find_mod_objs(modnm, False)]):
                templinhstr = automod_inh_templ.format(modname=modnm)
                spl[grp * 3 + 1] = templstr + templinhstr
            else:
                spl[grp * 3 + 1] = templstr

        #now re-assign the source string to the modified version
        source[0] = ''.join(spl)


def setup(app):
    # need autosummary for automodapi
    app.setup_extension('sphinx.ext.autosummary')
    app.add_directive('automodsumm', Automodsumm)

    #need inheritance-diagram for automod-diagram
    app.setup_extension('sphinx.ext.inheritance_diagram')
    app.add_directive('automod-diagram', Automoddiagram)

    app.connect('source-read', process_automodapi)
