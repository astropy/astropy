# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module's main purpose is to act as a script to create new versions
of erfa.pyx when ERFA is updated (or this generator is enhanced).

`Jinja2 <http://jinja.pocoo.org/>`_ must be installed for this
module/script to function.
"""

import re
import os.path


__all__ = ['Function']


ctype_to_dtype = {'double'       : "numpy.double",
                  'double *'     : "numpy.double",
                  'int'          : "numpy.int",
                  'int *'        : "numpy.int",
                  'int[4]'       : "numpy.dtype([('', 'i', (4,))])",
                  'double[2]'    : "numpy.dtype([('', 'd', (2,))])",
                  'double[3]'    : "numpy.dtype([('p', 'd', (3,))])",
                  'double[2][3]' : "numpy.dtype([('pv', 'd', (2,3))])",
                  'double[3][3]' : "numpy.dtype([('r', 'd', (3,3))])",
                  'eraASTROM *'  : "dt_eraASTROM",
                  'char *'       : "numpy.dtype('S1')",
                  }


class FunctionDoc(object):

    def __init__(self, doc):
        self.doc = doc.replace("**", "  ").replace("/*", "  ").replace("*/", "  ")
        self.__input = None
        self.__output = None

    @property
    def input(self):
        if self.__input is None:
            self.__input = []
            result = re.search("Given([^\n]*):\n(.+?)  \n", self.doc, re.DOTALL)
            if result is not None:
                __input = result.group(2)
                for i in __input.split("\n"):
                    arg_doc = ArgumentDoc(i)
                    if arg_doc.name is not None:
                        self.__input.append(arg_doc)
            result = re.search("Given and returned([^\n]*):\n(.+?)  \n", self.doc, re.DOTALL)
            if result is not None:
                __input = result.group(2)
                for i in __input.split("\n"):
                    arg_doc = ArgumentDoc(i)
                    if arg_doc.name is not None:
                        self.__input.append(arg_doc)
        return self.__input

    @property
    def output(self):
        if self.__output is None:
            self.__output = []
            result = re.search("Returned([^\n]*):\n(.+?)  \n", self.doc, re.DOTALL)
            if result is not None:
                __output = result.group(2)
                for i in __output.split("\n"):
                    arg_doc = ArgumentDoc(i)
                    if arg_doc.name is not None:
                        self.__output.append(arg_doc)
            result = re.search("Given and returned([^\n]*):\n(.+?)  \n", self.doc, re.DOTALL)
            if result is not None:
                __output = result.group(2)
                for i in __output.split("\n"):
                    arg_doc = ArgumentDoc(i)
                    if arg_doc.name is not None:
                        self.__output.append(arg_doc)
        return self.__output

    def __repr__(self):
        return self.doc


class ArgumentDoc(object):

    def __init__(self, doc):
        match = re.search("^ +([^ ]+)[ ]+([^ ]+)[ ]+(.+)", doc)
        if match is not None:
            self.name = match.group(1)
            self.type = match.group(2)
            self.doc = match.group(3)
        else:
            self.name = None
            self.type = None
            self.doc = None

    def __repr__(self):
        return "    {0:15} {1:15} {2}".format(self.name, self.type, self.doc)


class Argument(object):

    def __init__(self, definition, doc):
        self.__doc = doc
        self.__inout_state = None
        self.definition = definition.strip()
        if "*" in self.definition:
            self.ctype, self.name = self.definition.split("*", 1)
            self.ctype += "*"
        else:
            self.ctype, self.name = self.definition.rsplit(" ", 1)
            if "[" in self.name:
                self.name, arr = self.name.split("[", 1)
                self.ctype += ("["+arr)

    @property
    def inout_state(self):
        if self.__inout_state is None:
            self.__inout_state = ''
            for i in self.__doc.input:
                if self.name in i.name.split(','):
                    self.__inout_state = 'in'
            for o in self.__doc.output:
                if self.name in o.name.split(','):
                    if self.__inout_state == 'in':
                        self.__inout_state = 'inout'
                    else:
                        self.__inout_state = 'out'
        return self.__inout_state

    @property
    def ctype_ptr(self):
        if self.ctype[-1] == ']':
            return self.ctype.split('[')[0]+" *"
        elif self.ctype[:6] == 'const ':
            return self.ctype[6:]
        else:
            return self.ctype

    @property
    def dtype(self):
        return ctype_to_dtype[self.ctype]

    def __repr__(self):
        return "Argument('{0}', name='{1}', ctype='{2}', inout_state='{3}')".format(self.definition, self.name, self.ctype, self.inout_state)


class Return(object):

    def __init__(self, ctype, doc):
        self.name = 'ret'
        self.ctype = ctype
        self.inout_state = 'ret'
        self.ctype_ptr = ctype

    @property
    def dtype(self):
        return ctype_to_dtype[self.ctype]


class Function(object):

    def __init__(self, name, source_path):
        self.name = name
        self.pyname = name.split('era')[-1].lower()
        self.filename = name.split("era")[-1].lower()+".c"
        self.filepath = os.path.join(os.path.normpath(source_path), self.filename)
        pattern = "\n([^\n]+{0} ?\([^)]+\)).+?(/\*.+?\*/)".format(name)
        p = re.compile(pattern, flags=re.DOTALL|re.MULTILINE)
        with open(self.filepath) as f:
            search = p.search(f.read())
            self.cfunc = search.group(1)
            self.__doc = FunctionDoc(search.group(2))
        self.args = []
        for arg in re.search("\(([^)]+)\)", self.cfunc, flags=re.MULTILINE|re.DOTALL).group(1).split(','):
            self.args.append(Argument(arg, self.__doc))
        self.ret = re.search("^(.*){0}".format(name), self.cfunc).group(1).strip()
        if self.ret == 'double':
            self.args.append(Return(self.ret, self.__doc))

    def args_by_inout(self, inout_filter, prop=None, join=None):
        result = []
        for arg in self.args:
            if arg.inout_state in inout_filter.split('|'):
                if prop is None:
                    result.append(arg)
                else:
                    result.append(getattr(arg, prop))
        if join is not None:
            return join.join(result)
        else:
            return result


def main(srcdir):
    from jinja2 import Environment, PackageLoader

    #Prepare the jinja2 templating environment
    env = Environment(loader=PackageLoader('cython_numpy_auto', '.'))

    def prefix(a_list, pre):
        return [pre+'{0}'.format(an_element) for an_element in a_list]
    def postfix(a_list, post):
        return ['{0}'.format(an_element)+post for an_element in a_list]
    def surround(a_list, pre, post):
        return [pre+'{0}'.format(an_element)+post for an_element in a_list]
    env.filters['prefix'] = prefix
    env.filters['postfix'] = postfix
    env.filters['surround'] = surround

    erfa_pyx_in = env.get_template('erfa.pyx.in')

    #Extract all the ERFA function names from erfa.h
    with open(srcdir + "/erfa.h", "r") as f:

        erfa_h = f.read()

        funcs = []
        section_subsection_functions = re.findall('/\* (\w*)/(\w*) \*/\n(.*?)\n\n',
                                                 erfa_h, flags=re.DOTALL|re.MULTILINE)
        for section, subsection, functions in section_subsection_functions:
            print("{0}.{1}".format(section, subsection))
            if section == "Astronomy":
                func_names = re.findall(' (\w+)\(.*?\);', functions, flags=re.DOTALL)
                for name in func_names:
                    print("{0}.{1}.{2}...".format(section, subsection, name))
                    funcs.append(Function(name, srcdir))
        print("Done!")
        #Render the template and save
        erfa_pyx = erfa_pyx_in.render(funcs=funcs)
        with open("erfa.pyx", "w") as f:
            f.write(erfa_pyx)
