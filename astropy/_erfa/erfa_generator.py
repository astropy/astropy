# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module's main purpose is to act as a script to create new versions
of erfa.c when ERFA is updated (or this generator is enhanced).

`Jinja2 <http://jinja.pocoo.org/>`_ must be installed for this
module/script to function.

Note that this does *not* currently automate the process of creating structs
or dtypes for those structs.  They should be added manually in the template file.
"""
from __future__ import absolute_import, division, print_function
# note that we do *not* use unicode_literals here, because that makes the
# generated code's strings have u'' in them on py 2.x

import re
import os.path
from collections import OrderedDict


ctype_to_dtype = {'double'     : "numpy.double",
                  'int'        : "numpy.intc",
                  'eraASTROM'  : "dt_eraASTROM",
                  'eraLDBODY'  : "dt_eraLDBODY",
                  'char'       : "numpy.dtype('S16')",
                  'const char' : "numpy.dtype('S16')",
                  }


NDIMS_REX = re.compile(re.escape("numpy.dtype([('fi0', '.*', <(.*)>)])").replace(r'\.\*','.*').replace(r'\<', '(').replace(r'\>',')'))


class FunctionDoc(object):

    def __init__(self, doc):
        self.doc = doc.replace("**", "  ").replace("/*\n", "").replace("*/", "")
        self.__input = None
        self.__output = None
        self.__ret_info = None

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

    @property
    def ret_info(self):
        if self.__ret_info is None:
            ret_info = []
            result = re.search("Returned \\(function value\\)([^\n]*):\n(.+?)  \n", self.doc, re.DOTALL)
            if result is not None:
                ret_info.append(ReturnDoc(result.group(2)))

            if len(ret_info) == 0:
                self.__ret_info = ''
            elif len(ret_info) == 1:
                self.__ret_info = ret_info[0]
            else:
                raise ValueError("Multiple C return sections found in this doc:\n" + self.doc)

        return self.__ret_info

    def __repr__(self):
        return self.doc.replace("  \n", "\n")


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
        self.doc = doc
        self.__inout_state = None
        self.ctype, ptr_name_arr = definition.strip().rsplit(" ", 1)
        if "*" == ptr_name_arr[0]:
            self.is_ptr = True
            name_arr = ptr_name_arr[1:]
        else:
            self.is_ptr = False
            name_arr = ptr_name_arr
        if "[]" in ptr_name_arr:
            self.is_ptr = True
            name_arr = name_arr[:-2]
        if "[" in name_arr:
            self.name, arr = name_arr.split("[", 1)
            self.shape = tuple([int(size) for size in arr[:-1].split("][")])
        else:
            self.name = name_arr
            self.shape = ()

    @property
    def inout_state(self):
        if self.__inout_state is None:
            self.__inout_state = ''
            for i in self.doc.input:
                if self.name in i.name.split(','):
                    self.__inout_state = 'in'
            for o in self.doc.output:
                if self.name in o.name.split(','):
                    if self.__inout_state == 'in':
                        self.__inout_state = 'inout'
                    else:
                        self.__inout_state = 'out'
        return self.__inout_state

    @property
    def ctype_ptr(self):
        if (self.is_ptr) | (len(self.shape)>0):
            return self.ctype+" *"
        else:
            return self.ctype

    @property
    def name_in_broadcast(self):
        if len(self.shape)>0:
            return "{0}_in[...{1}]".format(self.name, ",0"*len(self.shape))
        else:
            return "{0}_in".format(self.name)

    @property
    def name_out_broadcast(self):
        if len(self.shape)>0:
            return "{0}_out[...{1}]".format(self.name, ",0"*len(self.shape))
        else:
            return "{0}_out".format(self.name)

    @property
    def dtype(self):
        return ctype_to_dtype[self.ctype]

    @property
    def ndim(self):
        return len(self.shape)

    @property
    def cshape(self):
        return ''.join(['[{0}]'.format(s) for s in self.shape])

    @property
    def name_for_call(self):
        if self.is_ptr:
            return '_'+self.name
        else:
            return '*_'+self.name

    def __repr__(self):
        return "Argument('{0}', name='{1}', ctype='{2}', inout_state='{3}')".format(self.definition, self.name, self.ctype, self.inout_state)


class ReturnDoc(object):

    def __init__(self, doc):
        self.doc = doc

        self.infoline = doc.split('\n')[0].strip()
        self.type = self.infoline.split()[0]
        self.descr = self.infoline.split()[1]

        if self.descr.startswith('status'):
            self.statuscodes = statuscodes = {}

            code = None
            for line in doc[doc.index(':')+1:].split('\n'):
                ls = line.strip()
                if ls != '':
                    if ' = ' in ls:
                        code, msg = ls.split(' = ')
                        if code != 'else':
                            code = int(code)
                        statuscodes[code] = msg
                    elif code is not None:
                        statuscodes[code] += ls
        else:
            self.statuscodes = None

    def __repr__(self):
        return "Return value, type={0:15}, {1}, {2}".format(self.type, self.descr, self.doc)


class Return(object):

    def __init__(self, ctype, doc):
        self.name = 'c_retval'
        self.name_out_broadcast = self.name+"_out"
        self.inout_state = 'stat' if ctype == 'int' else 'ret'
        self.ctype = ctype
        self.ctype_ptr = ctype
        self.shape = ()
        self.doc = doc

    def __repr__(self):
        return "Return(name='{0}', ctype='{1}', inout_state='{2}')".format(self.name, self.ctype, self.inout_state)

    @property
    def dtype(self):
        return ctype_to_dtype[self.ctype]

    @property
    def nd_dtype(self):
        """
        This if the return type has a multi-dimensional output, like
        double[3][3]
        """
        return "'fi0'" in self.dtype

    @property
    def doc_info(self):
        return self.doc.ret_info


class Function(object):
    """
    A class representing a C function.

    Parameters
    ----------
    name : str
        The name of the function
    source_path : str
        Either a directory, which means look for the function in a
        stand-alone file (like for the standard ERFA distribution), or a
        file, which means look for the function in that file (as for the
        astropy-packaged single-file erfa.c).
    match_line : str, optional
        If given, searching of the source file will skip until it finds
        a line matching this string, and start from there.
    """

    def __init__(self, name, source_path, match_line=None):
        self.name = name
        self.pyname = name.split('era')[-1].lower()
        self.filename = self.pyname+".c"
        if os.path.isdir(source_path):
            self.filepath = os.path.join(os.path.normpath(source_path), self.filename)
        else:
            self.filepath = source_path

        with open(self.filepath) as f:
            if match_line:
                line = f.readline()
                while line != '':
                    if line.startswith(match_line):
                        filecontents = '\n' + line + f.read()
                        break
                    line = f.readline()
                else:
                    msg = ('Could not find the match_line "{0}" in '
                           'the source file "{1}"')
                    raise ValueError(msg.format(match_line, self.filepath))
            else:
                filecontents = f.read()

        pattern = "\n([^\n]+{0} ?\([^)]+\)).+?(/\*.+?\*/)".format(name)
        p = re.compile(pattern, flags=re.DOTALL|re.MULTILINE)

        search = p.search(filecontents)
        self.cfunc = " ".join(search.group(1).split())
        self.doc = FunctionDoc(search.group(2))

        self.args = []
        for arg in re.search("\(([^)]+)\)", self.cfunc).group(1).split(', '):
            self.args.append(Argument(arg, self.doc))
        self.ret = re.search("^(.*){0}".format(name), self.cfunc).group(1).strip()
        if self.ret != 'void':
            self.args.append(Return(self.ret, self.doc))

    def args_by_inout(self, inout_filter, prop=None, join=None):
        """
        Gives all of the arguments and/or returned values, depending on whether
        they are inputs, outputs, etc.

        The value for `inout_filter` should be a string containing anything
        that arguments' `inout_state` attribute produces.  Currently, that can be:

          * "in" : input
          * "out" : output
          * "inout" : something that's could be input or output (e.g. a struct)
          * "ret" : the return value of the C function
          * "stat" : the return value of the C function if it is a status code

        It can also be a "|"-separated string giving inout states to OR
        together.
        """
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

    def __repr__(self):
        return "Function(name='{0}', pyname='{1}', filename='{2}', filepath='{3}')".format(self.name, self.pyname, self.filename, self.filepath)


class Constant(object):

    def __init__(self, name, value, doc):
        self.name = name.replace("ERFA_","")
        self.value = value.replace("ERFA_","")
        self.doc = doc


def main(srcdir, outfn, templateloc, verbose=True):
    from jinja2 import Environment, FileSystemLoader

    if verbose:
        print_ = lambda *args, **kwargs: print(*args, **kwargs)
    else:
        print_ = lambda *args, **kwargs: None

    #Prepare the jinja2 templating environment
    env = Environment(loader=FileSystemLoader(templateloc))

    def prefix(a_list, pre):
        return [pre+'{0}'.format(an_element) for an_element in a_list]
    def postfix(a_list, post):
        return ['{0}'.format(an_element)+post for an_element in a_list]
    def surround(a_list, pre, post):
        return [pre+'{0}'.format(an_element)+post for an_element in a_list]
    env.filters['prefix'] = prefix
    env.filters['postfix'] = postfix
    env.filters['surround'] = surround

    erfa_c_in = env.get_template('core.c.templ')
    erfa_py_in = env.get_template('core.py.templ')

    #Extract all the ERFA function names from erfa.h
    if os.path.isdir(srcdir):
        erfahfn = os.path.join(srcdir, 'erfa.h')
        multifilserc = True
    else:
        erfahfn = os.path.join(os.path.split(srcdir)[0], 'erfa.h')
        multifilserc = False

    with open(erfahfn, "r") as f:
        erfa_h = f.read()

    funcs = OrderedDict()
    section_subsection_functions = re.findall('/\* (\w*)/(\w*) \*/\n(.*?)\n\n',
                                              erfa_h, flags=re.DOTALL|re.MULTILINE)
    for section, subsection, functions in section_subsection_functions:
        print_("{0}.{1}".format(section, subsection))
        if section == "Astronomy":
            func_names = re.findall(' (\w+)\(.*?\);', functions, flags=re.DOTALL)
            for name in func_names:
                print_("{0}.{1}.{2}...".format(section, subsection, name))
                if multifilserc:
                    # easy because it just looks in the file itself
                    funcs[name] = Function(name, srcdir)
                else:
                    # Have to tell it to look for a declaration matching
                    # the start of the header declaration, otherwise it
                    # might find a *call* of the function instead of the
                    # definition
                    for line in functions.split('\n'):
                        if name in line:
                            # [:-1] is to remove trailing semicolon, and
                            # splitting on '(' is because the header and
                            # C files don't necessarily have to match
                            # argument names and line-breaking or
                            # whitespace
                            match_line = line[:-1].split('(')[0]
                            funcs[name] = Function(name, srcdir, match_line)
                            break
                    else:
                        raise ValueError("A name for a C file wasn't "
                                         "found in the string that "
                                         "spawned it.  This should be "
                                         "impossible!")

    funcs = list(funcs.values())

    #Extract all the ERFA constants from erfam.h
    erfamhfn = os.path.join(srcdir, 'erfam.h')
    with open(erfamhfn, 'r') as f:
        erfa_m_h = f.read()
    constants = []
    for chunk in erfa_m_h.split("\n\n"):
        result = re.findall("#define (ERFA_\w+?) (.+?)$", chunk, flags=re.DOTALL|re.MULTILINE)
        if result:
            doc = re.findall("/\* (.+?) \*/\n", chunk, flags=re.DOTALL)
            for (name, value) in result:
                constants.append(Constant(name, value, doc))

    print_("Rendering template")
    erfa_c = erfa_c_in.render(funcs=funcs)
    erfa_py = erfa_py_in.render(funcs=funcs, constants=constants)

    if outfn is not None:
        outfn_c = os.path.splitext(outfn)[0] + ".c"
        print_("Saving to", outfn, 'and', outfn_c)
        with open(outfn, "w") as f:
            f.write(erfa_py)
        with open(outfn_c, "w") as f:
            f.write(erfa_c)

    print_("Done!")

    return erfa_c, erfa_py, funcs

DEFAULT_ERFA_LOC = os.path.join(os.path.split(__file__)[0],
                                '../../cextern/erfa')
DEFAULT_TEMPLATE_LOC = os.path.split(__file__)[0]

if __name__ == '__main__':
    from argparse import ArgumentParser

    ap = ArgumentParser()
    ap.add_argument('srcdir', default=DEFAULT_ERFA_LOC, nargs='?',
                    help='Directory where the ERFA c and header files '
                         'can be found or to a single erfa.c file '
                         '(which must be in the same directory as '
                         'erfa.h). Defaults to the builtin astropy '
                         'erfa: "{0}"'.format(DEFAULT_ERFA_LOC))
    ap.add_argument('-o', '--output', default='core.py',
                    help='The output filename.  This is the name for only the '
                         'pure-python output, the C part will have the '
                         'same name but with a ".c" extension.')
    ap.add_argument('-t', '--template-loc',
                    default=DEFAULT_TEMPLATE_LOC,
                    help='the location where the "core.c.templ" '
                         'template can be found.')
    ap.add_argument('-q', '--quiet', action='store_false', dest='verbose',
                    help='Suppress output normally printed to stdout.')

    args = ap.parse_args()
    main(args.srcdir, args.output, args.template_loc)
