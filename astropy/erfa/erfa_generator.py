#!/usr/bin/env python

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module's main purpose is to act as a script to create new versions
of core.pyx when ERFA is updated (or this generator is enhanced).

`Jinja2 <http://jinja.pocoo.org/>`_ must be installed for this
module/script to function.

Note that this does *not* currently automate the process of creating structs or
dtypes for those structs.  They should be added manually in the template file.
"""

from __future__ import absolute_import, division, print_function
# note that we do *not* use unicode_literals here, because that makes the
# generated code's strings have u'' in them on py 2.x

import logging
import json
import re
import os.path
import sys

from argparse import ArgumentParser

from jinja2 import Environment, FileSystemLoader


log = logging.getLogger(__name__)
log.addHandler(logging.StreamHandler(sys.stdout))


ctype_to_dtype = {
    'double'     : "numpy.double",
    'int'        : "numpy.intc",
    'eraASTROM'  : "dt_eraASTROM",
    'eraLDBODY'  : "dt_eraLDBODY",
    'char'       : "numpy.dtype('S16')",
    'const char' : "numpy.dtype('S16')",
}


DEFAULT_ERFA_LOC = os.path.join(os.path.dirname(__file__), os.pardir,
                                os.pardir, 'cextern', 'erfa')
DEFAULT_TEMPLATE_LOC = os.path.join(os.path.dirname(__file__))


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
        self.definition = definition
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

    def to_dict(self):
        """
        Represent this `Argument` as a `dict` for serialization purposes.
        """

        ret = {'name': self.name, 'ctype': self.ctype,
               'direction': self.inout_state}
        if self.shape:
            ret['shape'] = self.shape

        return ret

    def __repr__(self):
        return ("Argument('{0}', name='{1}', ctype='{2}', "
                "inout_state='{3}')").format(self.definition, self.name,
                                             self.ctype, self.inout_state)


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
        return "Return value, type={0:15}, {1}, {2}".format(self.type,
                                                            self.descr,
                                                            self.doc)


class Return(object):

    def __init__(self, ctype, doc):
        self.name = 'c_statval' if ctype == 'int' else 'c_retval'
        self.name_out_broadcast = self.name+"_out"
        self.inout_state = 'stat' if ctype == 'int' else 'ret'
        self.ctype = ctype
        self.ctype_ptr = ctype
        self.shape = ()
        self.doc = doc

    def __repr__(self):
        return "Return(name='{0}', ctype='{1}', inout_state='{2}')".format(
                self.name, self.ctype, self.inout_state)

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
            self.filepath = os.path.join(os.path.normpath(source_path),
                                         self.filename)
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

    def to_dict(self):
        """
        Convert this `Function` to a `dict` representation for serialization.

        Does not dump the function's pyname, as this will be used as the key
        to this function in the dictionary of all ERFA functions.
        """

        ret = {'name': self.name, 'doc': str(self.doc)}
        args = []
        for arg in self.args:
            if isinstance(arg, Return):
                ret['return'] = arg.ctype
                if arg.doc_info and arg.doc_info.statuscodes:
                    ret['statuscodes'] = arg.doc_info.statuscodes
                continue
            args.append(arg.to_dict())

        if args:
            ret['arguments'] = args

        return ret

    def __repr__(self):
        return ("Function(name='{0}', pyname='{1}', filename='{2}', "
                "filepath='{3}')").format(self.name, self.pyname,
                                          self.filename, self.filepath)


class Constant(object):

    def __init__(self, name, value, doc):
        self.name = name.replace("ERFA_","")
        self.value = value.replace("ERFA_","")
        self.doc = doc


def extract_erfa_functions(src):
    """
    Extract all ERFA function declarations from erfa.h and the associated
    source file(s).
    """

    if os.path.isdir(src):
        erfa_h_filename = os.path.join(src, 'erfa.h')
        multi_file = True
    else:
        erfa_h_filename = os.path.join(os.path.dirname(src), 'erfa.h')
        multi_file = False

    with open(erfa_h_filename) as f:
        erfa_h = f.read()

    funcs = []
    section_subsection_functions = re.findall('/\* (\w*)/(\w*) \*/\n(.*?)\n\n',
                                              erfa_h, flags=re.DOTALL|re.MULTILINE)

    for section, subsection, functions in section_subsection_functions:
        log.debug("{0}.{1}".format(section, subsection))
        if section != "Astronomy":
            # For now we omit all functions not in the "Astronomy" section of
            # erfa.h
            continue

        func_names = re.findall(' (\w+)\(.*?\);', functions, flags=re.DOTALL)
        for name in func_names:
            log.debug("{0}.{1}.{2}...".format(section, subsection, name))
            if multi_file:
                # easy because it just looks in the file itself
                funcs.append(Function(name, src))
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
                        funcs.append(Function(name, src, match_line))
                        break
                else:
                    raise ValueError(
                        "A name for a C file wasn't found in the string that "
                        "spawned it.  This should be impossible!")

    return funcs


def extract_erfa_constants(src):
    """Extract all ERFA constants from erfam.h."""

    if os.path.isdir(src):
        erfam_h_filename = os.path.join(src, 'erfam.h')
    else:
        erfam_h_filename = os.path.join(os.path.dirname(src), 'erfam.h')

    with open(erfam_h_filename) as f:
        erfam_h = f.read()

    constants = []
    for chunk in erfam_h.split("\n\n"):
        result = re.findall("#define (ERFA_\w+?) (.+?)$", chunk,
                            flags=re.DOTALL|re.MULTILINE)
        if result:
            doc = re.findall("/\* (.+?) \*/\n", chunk, flags=re.DOTALL)
            for name, value in result:
                constants.append(Constant(name, value, doc))

    return constants


def generate_erfa_pyx(env, funcs, stream=None):
    """
    Generate the core.pyx file from the given template and list of functions.
    """

    log.debug("Rendering template")

    def prefix(items, pre):
        templ = '{0}{{0}}'.format(pre)
        return map(templ.format, items)

    def postfix(items, post):
        templ = '{{0}}{0}'.format(post)
        return map(templ.format, items)

    def surround(items, pre, post):
        templ = '{0}{{0}}{1}'.format(pre, post)
        return map(templ.format, items)

    env.filters['prefix'] = prefix
    env.filters['postfix'] = postfix
    env.filters['surround'] = surround

    erfa_pyx_in = env.get_template('core.pyx.templ')

    if stream is None:
        return erfa_pyx_in.render(funcs=funcs)
    else:
        erfa_pyx_in.stream(funcs=funcs).dump(stream)


def generate_erfa_json(funcs, stream=None):
    erfa_json = {func.pyname: func.to_dict() for func in funcs}

    if stream is None:
        return json.dumps(erfa_json, indent=2)
    else:
        json.dump(erfa_json, stream, indent=2)


def generate_erfa_constants(env, constants, stream=None):
    templ = env.get_template('constants.py.templ')

    if stream is None:
        return templ.render(constants=constants)
    else:
        templ.stream(constants=constants).dump(stream)


def write_erfa_sources(src=DEFAULT_ERFA_LOC, templates=DEFAULT_TEMPLATE_LOC,
                       out='.', verbose=False):
    """
    Generate and write the core.pyx and erfa.json sources.
    """

    if verbose:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)

    funcs = extract_erfa_functions(src)
    constants = extract_erfa_constants(src)

    erfa_pyx_filename = os.path.join(out, 'core.pyx')
    erfa_json_filename = os.path.join(out, 'erfa.json')
    constants_py_filename = os.path.join(out, 'constants.py')

    #Prepare the jinja2 templating environment
    env = Environment(loader=FileSystemLoader(templates))

    generate_erfa_pyx(env, funcs, stream=open(erfa_pyx_filename, 'w'))
    generate_erfa_json(funcs, stream=open(erfa_json_filename, 'w'))
    generate_erfa_constants(env, constants,
                            stream=open(constants_py_filename, 'w'))


def main(argv):
    ap = ArgumentParser()
    ap.add_argument('src', default=DEFAULT_ERFA_LOC, nargs='?',
                    help='Directory where the ERFA c and header files '
                         'can be found or to a single erfa.c file '
                         '(which must be in the same directory as '
                         'erfa.h). Defaults to the builtin astropy '
                         'erfa: "{0}"'.format(DEFAULT_ERFA_LOC))
    ap.add_argument('-t', '--templates', default=DEFAULT_TEMPLATE_LOC,
                    help='Path to look for Jinja2 templates for core.pyx '
                         'and constants.py')
    ap.add_argument('-o', '--output', default='.',
                    help='Directory to which generated sources should be '
                         'output (default: .)')
    ap.add_argument('-q', '--quiet', action='store_false', dest='verbose',
                    help='Suppress output normally printed to stdout.')

    args = ap.parse_args(argv)

    write_erfa_sources(args.src, args.templates, args.output,
                       verbose=args.verbose)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
