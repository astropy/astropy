# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module's main purpose is to act as a script to create new versions
of ufunc.c when ERFA is updated (or this generator is enhanced).

`Jinja2 <http://jinja.pocoo.org/>`_ must be installed for this
module/script to function.

Note that this does *not* currently automate the process of creating structs
or dtypes for those structs.  They should be added manually in the template file.
"""

import re
import os.path
from collections import OrderedDict

DEFAULT_ERFA_LOC = os.path.join(os.path.split(__file__)[0],
                                '../../cextern/erfa')
DEFAULT_TEMPLATE_LOC = os.path.split(__file__)[0]

NDIMS_REX = re.compile(re.escape("numpy.dtype([('fi0', '.*', <(.*)>)])").replace(r'\.\*', '.*').replace(r'\<', '(').replace(r'\>', ')'))


class FunctionDoc:

    def __init__(self, doc):
        self.doc = doc.replace("**", "  ").replace("/*\n", "").replace("*/", "")
        self.__input = None
        self.__output = None
        self.__ret_info = None

    def _get_arg_doc_list(self, doc_lines):
        """Parse input/output doc section lines, getting arguments from them.

        Ensure all elements of eraASTROM and eraLDBODY are left out, as those
        are not input or output arguments themselves.  Also remove the nb
        argument in from of eraLDBODY, as we infer nb from the python array.
        """
        doc_list = []
        skip = []
        for d in doc_lines:
            arg_doc = ArgumentDoc(d)
            if arg_doc.name is not None:
                if skip:
                    if skip[0] == arg_doc.name:
                        skip.pop(0)
                        continue
                    else:
                        raise RuntimeError("We whould be skipping {} "
                                           "but {} encountered."
                                           .format(skip[0], arg_doc.name))

                if arg_doc.type.startswith('eraLDBODY'):
                    # Special-case LDBODY: for those, the previous argument
                    # is always the number of bodies, but we don't need it
                    # as an input argument for the ufunc since we're going
                    # to determine this from the array itself. Also skip
                    # the description of its contents; those are not arguments.
                    doc_list.pop()
                    skip = ['bm', 'dl', 'pv']
                elif arg_doc.type.startswith('eraASTROM'):
                    # Special-case ASTROM: need to skip the description
                    # of its contents; those are not arguments.
                    skip = ['pmt', 'eb', 'eh', 'em', 'v', 'bm1',
                            'bpn', 'along', 'xpl', 'ypl', 'sphi',
                            'cphi', 'diurab', 'eral', 'refa', 'refb']

                doc_list.append(arg_doc)

        return doc_list

    @property
    def input(self):
        if self.__input is None:
            self.__input = []
            for regex in ("Given([^\n]*):.*?\n(.+?)  \n",
                          "Given and returned([^\n]*):\n(.+?)  \n"):
                result = re.search(regex, self.doc, re.DOTALL)
                if result is not None:
                    doc_lines = result.group(2).split("\n")
                    self.__input += self._get_arg_doc_list(doc_lines)

        return self.__input

    @property
    def output(self):
        if self.__output is None:
            self.__output = []
            for regex in ("Given and returned([^\n]*):\n(.+?)  \n",
                          "Returned([^\n]*):.*?\n(.+?)  \n"):
                result = re.search(regex, self.doc, re.DOTALL)
                if result is not None:
                    doc_lines = result.group(2).split("\n")
                    self.__output += self._get_arg_doc_list(doc_lines)

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


class ArgumentDoc:

    def __init__(self, doc):
        match = re.search("^ +([^ ]+)[ ]+([^ ]+)[ ]+(.+)", doc)
        if match is not None:
            self.name = match.group(1)
            if self.name.startswith('*'):  # easier than getting the regex to behave...
                self.name = self.name.replace('*', '')
            self.type = match.group(2)
            self.doc = match.group(3)
        else:
            self.name = None
            self.type = None
            self.doc = None

    def __repr__(self):
        return f"    {self.name:15} {self.type:15} {self.doc}"


class Variable:
    """Properties shared by Argument and Return."""
    @property
    def npy_type(self):
        """Predefined type used by numpy ufuncs to indicate a given ctype.

        Eg., NPY_DOUBLE for double.
        """
        return "NPY_" + self.ctype.upper()

    @property
    def dtype(self):
        """Name of dtype corresponding to the ctype.

        Specifically,
        double : dt_double
        int : dt_int
        double[3]: dt_vector
        double[2][3] : dt_pv
        double[2] : dt_pvdpv
        double[3][3] : dt_matrix
        int[4] : dt_ymdf | dt_hmsf | dt_dmsf, depding on name
        eraASTROM: dt_eraASTROM
        eraLDBODY: dt_eraLDBODY
        char : dt_sign
        char[] : dt_type

        The corresponding dtypes are defined in ufunc.c, where they are
        used for the loop definitions.  In core.py, they are also used
        to view-cast regular arrays to these structured dtypes.
        """
        if self.ctype == 'const char':
            return 'dt_type'
        elif self.ctype == 'char':
            return 'dt_sign'
        elif self.ctype == 'int' and self.shape == (4,):
            return 'dt_' + self.name[1:]
        elif self.ctype == 'double' and self.shape == (3,):
            return 'dt_double'
        elif self.ctype == 'double' and self.shape == (2, 3):
            return 'dt_pv'
        elif self.ctype == 'double' and self.shape == (2,):
            return 'dt_pvdpv'
        elif self.ctype == 'double' and self.shape == (3, 3):
            return 'dt_double'
        elif not self.shape:
            return 'dt_' + self.ctype
        else:
            raise ValueError("ctype {} with shape {} not recognized."
                             .format(self.ctype, self.shape))

    @property
    def view_dtype(self):
        """Name of dtype corresponding to the ctype for viewing back as array.

        E.g., dt_double for double, dt_double33 for double[3][3].

        The types are defined in core.py, where they are used for view-casts
        of structured results as regular arrays.
        """
        if self.ctype == 'const char':
            return 'dt_bytes12'
        elif self.ctype == 'char':
            return 'dt_bytes1'
        else:
            raise ValueError('Only char ctype should need view back!')

    @property
    def ndim(self):
        return len(self.shape)

    @property
    def size(self):
        size = 1
        for s in self.shape:
            size *= s
        return size

    @property
    def cshape(self):
        return ''.join([f'[{s}]' for s in self.shape])

    @property
    def signature_shape(self):
        if self.ctype == 'eraLDBODY':
            return '(n)'
        elif self.ctype == 'double' and self.shape == (3,):
            return '(3)'
        elif self.ctype == 'double' and self.shape == (3, 3):
            return '(3, 3)'
        else:
            return '()'


class Argument(Variable):

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
    def name_for_call(self):
        """How the argument should be used in the call to the ERFA function.

        This takes care of ensuring that inputs are passed by value,
        as well as adding back the number of bodies for any LDBODY argument.
        The latter presumes that in the ufunc inner loops, that number is
        called 'nb'.
        """
        if self.ctype == 'eraLDBODY':
            assert self.name == 'b'
            return 'nb, _' + self.name
        elif self.is_ptr:
            return '_'+self.name
        else:
            return '*_'+self.name

    def __repr__(self):
        return f"Argument('{self.definition}', name='{self.name}', ctype='{self.ctype}', inout_state='{self.inout_state}')"


class ReturnDoc:

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
        return f"Return value, type={self.type:15}, {self.descr}, {self.doc}"


class Return(Variable):

    def __init__(self, ctype, doc):
        self.name = 'c_retval'
        self.inout_state = 'stat' if ctype == 'int' else 'ret'
        self.ctype = ctype
        self.shape = ()
        self.doc = doc

    def __repr__(self):
        return f"Return(name='{self.name}', ctype='{self.ctype}', inout_state='{self.inout_state}')"

    @property
    def doc_info(self):
        return self.doc.ret_info


class Function:
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

        pattern = fr"\n([^\n]+{name} ?\([^)]+\)).+?(/\*.+?\*/)"
        p = re.compile(pattern, flags=re.DOTALL | re.MULTILINE)

        search = p.search(filecontents)
        self.cfunc = " ".join(search.group(1).split())
        self.doc = FunctionDoc(search.group(2))

        self.args = []
        for arg in re.search(r"\(([^)]+)\)", self.cfunc).group(1).split(', '):
            self.args.append(Argument(arg, self.doc))
        self.ret = re.search(f"^(.*){name}", self.cfunc).group(1).strip()
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

    @property
    def user_dtype(self):
        """The non-standard dtype, if any, needed by this function's ufunc.

        This would be any structured array for any input or output, but
        we give preference to LDBODY, since that also decides that the ufunc
        should be a generalized ufunc.
        """
        user_dtype = None
        for arg in self.args_by_inout('in|inout|out'):
            if arg.ctype == 'eraLDBODY':
                return arg.dtype
            elif user_dtype is None and arg.dtype not in ('dt_double',
                                                          'dt_int'):
                user_dtype = arg.dtype

        return user_dtype

    @property
    def signature(self):
        """Possible signature, if this function should be a gufunc."""
        if all(arg.signature_shape == '()'
               for arg in self.args_by_inout('in|inout|out')):
            return None

        return '->'.join(
            [','.join([arg.signature_shape for arg in args])
             for args in (self.args_by_inout('in|inout'),
                          self.args_by_inout('inout|out|ret|stat'))])

    @property
    def python_call(self):
        outnames = [arg.name for arg in self.args_by_inout('inout|out|stat|ret')]
        argnames = [arg.name for arg in self.args_by_inout('in|inout')]
        return '{out} = {func}({args})'.format(out=', '.join(outnames),
                                               func='ufunc.' + self.pyname,
                                               args=', '.join(argnames))

    def __repr__(self):
        return f"Function(name='{self.name}', pyname='{self.pyname}', filename='{self.filename}', filepath='{self.filepath}')"


class Constant:

    def __init__(self, name, value, doc):
        self.name = name.replace("ERFA_", "")
        self.value = value.replace("ERFA_", "")
        self.doc = doc


class ExtraFunction(Function):
    """
    An "extra" function - e.g. one not following the SOFA/ERFA standard format.

    Parameters
    ----------
    cname : str
        The name of the function in C
    prototype : str
        The prototype for the function (usually derived from the header)
    pathfordoc : str
        The path to a file that contains the prototype, with the documentation
        as a multiline string *before* it.
    """

    def __init__(self, cname, prototype, pathfordoc):
        self.name = cname
        self.pyname = cname.split('era')[-1].lower()
        self.filepath, self.filename = os.path.split(pathfordoc)

        self.prototype = prototype.strip()
        if prototype.endswith('{') or prototype.endswith(';'):
            self.prototype = prototype[:-1].strip()

        incomment = False
        lastcomment = None
        with open(pathfordoc, 'r') as f:
            for l in f:
                if incomment:
                    if l.lstrip().startswith('*/'):
                        incomment = False
                        lastcomment = ''.join(lastcomment)
                    else:
                        if l.startswith('**'):
                            l = l[2:]
                        lastcomment.append(l)
                else:
                    if l.lstrip().startswith('/*'):
                        incomment = True
                        lastcomment = []
                    if l.startswith(self.prototype):
                        self.doc = lastcomment
                        break
            else:
                raise ValueError('Did not find prototype {} in file '
                                 '{}'.format(self.prototype, pathfordoc))

        self.args = []
        argset = re.search(fr"{self.name}\(([^)]+)?\)",
                           self.prototype).group(1)
        if argset is not None:
            for arg in argset.split(', '):
                self.args.append(Argument(arg, self.doc))
        self.ret = re.match(f"^(.*){self.name}",
                            self.prototype).group(1).strip()
        if self.ret != 'void':
            self.args.append(Return(self.ret, self.doc))

    def __repr__(self):
        r = super().__repr__()
        if r.startswith('Function'):
            r = 'Extra' + r
        return r


def main(srcdir=DEFAULT_ERFA_LOC, outfn='core.py', ufuncfn='ufunc.c',
         templateloc=DEFAULT_TEMPLATE_LOC, extra='erfa_additions.h',
         verbose=True):
    from jinja2 import Environment, FileSystemLoader

    if verbose:
        print_ = lambda *args, **kwargs: print(*args, **kwargs)
    else:
        print_ = lambda *args, **kwargs: None

    # Prepare the jinja2 templating environment
    env = Environment(loader=FileSystemLoader(templateloc))

    def prefix(a_list, pre):
        return [pre+f'{an_element}' for an_element in a_list]

    def postfix(a_list, post):
        return [f'{an_element}'+post for an_element in a_list]

    def surround(a_list, pre, post):
        return [pre+f'{an_element}'+post for an_element in a_list]
    env.filters['prefix'] = prefix
    env.filters['postfix'] = postfix
    env.filters['surround'] = surround

    erfa_c_in = env.get_template(ufuncfn + '.templ')
    erfa_py_in = env.get_template(outfn + '.templ')

    # Extract all the ERFA function names from erfa.h
    if os.path.isdir(srcdir):
        erfahfn = os.path.join(srcdir, 'erfa.h')
        multifilserc = True
    else:
        erfahfn = os.path.join(os.path.split(srcdir)[0], 'erfa.h')
        multifilserc = False

    with open(erfahfn, "r") as f:
        erfa_h = f.read()
        print_("read erfa header")
    if extra:
        with open(os.path.join(templateloc or '.', extra), "r") as f:
            erfa_h += f.read()
        print_("read extra header")

    funcs = OrderedDict()
    section_subsection_functions = re.findall(
        r'/\* (\w*)/(\w*) \*/\n(.*?)\n\n', erfa_h,
        flags=re.DOTALL | re.MULTILINE)
    for section, subsection, functions in section_subsection_functions:
        print_(f"{section}.{subsection}")
        # Right now, we compile everything, but one could be more selective.
        # In particular, at the time of writing (2018-06-11), what was
        # actually require for astropy was not quite everything, but:
        # ((section == 'Extra')
        #       or (section == "Astronomy")
        #       or (subsection == "AngleOps")
        #       or (subsection == "SphericalCartesian")
        #       or (subsection == "MatrixVectorProducts")
        #       or (subsection == 'VectorOps'))
        if True:

            func_names = re.findall(r' (\w+)\(.*?\);', functions,
                                    flags=re.DOTALL)
            for name in func_names:
                print_(f"{section}.{subsection}.{name}...")
                if multifilserc:
                    # easy because it just looks in the file itself
                    cdir = (srcdir if section != 'Extra' else
                            templateloc or '.')
                    funcs[name] = Function(name, cdir)
                else:
                    # Have to tell it to look for a declaration matching
                    # the start of the header declaration, otherwise it
                    # might find a *call* of the function instead of the
                    # definition
                    for line in functions.split(r'\n'):
                        if name in line:
                            # [:-1] is to remove trailing semicolon, and
                            # splitting on '(' is because the header and
                            # C files don't necessarily have to match
                            # argument names and line-breaking or
                            # whitespace
                            match_line = line[:-1].split('(')[0]
                            funcs[name] = Function(name, cdir, match_line)
                            break
                    else:
                        raise ValueError("A name for a C file wasn't "
                                         "found in the string that "
                                         "spawned it.  This should be "
                                         "impossible!")

    funcs = funcs.values()

    # Extract all the ERFA constants from erfam.h
    erfamhfn = os.path.join(srcdir, 'erfam.h')
    with open(erfamhfn, 'r') as f:
        erfa_m_h = f.read()
    constants = []
    for chunk in erfa_m_h.split("\n\n"):
        result = re.findall(r"#define (ERFA_\w+?) (.+?)$", chunk,
                            flags=re.DOTALL | re.MULTILINE)
        if result:
            doc = re.findall(r"/\* (.+?) \*/\n", chunk, flags=re.DOTALL)
            for (name, value) in result:
                constants.append(Constant(name, value, doc))

    # TODO: re-enable this when const char* return values and
    #       non-status code integer rets are possible
    # #Add in any "extra" functions from erfaextra.h
    # erfaextrahfn = os.path.join(srcdir, 'erfaextra.h')
    # with open(erfaextrahfn, 'r') as f:
    #     for l in f:
    #         ls = l.strip()
    #         match = re.match('.* (era.*)\(', ls)
    #         if match:
    #             print_("Extra:  {0} ...".format(match.group(1)))
    #             funcs.append(ExtraFunction(match.group(1), ls, erfaextrahfn))

    print_("Rendering template")
    erfa_c = erfa_c_in.render(funcs=funcs)
    erfa_py = erfa_py_in.render(funcs=funcs, constants=constants)

    if outfn is not None:
        print_("Saving to", outfn, 'and', ufuncfn)
        with open(os.path.join(templateloc, outfn), "w") as f:
            f.write(erfa_py)
        with open(os.path.join(templateloc, ufuncfn), "w") as f:
            f.write(erfa_c)

    print_("Done!")

    return erfa_c, erfa_py, funcs


if __name__ == '__main__':
    from argparse import ArgumentParser

    ap = ArgumentParser()
    ap.add_argument('srcdir', default=DEFAULT_ERFA_LOC, nargs='?',
                    help='Directory where the ERFA c and header files '
                         'can be found or to a single erfa.c file '
                         '(which must be in the same directory as '
                         'erfa.h). Defaults to the builtin astropy '
                         'erfa: "{}"'.format(DEFAULT_ERFA_LOC))
    ap.add_argument('-o', '--output', default='core.py',
                    help='The output filename for the pure-python output.')
    ap.add_argument('-u', '--ufunc', default='ufunc.c',
                    help='The output filename for the ufunc .c output')
    ap.add_argument('-t', '--template-loc',
                    default=DEFAULT_TEMPLATE_LOC,
                    help='the location where the "core.py.templ" and '
                         '"ufunc.c.templ templates can be found.')
    ap.add_argument('-x', '--extra',
                    default='erfa_additions.h',
                    help='header file for any extra files in the template '
                         'location that should be included.')
    ap.add_argument('-q', '--quiet', action='store_false', dest='verbose',
                    help='Suppress output normally printed to stdout.')

    args = ap.parse_args()
    main(args.srcdir, args.output, args.ufunc, args.template_loc,
         args.extra)
