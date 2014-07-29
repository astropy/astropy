from code_analyzer import *


def generate(era_func_names, source_path):

    era_funcs = [Function(name, source_path) for name in era_func_names]
    
    wrapper = []
    wrapper.append("import numpy as np")
    wrapper.append("cimport numpy as np")
    wrapper.append("")
    wrapper.append("np.import_array()")
    wrapper.append("")
    wrapper.append("__all__ = [{0}]".format(", ".join(["'{0}'".format(func.pyname) for func in era_funcs]+["'dt_eraASTROM'"])))
    wrapper.append("")
    wrapper.append("cdef extern from 'erfa.h':")
    wrapper.append("    struct eraASTROM:")
    wrapper.append("        pass")
    for func in era_funcs:
        wrapper.append("    {0}".format(func.cfunc))
    wrapper.append("")
    wrapper.append("dt_eraASTROM = np.dtype([('pmt','d'),('eb','d',(3,)),('eh','d',(3,)),('em','d'),('v','d',(3,)),('bm1 ','d'),('bpn','d',(3,3)),('along','d'),('phi','d'),('xpl','d'),('ypl','d'),('sphi','d'),('cphi','d'),('diurab','d'),('eral','d'),('refa','d'),('refb','d')], align=True)")
    wrapper.append("")
    for func in era_funcs:
        wrapper.append("#=== {0} ===".format(func.name))
        wrapper.append("")
        wrapper.append("def {0}({1}):".format(func.pyname, func.args_by_inout('in|inout','name',', ')))
        wrapper.append("    ")
        wrapper.append("    shape = np.broadcast({0}).shape".format(func.args_by_inout('in|inout','name',', ')))
        for arg in func.args_by_inout('out'):
            wrapper.append("    {0}_out = np.empty(shape, dtype={1})".format(arg.name, arg.dtype))
        for arg in func.args_by_inout('inout'):
            wrapper.append("    {0}_out = np.array(np.broadcast_arrays({1})[{2}], dtype={3})".format(arg.name, func.args_by_inout('in|inout','name',', '), func.args.index(arg), arg.dtype))
        wrapper.append("    ")
        wrapper.append("    cdef np.broadcast it = np.broadcast({0})".format(', '.join([arg.name if arg.is_in else (arg.name+"_out") for arg in func.args_by_inout('in|inout|out')])))
        wrapper.append("    ")
        for arg in func.args_by_inout('in|out|inout'):
            wrapper.append("    cdef {0} _{1}".format(arg.ctype_ptr, arg.name))
        wrapper.append("    ")
        wrapper.append("    while np.PyArray_MultiIter_NOTDONE(it):")
        wrapper.append("        ")
        for arg in func.args_by_inout('in|out|inout'):
            if arg.ctype_ptr[-1] == '*':
                wrapper.append("        _{0} = (<{1}>np.PyArray_MultiIter_DATA(it, {2}))".format(arg.name, arg.ctype_ptr, func.args.index(arg)))
            else:
                wrapper.append("        _{0} = (<{1}*>np.PyArray_MultiIter_DATA(it, {2}))[0]".format(arg.name, arg.ctype_ptr, func.args.index(arg)))
        wrapper.append("        ")
        wrapper.append("        {0}({1})".format(func.name, ", ".join(["_"+name for name in func.args_by_inout('in|out|inout','name')])))
        wrapper.append("        ")
        wrapper.append("        np.PyArray_MultiIter_NEXT(it)")
        wrapper.append("    ")
        wrapper.append("    return {0}".format(', '.join([arg+"_out" for arg in func.args_by_inout('out|inout','name')])))
        wrapper.append("")
    
    with open("./erfa.pyx","w") as f:
        f.write("\n".join(wrapper))


if __name__ == '__main__':
    era_func_names = ["eraAtco13", "eraD2dtf", "eraAper"]
    source_path = "../../erfa/src"
    generate(era_func_names, source_path)
