from jinja2 import Environment, PackageLoader
from cython_numpy_auto.code_analyzer import Function

env = Environment(loader=PackageLoader('cython_numpy_auto', '.'))
erfa_pyx_in = env.get_template('erfa.pyx.in')

era_func_names = ["eraAtco13", "eraD2dtf", "eraAper"]
source_path = "../../erfa/src"

erfa_pyx = erfa_pyx_in.render(funcs=[Function(name, source_path) for name in era_func_names])

with open("erfa.pyx", "w") as f:
    f.write(erfa_pyx)
