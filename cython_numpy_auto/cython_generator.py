from jinja2 import Environment, PackageLoader
from cython_numpy_auto.code_analyzer import Function

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

era_func_names = ["eraAtco13", "eraD2dtf", "eraAper"]
source_path = "../../erfa/src"

erfa_pyx = erfa_pyx_in.render(funcs=[Function(name, source_path) for name in era_func_names])

with open("erfa.pyx", "w") as f:
    f.write(erfa_pyx)
