import re
from jinja2 import Environment, PackageLoader
from cython_numpy_auto.code_analyzer import Function

ERFA_SOURCES = "../../erfa/src"

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
with open(ERFA_SOURCES+"/erfa.h", "r") as f:
    func_names = re.findall(' (\w+)\(.*?\);', f.read(), flags=re.DOTALL)
    funcs = []
    for name in func_names:
        print("Parsing {0}...".format(name))
        funcs.append(Function(name, ERFA_SOURCES))
        print("Done!")
    #Render the template and save
    erfa_pyx = erfa_pyx_in.render(funcs=funcs)
    with open("erfa.pyx", "w") as f:
        f.write(erfa_pyx)
