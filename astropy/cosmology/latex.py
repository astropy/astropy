import re                 # regex operations
import os                 # file removal
import tempfile           # creating temporary files
from astropy.cosmology import FlatLambdaCDM  # the cosmology class we build

def read_latex(file_path):
    # Open the .tex file and read all lines
    f = open(file_path, 'r')
    lines = f.readlines()
    f.close()

    params = {}  # store parsed parameters

    # regex: (name) = (number)
    pattern = re.compile(r"([A-Za-z0-9_]+)\s*=\s*([0-9.+\-Ee]+)")
    for line in lines:
        m = pattern.search(line)
        if m:
            key = m.group(1)           # parameter name
            val = float(m.group(2))    # numeric value
            params[key] = val

    # default values if missing
    H0_val = params.get('H0', 70.0)
    Om0_val = params.get('Om0', 0.3)

    # build and return a FlatLambdaCDM object
    cosm = FlatLambdaCDM(H0=H0_val, Om0=Om0_val)
    return cosm

def write_latex(cosmology, file_path):
    # prepare lines for a minimal .tex document
    lines = []
    h0 = cosmology.H0.value     # extract H0
    om0 = cosmology.Om0         # extract Om0
    ow0 = 1.0 - om0             # compute Ode0

    lines.append("\\documentclass{article}\n")
    lines.append("\\begin{document}\n")
    lines.append("\\section*{Cosmology Parameters}\n")
    lines.append(f"H0 = {h0:.2f}\n")      # write H0
    lines.append(f"Om0 = {om0:.4f}\n")    # write Om0
    lines.append(f"Ode0 = {ow0:.4f}\n")   # write Ode0
    lines.append("\\end{document}\n")

    # write to disk
    f = open(file_path, 'w')
    f.writelines(lines)
    f.close()

def _create_temp_latex(content):
    # helper: write raw string to a temp .tex file
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix='.tex', mode='w')
    tmp.write(content)
    tmp.close()
    return tmp.name

def test_read_latex():
    # make a small LaTeX document, read it, then assert values
    sample = (
        "\\documentclass{article}\n"
        "\\begin{document}\n"
        "H0 = 65.0\n"
        "Om0 = 0.25\n"
        "\\end{document}\n"
    )
    path = _create_temp_latex(sample)
    cosmo = read_latex(path)
    os.remove(path)
    assert abs(cosmo.H0.value - 65.0) < 1e-6
    assert abs(cosmo.Om0 - 0.25) < 1e-6

def test_write_and_read_roundtrip():
    # write a FlatLambdaCDM to .tex and read it back
    original = FlatLambdaCDM(H0=68.0, Om0=0.32)
    tmpfile = tempfile.NamedTemporaryFile(delete=False, suffix='.tex')
    tmpfile.close()

    write_latex(original, tmpfile.name)  # dump to file
    loaded = read_latex(tmpfile.name)    # parse it back
    os.remove(tmpfile.name)

    assert abs(loaded.H0.value - 68.0) < 1e-6
    assert abs(loaded.Om0 - 0.32) < 1e-6

def read_latex_from_string(latex_string):
    # convenience: parse a raw LaTeX string without manual files
    tmp_path = _create_temp_latex(latex_string)
    result = read_latex(tmp_path)
    os.remove(tmp_path)
    return result

def test_read_latex_from_string_function():
    # test the string‐based reader
    latex_str = "H0 = 72.0\nOm0 = 0.28\n"
    cosmo = read_latex_from_string(latex_str)
    assert abs(cosmo.H0.value - 72.0) < 1e-6
    assert abs(cosmo.Om0 - 0.28) < 1e-6

def example_usage():
    # example script you can run directly
    content = (
        "\\documentclass{article}\n"
        "\\begin{document}\n"
        "H0 = 75.0\n"
        "Om0 = 0.27\n"
        "Ode0 = 0.73\n"
        "\\end{document}\n"
    )
    temp_path = _create_temp_latex(content)
    cosmo = read_latex(temp_path)
    os.remove(temp_path)
    print("Loaded cosmology:")
    print("  H0 =", cosmo.H0.value)
    print("  Om0 =", cosmo.Om0)

if __name__ == "__main__":
    example_usage()
