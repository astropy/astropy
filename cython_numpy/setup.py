from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy as np

setup(
    name = "erfa",
    ext_modules = cythonize([
        Extension(name="erfa",
                  sources=["erfa.pyx"],
                  libraries=['erfa'],
                  include_dirs = [np.get_include(), '/opt/local/include'],
                  library_dirs = ['/opt/local/lib'],
                  language="c",)
    ])
)
