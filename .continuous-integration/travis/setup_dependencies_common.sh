#!/bin/bash -x

hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda info -a

# CONDA
conda create --yes -n test -c astropy-ci-extras python=$PYTHON_VERSION pip
source activate test

# --no-use-wheel requirement is temporary due to
# https://github.com/astropy/astropy/issues/4180
# and may be removed once the upstream fix is in place
export PIP_INSTALL="pip install --no-use-wheel"

# EGG_INFO
if [[ $SETUP_CMD == egg_info ]]
then
  return  # no more dependencies needed
fi

# PEP8
if [[ $MAIN_CMD == pep8* ]]
then
  $PIP_INSTALL pep8
  return  # no more dependencies needed
fi

# CORE DEPENDENCIES
conda install --yes pytest Cython jinja2 psutil

# NUMPY
if [[ $NUMPY_VERSION == dev ]]
then
  $PIP_INSTALL git+http://github.com/numpy/numpy.git
  export CONDA_INSTALL="conda install --yes python=$PYTHON_VERSION"
else
  conda install --yes numpy=$NUMPY_VERSION
  export CONDA_INSTALL="conda install --yes python=$PYTHON_VERSION numpy=$NUMPY_VERSION"
fi

# Now set up shortcut to conda install command to make sure the Python and Numpy
# versions are always explicitly specified.

# OPTIONAL DEPENDENCIES
if $OPTIONAL_DEPS
then
  $CONDA_INSTALL scipy h5py matplotlib pyyaml scikit-image pandas
  $PIP_INSTALL beautifulsoup4
fi

# DOCUMENTATION DEPENDENCIES
# build_sphinx needs sphinx as well as matplotlib and wcsaxes (for plot_directive). 
# Note that this matplotlib will *not* work with py 3.x, but our sphinx build is
# currently 2.7, so that's fine
if [[ $SETUP_CMD == build_sphinx* ]]
then
  $CONDA_INSTALL Sphinx=1.2.2 Pygments matplotlib
  $PIP_INSTALL wcsaxes
fi

# COVERAGE DEPENDENCIES
# cpp-coveralls must be installed first.  It installs two identical
# scripts: 'cpp-coveralls' and 'coveralls'.  The latter will overwrite
# the script installed by 'coveralls', unless it's installed first.
if [[ $SETUP_CMD == 'test --coverage' ]]
then
  $PIP_INSTALL cpp-coveralls;
  $PIP_INSTALL coverage coveralls;
fi
