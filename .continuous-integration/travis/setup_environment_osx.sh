#!/bin/bash

# Install conda
# http://conda.pydata.org/docs/travis.html#the-travis-yml-file
wget http://repo.continuum.io/miniconda/Miniconda-3.7.3-MacOSX-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"

# Install Python dependencies
source "$( dirname "${BASH_SOURCE[0]}" )"/setup_dependencies_common.sh
