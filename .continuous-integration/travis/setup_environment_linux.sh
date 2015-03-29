#!/bin/bash

# Install conda
wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b
export PATH=/home/travis/miniconda/bin:$PATH
conda update --yes conda

# Install non-Python dependencies for documentation
if [[ $SETUP_CMD == build_sphinx* ]]
then
  sudo apt-get update
  sudo apt-get install graphviz texlive-latex-extra dvipng
fi

# Install Python dependencies
source "$( dirname "${BASH_SOURCE[0]}" )"/setup_dependencies_common.sh
