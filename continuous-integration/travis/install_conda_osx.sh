#!/bin/bash

wget http://repo.continuum.io/miniconda/Miniconda-3.7.3-MacOSX-x86_64.sh -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b
export PATH=/Users/travis/miniconda/bin:$PATH
conda update --yes conda
