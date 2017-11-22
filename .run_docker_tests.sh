#!/bin/bash

docker info

cat << EOF | docker run -i \
                        -v ${PWD}:/astropy_src \
                        -a stdin -a stdout -a stderr \
                        astropy/astropy-py35-32bit-test-env:1.10 \
                        bash || exit $?

cd /astropy_src

echo "Output of uname -m:"
uname -m

echo "Output of sys.maxsize in Python:"
python3 -c 'import sys; print(sys.maxsize)'

easy_install-3.5 pytest-xdist

PYTHONHASHSEED=42 python3 setup.py test --parallel=4

EOF
