#!/bin/bash

docker info

cat << EOF | docker run -i \
                        -v ${PWD}:/astropy_src \
                        -a stdin -a stdout -a stderr \
                        astropy/astropy-32bit-test-env:1.9 \
                        bash || exit $?

cd /astropy_src

echo "Output of uname -m:"
uname -m

echo "Output of sys.maxsize in Python:"
python -c 'import sys; print(sys.maxsize)'

easy_install pytest pytest-xdist

python setup.py test --parallel=4

EOF
