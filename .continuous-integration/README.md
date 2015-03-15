About
=====

This directory contains a set of scripts that are used by the 
``.travis.yml`` and ``appveyor.yml`` files for the 
[Travis](http://travis-ci.org) and [AppVeyor](http://www.appveyor.com/) 
services respectively.

The scripts include:

* ``appveyor/install-miniconda.ps1`` - set up conda on Windows
* ``appveyor/windows_sdk.cmd`` - set up the compiler environment on Windows
* ``travis/setup_dependencies_common.sh`` - set up conda packages on Linux and MacOS X
* ``travis/setup_environment_linux.sh`` - set up conda and non-Python dependencies on Linux
* ``travis/setup_environment_osx.sh`` - set up conda and non-Python dependencies on MacOS X
