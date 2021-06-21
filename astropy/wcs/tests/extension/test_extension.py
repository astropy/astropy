# Licensed under a 3-clause BSD style license - see LICENSE.rst


import os
import subprocess
import sys

import pytest

import astropy


def test_wcsapi_extension(tmpdir):
    # Test that we can build a simple C extension with the astropy.wcs C API

    build_dir = tmpdir.mkdir('build').strpath
    install_dir = tmpdir.mkdir('install').strpath
    record_file = os.path.join(build_dir, 'record.txt')

    setup_path = os.path.dirname(__file__)
    astropy_path = os.path.dirname(astropy.__path__[0])

    env = os.environ.copy()
    paths = [install_dir, astropy_path]
    if env.get('PYTHONPATH'):
        paths.append(env.get('PYTHONPATH'))
    env['PYTHONPATH'] = os.pathsep.join(paths)

    # Build the extension
    # This used to use subprocess.check_call, but on Python 3.4 there was
    # a mysterious Heisenbug causing this to fail with a non-zero exit code
    # *unless* the output is redirected.  This bug also did not occur in an
    # interactive session, so it likely had something to do with pytest's
    # output capture
    #
    # --single-version-externally-managed --record is essentially equivalent to
    # what pip does; otherwise setuptools' `setup.py install` will try to build
    # and install an egg (deprecated behavior, but still the default for
    # backwards-compat).
    # In the future we might change this to just run pip.
    p = subprocess.Popen([sys.executable, 'setup.py', 'build',
                          f'--build-base={build_dir}', 'install',
                          f'--install-lib={install_dir}',
                          '--single-version-externally-managed',
                          f'--record={record_file}'],
                          cwd=setup_path, env=env,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Whether the process fails or not this isn't likely to produce a great
    # deal of output so communicate should be fine in almost all cases
    stdout, stderr = p.communicate()

    try:
        stdout, stderr = stdout.decode('utf8'), stderr.decode('utf8')
    except UnicodeDecodeError:
        # Don't try to guess about encoding; just display the text
        stdout, stderr = stdout.decode('latin1'), stderr.decode('latin1')

    # If compilation fails, we can skip this test, since the
    # dependencies necessary to compile an extension may be missing.
    # If it passes, however, we want to continue and ensure that the
    # extension created is actually usable.  However, if we're on
    # continuous integration setup, we
    # don't want to ever skip, because having it fail in that
    # environment probably indicates something more serious that we
    # want to know about.
    if (not ('CI' in os.environ or
             'CONTINUOUS_INTEGRATION' in os.environ) and
            p.returncode):
        pytest.skip("system unable to compile extensions")
        return

    assert p.returncode == 0, (
        "setup.py exited with non-zero return code {}\n"
        "stdout:\n\n{}\n\nstderr:\n\n{}\n".format(
            p.returncode, stdout, stderr))

    code = """
    import sys
    import wcsapi_test
    sys.exit(wcsapi_test.test())
    """

    code = code.strip().replace('\n', '; ')

    # Import and run the extension
    subprocess.check_call([sys.executable, '-c', code], env=env)
