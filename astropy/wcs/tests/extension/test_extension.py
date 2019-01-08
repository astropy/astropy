# Licensed under a 3-clause BSD style license - see LICENSE.rst


import os
import subprocess
import sys

import pytest


def test_wcsapi_extension(tmpdir):
    # Test that we can build a simple C extension with the astropy.wcs C API

    build_dir = tmpdir.mkdir('build').strpath
    install_dir = tmpdir.mkdir('install').strpath

    setup_path = os.path.dirname(__file__)
    astropy_path = os.path.abspath(
        os.path.join(setup_path, '..', '..', '..', '..'))

    env = os.environ.copy()
    paths = [install_dir, astropy_path]
    if env.get('PYTHONPATH'):
        paths.append(env.get('PYTHONPATH'))
    env[str('PYTHONPATH')] = str(os.pathsep.join(paths))

    # Build the extension
    # This used to use subprocess.check_call, but on Python 3.4 there was
    # a mysterious Heisenbug causing this to fail with a non-zero exit code
    # *unless* the output is redirected.  This bug also did not occur in an
    # interactive session, so it likely had something to do with pytest's
    # output capture
    p = subprocess.Popen([sys.executable, 'setup.py', 'build',
                          '--build-base={0}'.format(build_dir), 'install',
                          '--install-lib={0}'.format(install_dir),
                          astropy_path], cwd=setup_path, env=env,
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
    # Travis-CI, or another generic continuous integration setup, we
    # don't want to ever skip, because having it fail in that
    # environment probably indicates something more serious that we
    # want to know about.
    if (not (str('CI') in os.environ or
             str('TRAVIS') in os.environ or
             str('CONTINUOUS_INTEGRATION') in os.environ) and
        p.returncode):
        pytest.skip("system unable to compile extensions")
        return

    assert p.returncode == 0, (
        "setup.py exited with non-zero return code {0}\n"
        "stdout:\n\n{1}\n\nstderr:\n\n{2}\n".format(
            p.returncode, stdout, stderr))

    code = """
    import sys
    import wcsapi_test
    sys.exit(wcsapi_test.test())
    """

    code = code.strip().replace('\n', '; ')

    # Import and run the extension
    subprocess.check_call([sys.executable, '-c', code], env=env)
