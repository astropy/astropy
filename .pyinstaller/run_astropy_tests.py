import os
import shutil
import sys

import erfa  # noqa
import matplotlib
import pytest

import astropy  # noqa

if len(sys.argv) == 3 and sys.argv[1] == '--astropy-root':
    ROOT = sys.argv[2]
else:
    # Make sure we don't allow any arguments to be passed - some tests call
    # sys.executable which becomes this script when producing a pyinstaller
    # bundle, but we should just error in this case since this is not the
    # regular Python interpreter.
    if len(sys.argv) > 1:
        print("Extra arguments passed, exiting early")
        sys.exit(1)

for root, dirnames, files in os.walk(os.path.join(ROOT, 'astropy')):

    # NOTE: we can't simply use
    # test_root = root.replace('astropy', 'astropy_tests')
    # as we only want to change the one which is for the module, so instead
    # we search for the last occurrence and replace that.
    pos = root.rfind('astropy')
    test_root = root[:pos] + 'astropy_tests' + root[pos + 7:]

    # Copy over the astropy 'tests' directories and their contents
    for dirname in dirnames:
        final_dir = os.path.relpath(os.path.join(test_root, dirname), ROOT)
        # We only copy over 'tests' directories, but not astropy/tests (only
        # astropy/tests/tests) since that is not just a directory with tests.
        if dirname == 'tests' and not root.endswith('astropy'):
            shutil.copytree(os.path.join(root, dirname), final_dir, dirs_exist_ok=True)
        else:
            # Create empty __init__.py files so that 'astropy_tests' still
            # behaves like a single package, otherwise pytest gets confused
            # by the different conftest.py files.
            init_filename = os.path.join(final_dir, '__init__.py')
            if not os.path.exists(os.path.join(final_dir, '__init__.py')):
                os.makedirs(final_dir, exist_ok=True)
                with open(os.path.join(final_dir, '__init__.py'), 'w') as f:
                    f.write("#")
    # Copy over all conftest.py files
    for file in files:
        if file == 'conftest.py':
            final_file = os.path.relpath(os.path.join(test_root, file), ROOT)
            shutil.copy2(os.path.join(root, file), final_file)

# Add the top-level __init__.py file
with open(os.path.join('astropy_tests', '__init__.py'), 'w') as f:
    f.write("#")

# Remove test file that tries to import all sub-packages at collection time
os.remove(os.path.join('astropy_tests', 'utils', 'iers', 'tests', 'test_leap_second.py'))

# Remove convolution tests for now as there are issues with the loading of the C extension.
# FIXME: one way to fix this would be to migrate the convolution C extension away from using
# ctypes and using the regular extension mechanism instead.
shutil.rmtree(os.path.join('astropy_tests', 'convolution'))
os.remove(os.path.join('astropy_tests', 'modeling', 'tests', 'test_convolution.py'))
os.remove(os.path.join('astropy_tests', 'modeling', 'tests', 'test_core.py'))
os.remove(os.path.join('astropy_tests', 'visualization', 'tests', 'test_lupton_rgb.py'))

# FIXME: PIL minversion check does not work
os.remove(os.path.join('astropy_tests', 'visualization', 'wcsaxes', 'tests', 'test_misc.py'))
os.remove(os.path.join('astropy_tests', 'visualization', 'wcsaxes', 'tests', 'test_wcsapi.py'))

# FIXME: The following tests rely on the fully qualified name of classes which
# don't seem to be the same.
os.remove(os.path.join('astropy_tests', 'table', 'mixins', 'tests', 'test_registry.py'))

# Copy the top-level conftest.py
shutil.copy2(os.path.join(ROOT, 'astropy', 'conftest.py'),
             os.path.join('astropy_tests', 'conftest.py'))

# matplotlib hook in pyinstaller 5.0 and later no longer collects every backend, see
# https://github.com/pyinstaller/pyinstaller/issues/6760
matplotlib.use('svg')

# We skip a few tests, which are generally ones that rely on explicitly
# checking the name of the current module (which ends up starting with
# astropy_tests rather than astropy).

SKIP_TESTS = ['test_exception_logging_origin',
              'test_log',
              'test_configitem',
              'test_config_noastropy_fallback',
              'test_no_home',
              'test_path',
              'test_rename_path',
              'test_data_name_third_party_package',
              'test_pkg_finder',
              'test_wcsapi_extension',
              'test_find_current_module_bundle',
              'test_minversion',
              'test_imports',
              'test_generate_config',
              'test_generate_config2',
              'test_create_config_file',
              'test_download_parallel_fills_cache']

# Run the tests!
sys.exit(pytest.main(['astropy_tests',
                      '-k ' + ' and '.join('not ' + test for test in SKIP_TESTS)],
                     plugins=['pytest_astropy.plugin',
                              'pytest_doctestplus.plugin',
                              'pytest_openfiles.plugin',
                              'pytest_remotedata.plugin',
                              'pytest_astropy_header.display']))
