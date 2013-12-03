"""
Implements the wrapper for the Astropy test runner in the form of the
``./setup.py test`` distutils command.
"""


import os
import shutil
import tempfile

from setuptools import Command

from .runner import TestRunner
from ..utils.compat.funcsigs import signature


def _fix_user_options(options):
    """
    This is for Python 2.x and 3.x compatibility.  distutils expects Command
    options to all be byte strings on Python 2 and Unicode strings on Python 3.
    """

    def to_str_or_none(x):
        if x is None:
            return None
        return str(x)

    return [tuple(to_str_or_none(x) for x in y) for y in options]


class AstropyTest(Command, object):
    description = 'Run the tests for this package'

    user_options = [
        ('package=', 'P',
         "The name of a specific package to test, e.g. 'io.fits' or 'utils'.  "
         "If nothing is specified, all default tests are run."),
        ('test-path=', 't',
         'Specify a test location by path.  If a relative path to a  .py file, '
         'it is relative to the built package, so e.g., a  leading "astropy/" '
         'is necessary.  If a relative  path to a .rst file, it is relative to '
         'the directory *below* the --docs-path directory, so a leading '
         '"docs/" is usually necessary.  May also be an absolute path.'),
        ('verbose-results', 'V',
         'Turn on verbose output from pytest.'),
        ('plugins=', 'p',
         'Plugins to enable when running pytest.'),
        ('pastebin=', 'b',
         "Enable pytest pastebin output. Either 'all' or 'failed'."),
        ('args=', 'a',
         'Additional arguments to be passed to pytest.'),
        ('remote-data', 'R', 'Run tests that download remote data.'),
        ('pep8', '8',
         'Enable PEP8 checking and disable regular tests. '
         'Requires the pytest-pep8 plugin.'),
        ('pdb', 'd',
         'Start the interactive Python debugger on errors.'),
        ('coverage', 'c',
         'Create a coverage report. Requires the coverage package.'),
        ('open-files', 'o', 'Fail if any tests leave files open.  Requires the '
         'psutil package.'),
        ('parallel=', 'j',
         'Run the tests in parallel on the specified number of '
         'CPUs.  If negative, all the cores on the machine will be '
         'used.  Requires the pytest-xdist plugin.'),
        ('docs-path=', None,
         'The path to the documentation .rst files.  If not provided, and '
         'the current directory contains a directory called "docs", that '
         'will be used.'),
        ('skip-docs', None,
         "Don't test the documentation .rst files."),
        ('repeat=', None,
         'How many times to repeat each test (can be used to check for '
         'sporadic failures).'),
        ('temp-root=', None,
         'The root directory in which to create the temporary testing files. '
         'If unspecified the system default is used (e.g. /tmp) as explained '
         'in the documentation for tempfile.mkstemp.')
    ]

    user_options = _fix_user_options(user_options)

    # Maps arguments to TestRunner.run_tests to the attributes of this class
    # that hold the value to pass for that argument
    # (Currently only one option that has a different name, due to clash with
    # a standard distutils command option.)
    runner_option_map = {'verbose': 'verbose_results'}

    package_name = ''

    def initialize_options(self):
        self.package = None
        self.test_path = None
        self.verbose_results = False
        self.plugins = None
        self.pastebin = None
        self.args = None
        self.remote_data = False
        self.pep8 = False
        self.pdb = False
        self.coverage = False
        self.open_files = False
        self.parallel = 0
        self.docs_path = None
        self.skip_docs = False
        self.repeat = None

        # temp_root is the root directory for all temp files/directories (i.e.
        # /tmp), though this may be changed with the temp-root option (e.g. if
        # the user's /tmp is too small or slow the user can specify a different
        # location for temp files created during testing).
        # temp_dir on the other hand is a temp directory created under
        # temp_root for making a test build
        # Finally, testing_path is the path we will actually change directories
        # into to run the tests--typically this is a path under temp_dir; this
        # is not user-specified and is set by _build_temp_install
        self.temp_root = None
        self.temp_dir = None
        self.testing_path = None

    def finalize_options(self):
        # Normally we would validate the options here, but that's handled in
        # run_tests
        pass

    def run(self):
        """
        Run the tests!
        """

        # Ensure there is a doc path
        if self.docs_path is None:
            if os.path.exists('docs'):
                self.docs_path = os.path.abspath('docs')

        # Build a testing install of the package
        self._build_temp_install()

        runner = TestRunner(os.path.join(self.testing_path,
                                         self.package_name))
        runner_sig = signature(runner.run_tests_in_subprocess)
        kwargs = dict((arg,
                       getattr(self, self.runner_option_map.get(arg, arg)))
                      for arg in runner_sig.parameters)

        # Run everything in a try: finally: so that the tmp dir gets deleted.
        try:
            retcode = runner.run_tests_in_subprocess(**kwargs)
        finally:
            # Remove temporary directory
            shutil.rmtree(self.temp_dir)

        raise SystemExit(retcode)

    def _build_temp_install(self):
        """
        Build the package and copy the build to a temporary directory for
        the purposes of testing this avoids creating pyc and __pycache__
        directories inside the build directory
        """

        self.reinitialize_command('build', inplace=True)
        self.run_command('build')
        build_cmd = self.get_finalized_command('build')
        new_path = os.path.abspath(build_cmd.build_lib)

        # On OSX the default path for temp files is under /var, but in most
        # cases on OSX /var is actually a symlink to /private/var; ensure we
        # dereference that link, because py.test is very sensitive to relative
        # paths...
        temp_dir = tempfile.mkdtemp(prefix=self.package_name + '-test-',
                                    dir=self.temp_root)
        self.temp_dir = os.path.realpath(temp_dir)
        self.testing_path = os.path.join(self.temp_dir,
                                         os.path.basename(new_path))
        shutil.copytree(new_path, self.testing_path)

        new_docs_path = os.path.join(self.temp_dir,
                                     os.path.basename(self.docs_path))
        shutil.copytree(self.docs_path, new_docs_path)
        self.docs_path = new_docs_path

        shutil.copy('setup.cfg', self.temp_dir)
