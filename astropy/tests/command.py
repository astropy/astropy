"""Implements the wrapper for the Astropy test runner.

This is for backward-compatibility for other downstream packages and can be removed
once astropy-helpers has reached end-of-life.

"""

import os
import shutil
import stat
import subprocess
import sys
import tempfile
from contextlib import contextmanager

from setuptools import Command

from astropy.logger import log
from astropy.utils.decorators import deprecated


@contextmanager
def _suppress_stdout():
    """
    A context manager to temporarily disable stdout.

    Used later when installing a temporary copy of astropy to avoid a
    very verbose output.
    """
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout


@deprecated("6.0")
class FixRemoteDataOption(type):
    """
    This metaclass is used to catch cases where the user is running the tests
    with --remote-data. We've now changed the --remote-data option so that it
    takes arguments, but we still want --remote-data to work as before and to
    enable all remote tests. With this metaclass, we can modify sys.argv
    before setuptools try to parse the command-line options.
    """

    def __init__(cls, name, bases, dct):
        try:
            idx = sys.argv.index("--remote-data")
        except ValueError:
            pass
        else:
            sys.argv[idx] = "--remote-data=any"

        try:
            idx = sys.argv.index("-R")
        except ValueError:
            pass
        else:
            sys.argv[idx] = "-R=any"

        return super().__init__(name, bases, dct)


@deprecated("6.0")
class AstropyTest(Command, metaclass=FixRemoteDataOption):
    description = "Run the tests for this package"

    user_options = [
        (
            "package=",
            "P",
            "The name of a specific package to test, e.g. 'io.fits' or 'utils'. "
            "Accepts comma separated string to specify multiple packages. "
            "If nothing is specified, all default tests are run.",
        ),
        (
            "test-path=",
            "t",
            "Specify a test location by path.  If a relative path to a  .py file, "
            'it is relative to the built package, so e.g., a  leading "astropy/" '
            "is necessary.  If a relative  path to a .rst file, it is relative to "
            "the directory *below* the --docs-path directory, so a leading "
            '"docs/" is usually necessary.  May also be an absolute path.',
        ),
        ("verbose-results", "V", "Turn on verbose output from pytest."),
        ("plugins=", "p", "Plugins to enable when running pytest."),
        ("pastebin=", "b", "Enable pytest pastebin output. Either 'all' or 'failed'."),
        ("args=", "a", "Additional arguments to be passed to pytest."),
        (
            "remote-data=",
            "R",
            "Run tests that download remote data. Should be "
            "one of none/astropy/any (defaults to none).",
        ),
        (
            "pep8",
            "8",
            "Enable PEP8 checking and disable regular tests. "
            "Requires the pytest-pep8 plugin.",
        ),
        ("pdb", "d", "Start the interactive Python debugger on errors."),
        ("coverage", "c", "Create a coverage report. Requires the coverage package."),
        (
            "parallel=",
            "j",
            "Run the tests in parallel on the specified number of "
            'CPUs.  If "auto", all the cores on the machine will be '
            "used.  Requires the pytest-xdist plugin.",
        ),
        (
            "docs-path=",
            None,
            "The path to the documentation .rst files.  If not provided, and "
            'the current directory contains a directory called "docs", that '
            "will be used.",
        ),
        ("skip-docs", None, "Don't test the documentation .rst files."),
        (
            "repeat=",
            None,
            "How many times to repeat each test (can be used to check for "
            "sporadic failures).",
        ),
        (
            "temp-root=",
            None,
            "The root directory in which to create the temporary testing files. "
            "If unspecified the system default is used (e.g. /tmp) as explained "
            "in the documentation for tempfile.mkstemp.",
        ),
        (
            "verbose-install",
            None,
            "Turn on terminal output from the installation of astropy in a "
            "temporary folder.",
        ),
        ("readonly", None, "Make the temporary installation being tested read-only."),
    ]

    package_name = ""

    def initialize_options(self):
        self.package = None
        self.test_path = None
        self.verbose_results = False
        self.plugins = None
        self.pastebin = None
        self.args = None
        self.remote_data = "none"
        self.pep8 = False
        self.pdb = False
        self.coverage = False
        self.parallel = 0
        self.docs_path = None
        self.skip_docs = False
        self.repeat = None
        self.temp_root = None
        self.verbose_install = False
        self.readonly = False

    def finalize_options(self):
        # Normally we would validate the options here, but that's handled in
        # run_tests
        pass

    def generate_testing_command(self):
        """
        Build a Python script to run the tests.
        """
        cmd_pre = ""  # Commands to run before the test function
        cmd_post = ""  # Commands to run after the test function

        if self.coverage:
            pre, post = self._generate_coverage_commands()
            cmd_pre += pre
            cmd_post += post

        set_flag = "import builtins; builtins._ASTROPY_TEST_ = True"

        cmd = (  # see _build_temp_install below
            "{cmd_pre}{0}; import {1.package_name}, sys; result = ("
            "{1.package_name}.test("
            "package={1.package!r}, "
            "test_path={1.test_path!r}, "
            "args={1.args!r}, "
            "plugins={1.plugins!r}, "
            "verbose={1.verbose_results!r}, "
            "pastebin={1.pastebin!r}, "
            "remote_data={1.remote_data!r}, "
            "pep8={1.pep8!r}, "
            "pdb={1.pdb!r}, "
            "parallel={1.parallel!r}, "
            "docs_path={1.docs_path!r}, "
            "skip_docs={1.skip_docs!r}, "
            "add_local_eggs_to_path=True, "
            "repeat={1.repeat!r})); "
            "{cmd_post}"
            "sys.exit(result)"
        )
        return cmd.format(set_flag, self, cmd_pre=cmd_pre, cmd_post=cmd_post)

    def run(self):
        """Run the tests!"""
        # Install the runtime dependencies.
        if self.distribution.install_requires:
            self.distribution.fetch_build_eggs(self.distribution.install_requires)

        # Ensure there is a doc path
        if self.docs_path is None:
            cfg_docs_dir = self.distribution.get_option_dict("build_docs").get(
                "source_dir", None
            )

            # Some affiliated packages use this.
            # See astropy/package-template#157
            if cfg_docs_dir is not None and os.path.exists(cfg_docs_dir[1]):
                self.docs_path = os.path.abspath(cfg_docs_dir[1])

            # fall back on a default path of "docs"
            elif os.path.exists("docs"):  # pragma: no cover
                self.docs_path = os.path.abspath("docs")

        # Build a testing install of the package
        self._build_temp_install()

        # Install the test dependencies
        # NOTE: we do this here after _build_temp_install because there is
        # a weird but which occurs if psutil is installed in this way before
        # astropy is built, Cython can have segmentation fault. Strange, eh?
        if self.distribution.tests_require:
            self.distribution.fetch_build_eggs(self.distribution.tests_require)

        # Copy any additional dependencies that may have been installed via
        # tests_requires or install_requires. We then pass the
        # add_local_eggs_to_path=True option to package.test() to make sure the
        # eggs get included in the path.
        if os.path.exists(".eggs"):
            shutil.copytree(".eggs", os.path.join(self.testing_path, ".eggs"))

        # This option exists so that we can make sure that the tests don't
        # write to an installed location.
        if self.readonly:
            log.info("changing permissions of temporary installation to read-only")
            self._change_permissions_testing_path(writable=False)

        # Run everything in a try: finally: so that the tmp dir gets deleted.
        try:
            # Construct this modules testing command
            cmd = self.generate_testing_command()

            # Run the tests in a subprocess--this is necessary since
            # new extension modules may have appeared, and this is the
            # easiest way to set up a new environment

            testproc = subprocess.Popen(
                [sys.executable, "-c", cmd], cwd=self.testing_path, close_fds=False
            )
            retcode = testproc.wait()
        except KeyboardInterrupt:
            import signal

            # If a keyboard interrupt is handled, pass it to the test
            # subprocess to prompt pytest to initiate its teardown
            testproc.send_signal(signal.SIGINT)
            retcode = testproc.wait()
        finally:
            # Remove temporary directory
            if self.readonly:
                self._change_permissions_testing_path(writable=True)
            shutil.rmtree(self.tmp_dir)

        raise SystemExit(retcode)

    def _build_temp_install(self):
        """
        Install the package and to a temporary directory for the purposes of
        testing. This allows us to test the install command, include the
        entry points, and also avoids creating pyc and __pycache__ directories
        inside the build directory.
        """
        # On OSX the default path for temp files is under /var, but in most
        # cases on OSX /var is actually a symlink to /private/var; ensure we
        # dereference that link, because pytest is very sensitive to relative
        # paths...

        tmp_dir = tempfile.mkdtemp(
            prefix=self.package_name + "-test-", dir=self.temp_root
        )
        self.tmp_dir = os.path.realpath(tmp_dir)

        log.info(f"installing to temporary directory: {self.tmp_dir}")

        # We now install the package to the temporary directory. We do this
        # rather than build and copy because this will ensure that e.g. entry
        # points work.
        self.reinitialize_command("install")
        install_cmd = self.distribution.get_command_obj("install")
        install_cmd.prefix = self.tmp_dir
        if self.verbose_install:
            self.run_command("install")
        else:
            with _suppress_stdout():
                self.run_command("install")

        # We now get the path to the site-packages directory that was created
        # inside self.tmp_dir
        install_cmd = self.get_finalized_command("install")
        self.testing_path = install_cmd.install_lib

        # Ideally, docs_path is set properly in run(), but if it is still
        # not set here, do not pretend it is, otherwise bad things happen.
        # See astropy/package-template#157
        if self.docs_path is not None:
            new_docs_path = os.path.join(
                self.testing_path, os.path.basename(self.docs_path)
            )
            shutil.copytree(self.docs_path, new_docs_path)
            self.docs_path = new_docs_path

        shutil.copy("pyproject.toml", self.testing_path)

    def _change_permissions_testing_path(self, writable=False):
        if writable:
            basic_flags = stat.S_IRUSR | stat.S_IWUSR
        else:
            basic_flags = stat.S_IRUSR
        for root, dirs, files in os.walk(self.testing_path):
            for dirname in dirs:
                os.chmod(os.path.join(root, dirname), basic_flags | stat.S_IXUSR)
            for filename in files:
                os.chmod(os.path.join(root, filename), basic_flags)

    def _generate_coverage_commands(self):
        """
        This method creates the post and pre commands if coverage is to be
        generated.
        """
        if self.parallel != 0:
            raise ValueError("--coverage can not be used with --parallel")

        try:
            import coverage  # noqa: F401
        except ImportError:
            raise ImportError(
                "--coverage requires that the coverage package is installed."
            )

        # Don't use get_pkg_data_filename here, because it
        # requires importing astropy.config and thus screwing
        # up coverage results for those packages.
        coveragerc = os.path.join(
            self.testing_path,
            self.package_name.replace(".", "/"),
            "tests",
            "coveragerc",
        )

        with open(coveragerc) as fd:
            coveragerc_content = fd.read()

        coveragerc_content = coveragerc_content.replace(
            "{packagename}", self.package_name.replace(".", "/")
        )
        tmp_coveragerc = os.path.join(self.tmp_dir, "coveragerc")
        with open(tmp_coveragerc, "wb") as tmp:
            tmp.write(coveragerc_content.encode("utf-8"))

        cmd_pre = (
            "import coverage; cov ="
            f' coverage.coverage(data_file=r"{os.path.abspath(".coverage")}",'
            f' config_file=r"{os.path.abspath(tmp_coveragerc)}"); cov.start();'
        )
        cmd_post = (
            "cov.stop(); from astropy.tests.helper import _save_coverage;"
            f' _save_coverage(cov, result, r"{os.path.abspath(".")}",'
            f' r"{os.path.abspath(self.testing_path)}");'
        )

        return cmd_pre, cmd_post
