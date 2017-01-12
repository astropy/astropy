from astropy.tests.runner import TestRunner as T_Runner
from astropy.tests.runner import TestRunnerBase as T_RunnerBase
from astropy.tests.runner import keyword
from astropy.tests.helper import pytest


def test_disable_kwarg():
    class no_remote_data(T_Runner):
        @keyword()
        def remote_data(self, remote_data, kwargs):
            return NotImplemented

    r = no_remote_data('.')
    with pytest.raises(TypeError):
        r.run_tests(remote_data='bob')


def test_wrong_kwarg():
    r = T_Runner('.')
    with pytest.raises(TypeError):
        r.run_tests(spam='eggs')


def test_invalid_kwarg():
    class bad_return(T_RunnerBase):
        @keyword()
        def remote_data(self, remote_data, kwargs):
            return 'bob'

    r = bad_return('.')
    with pytest.raises(TypeError):
        r.run_tests(remote_data='bob')


def test_new_kwarg():
    class Spam(T_RunnerBase):
        @keyword()
        def spam(self, spam, kwargs):
            return [spam]

    r = Spam('.')

    args = r._generate_args(spam='spam')

    assert ['spam'] == args


def test_priority():
    class Spam(T_RunnerBase):
        @keyword()
        def spam(self, spam, kwargs):
            return [spam]

        @keyword(priority=1)
        def eggs(self, eggs, kwargs):
            return [eggs]

    r = Spam('.')

    args = r._generate_args(spam='spam', eggs='eggs')

    assert ['eggs', 'spam'] == args


def test_docs():
    class Spam(T_RunnerBase):
        @keyword()
        def spam(self, spam, kwargs):
            """
            Spam Spam Spam
            """
            return [spam]

        @keyword()
        def eggs(self, eggs, kwargs):
            """
            eggs asldjasljd
            """
            return [eggs]

    r = Spam('.')
    assert "eggs" in r.run_tests.__doc__
    assert "Spam Spam Spam" in r.run_tests.__doc__
