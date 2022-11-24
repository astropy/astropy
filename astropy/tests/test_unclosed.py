"""
This file is an experiment to compare ResourceWarnings with pytest-openfiles
"""
import pytest


def test_closed_file(tmp_path):
    with open(tmp_path / "test1", "w") as fobj:
        pass


def test_unclosed_file(tmp_path):
    open(tmp_path / "test2", "w")


def test_evil_unclosed_file(tmp_path):
    fobj = open(tmp_path / "test3", "w")
    fobj = None


@pytest.fixture
def unclosed_fixture(tmp_path):
    return open(tmp_path / "test4", "w")


@pytest.fixture
def closed_fixture(tmp_path):
    fobj = open(tmp_path / "test5", "w")
    yield fobj
    fobj.close()


@pytest.fixture
def closed_fixture_context(tmp_path):
    with open(tmp_path / "test6", "w") as fobj:
        yield fobj


def test_closed_fixture(closed_fixture):
    pass


def test_unclosed_fixture(unclosed_fixture):
    pass


def test_closed_fixture_context(closed_fixture_context):
    pass
