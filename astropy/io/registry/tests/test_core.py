# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
from io import BytesIO

import pytest

import astropy.utils.data as data_utils
from astropy.io.registry.core import UnifiedInputRegistry, _expand_user_in_args


class DummyData:
    pass


def test_expand_user_in_args_only_expands_existing_parent(monkeypatch, tmp_path):
    home = tmp_path / "home"
    existing_parent = home / "existing"
    existing_parent.mkdir(parents=True)
    monkeypatch.setenv("HOME", str(home))

    expanded = _expand_user_in_args(("~/existing/data.ecsv", "sentinel"))
    assert expanded == (str(existing_parent / "data.ecsv"), "sentinel")

    unchanged = _expand_user_in_args(("~/missing/data.ecsv", "sentinel"))
    assert unchanged == ("~/missing/data.ecsv", "sentinel")


def test_read_uses_normalized_path_and_closes_fileobj_on_reader_error(
    monkeypatch, tmp_path
):
    registry = UnifiedInputRegistry()
    seen = {}

    class DummyContext:
        def __enter__(self):
            seen["fileobj"] = BytesIO(b"payload")
            return seen["fileobj"]

        def __exit__(self, exc_type, exc, tb):
            seen["exit"] = (exc_type, exc, tb)

    def fake_get_readable_fileobj(path, encoding="binary", cache=False):
        seen["path"] = path
        seen["encoding"] = encoding
        seen["cache"] = cache
        return DummyContext()

    def identifier(origin, path, fileobj, *args, **kwargs):
        seen["identifier"] = (origin, path, fileobj)
        return True

    def reader(fileobj):
        seen["reader_arg"] = fileobj
        raise RuntimeError("boom")

    registry.register_identifier("dummy", DummyData, identifier)
    registry.register_reader("dummy", DummyData, reader)
    monkeypatch.setattr(data_utils, "get_readable_fileobj", fake_get_readable_fileobj)

    path = tmp_path / "input.dat"

    with pytest.raises(RuntimeError, match="boom"):
        registry.read(DummyData, path, cache=True)

    assert seen["path"] == os.fspath(path)
    assert isinstance(seen["path"], str)
    assert seen["encoding"] == "binary"
    assert seen["cache"] is True
    assert seen["identifier"][0] == "read"
    assert seen["identifier"][1] == os.fspath(path)
    assert seen["reader_arg"] is seen["identifier"][2]
    assert seen["exit"][0] is RuntimeError


def test_read_falls_back_to_path_when_fileobj_probe_fails(monkeypatch, tmp_path):
    registry = UnifiedInputRegistry()
    seen = {}

    def fake_get_readable_fileobj(path, encoding="binary", cache=False):
        raise RuntimeError("cannot open")

    def identifier(origin, path, fileobj, *args, **kwargs):
        seen["identifier"] = (origin, path, fileobj)
        return True

    def reader(path):
        seen["reader_arg"] = path
        return DummyData()

    registry.register_identifier("dummy", DummyData, identifier)
    registry.register_reader("dummy", DummyData, reader)
    monkeypatch.setattr(data_utils, "get_readable_fileobj", fake_get_readable_fileobj)

    path = tmp_path / "input.dat"
    result = registry.read(DummyData, path)

    assert isinstance(result, DummyData)
    assert seen["identifier"] == ("read", os.fspath(path), None)
    assert seen["reader_arg"] == os.fspath(path)
