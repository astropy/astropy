import subprocess

import pytest


class _FakePopen:
    """
    Simulate graphviz subprocess:
    - If communicate() receives str (old buggy code), raise TypeError.
    - If receives bytes (fixed code), return binary stdout (e.g. PNG).
    """

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        self.returncode = 0

    def communicate(self, input=None):
        if isinstance(input, str):
            raise TypeError("a bytes-like object is required, not 'str'")
        return (b"PNGDATA", b"")


def test_to_dot_graph_binary_io_regression(tmp_path, monkeypatch):
    """
    This test FAILS on the old implementation because it calls:
        communicate(dotgraph_str)
    and/or writes bytes to a text file.

    It PASSES after the minimal fix:
        communicate(dotgraph.encode(...))
        open(savefn, "wb").write(stdout_bytes)
        decode stderr for error message
    """
    from astropy.coordinates.transformations.graph import TransformGraph

    monkeypatch.setattr(subprocess, "Popen", _FakePopen)

    g = TransformGraph()
    out = tmp_path / "graph.png"

    # Key: force subprocess branch (NOT the default "plain")
    g.to_dot_graph(savefn=str(out), savelayout="dot", saveformat="png")

    assert out.read_bytes() == b"PNGDATA"


def test_to_dot_graph_graphviz_error(tmp_path, monkeypatch):
    from astropy.coordinates.transformations.graph import TransformGraph

    class _FailPopen:
        def __init__(self, *args, **kwargs):
            self.returncode = 1

        def communicate(self, input=None):
            return (b"", b"graphviz failed")

    monkeypatch.setattr(subprocess, "Popen", _FailPopen)

    g = TransformGraph()
    out = tmp_path / "graph.png"

    with pytest.raises(OSError) as exc:
        g.to_dot_graph(savefn=str(out), savelayout="dot", saveformat="png")

    assert "graphviz failed" in str(exc.value)
