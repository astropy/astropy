"""Test utility functions."""

import sys

from ..util import StringIO

class CaptureStdout(object):
    """A simple context manager for redirecting stdout to a StringIO buffer."""

    def __init__(self):
        self.io = StringIO()

    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = self.io
        return self.io

    def __exit__(self, *args, **kwargs):
        sys.stdout = self._original_stdout
        self.io.close()
