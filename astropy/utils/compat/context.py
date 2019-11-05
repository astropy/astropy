import sys
import contextlib

__all__ = ['nullcontext']

if sys.version_info[:2] < (3, 7):
    @contextlib.contextmanager
    def nullcontext():
        yield None
else:
    from contextlib import nullcontext
