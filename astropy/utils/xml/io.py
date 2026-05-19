# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""I/O helpers for XML-related file writing."""

# STDLIB
import codecs
import contextlib
import gzip
import io
import os

__all__ = ["convert_to_writable_filelike"]


@contextlib.contextmanager
def convert_to_writable_filelike(fd, compressed=False):
    """Return a writable file-like object for streaming output."""
    if isinstance(fd, str):
        fd = os.path.expanduser(fd)
        if fd.endswith(".gz") or compressed:
            with gzip.GzipFile(filename=fd, mode="wb") as real_fd:
                encoded_fd = io.TextIOWrapper(real_fd, encoding="utf8")
                yield encoded_fd
                encoded_fd.flush()
                real_fd.flush()
                return
        else:
            with open(fd, "w", encoding="utf8") as real_fd:
                yield real_fd
                return
    elif hasattr(fd, "write"):
        assert callable(fd.write)

        if compressed:
            fd = gzip.GzipFile(fileobj=fd, mode="wb")

        needs_wrapper = False
        try:
            fd.write("")
        except TypeError:
            needs_wrapper = True

        if not hasattr(fd, "encoding") or fd.encoding is None:
            needs_wrapper = True

        if needs_wrapper:
            yield codecs.getwriter("utf-8")(fd)
        else:
            yield fd

        fd.flush()
        if isinstance(fd, gzip.GzipFile):
            fd.close()

        return
    else:
        raise TypeError("Can not be coerced to writable file-like object")
