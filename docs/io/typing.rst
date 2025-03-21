********************************
I/O Typing (`astropy.io.typing`)
********************************


``astropy.io`` provides type annotations through the :mod:`astropy.io.typing` module.
These type annotations allow users to specify the expected types of variables, function
parameters, and return values when working with I/O. By using type annotations,
developers can improve code readability, catch potential type-related errors early, and
enable better code documentation and tooling support.

For example, the following function uses type annotations to specify that the
``filename`` parameter can be any type of path-like object (e.g. a string, byte-string,
or pathlib.Path object).

.. code-block:: python

    from astropy.io import fits
    from astropy.io.typing import PathLike

    def read_fits_file(filename: PathLike) -> fits.HDUList:
         return fits.open(filename)


The :mod:`astropy.io.typing` module also provides type aliases for file-like objects
that support reading and writing. The following example uses the
:class:`~astropy.io.typing.ReadableFileLike` type alias to specify that the ``fileobj``
parameter can be any file-like object that supports reading. Using a
:class:`~typing.TypeVar`, the return type of the function is specified to be the same
type as the file-like object can read.


.. code-block:: python

    from typing import TypeVar
    from astropy.io.typing import ReadableFileLike

    R = TypeVar('R')  # type of object returned by fileobj.read()

    def read_from_file(fileobj: ReadableFileLike[R]) -> R:
         """Reads from a file-like object and returns the result."""
         return fileobj.read()


Reference/API
=============

.. automodapi:: astropy.io.typing
