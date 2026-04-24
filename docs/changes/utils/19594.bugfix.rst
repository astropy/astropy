The dummy file object used by ``astropy.utils.misc.silence`` to replace
``sys.stdout``/``sys.stderr`` now implements ``flush()`` and ``isatty()``,
so importing libraries that probe the stream (e.g. IPython 9.13 at import
time) no longer raises ``AttributeError`` under ``silence``.
