import sys


# SAMP relies on sockets, which are not supported on pyodide/emscripten, so we
# skip all the tests here since most of them don't and can't work. In any case,
# astropy.samp is deprecated, so this will all be removed soon.

if sys.platform == "emscripten":
    collect_ignore_glob = ["*"]
