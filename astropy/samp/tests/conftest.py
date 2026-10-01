import sys

if sys.platform == 'emscripten':
    collect_ignore_glob = ['*']
