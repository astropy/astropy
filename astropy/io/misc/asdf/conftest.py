from pathlib import Path

from astropy.utils.introspection import minversion


def get_asdf_tests():
    # return a list of filenames for all ".py" files in this
    # directory and recursively in every sub directory. These
    # are the files that pytest will import while attempting
    # to find tests. This list is used below to ignore all of
    # these files if an incompatible version of ASDF is installed
    asdf_dir = Path(__file__).parent.resolve()
    paths = Path(asdf_dir).rglob("*.py")
    return [str(p.relative_to(asdf_dir)) for p in paths]


collect_ignore = get_asdf_tests()
try:
    import asdf
except ImportError:
    pass
else:
    if not minversion(asdf, "3.0.0.dev"):
        collect_ignore = []
