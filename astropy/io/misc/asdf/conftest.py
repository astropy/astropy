from pathlib import Path
from astropy.utils.introspection import minversion


def get_asdf_tests():
    asdf_dir = Path(__file__).parent.resolve()
    paths = Path(asdf_dir).rglob("test_*.py")

    return [str(p.relative_to(asdf_dir)) for p in paths]


collect_ignore = get_asdf_tests()
try:
    import asdf
except ImportError:
    pass
else:
    if not minversion(asdf, "3.0.0"):
        collect_ignore = []
