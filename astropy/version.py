from packaging.version import Version

try:
    from ._version import version
except ImportError:
    import warnings

    warnings.warn(
        f"could not determine {__name__.split('.')[0]} package version; "
        "this indicates a broken installation"
    )
    del warnings

    version = "0.0.0"


# We use Version to define major, minor, micro, but ignore any suffixes.
def split_version(version):
    pieces = [0, 0, 0]

    try:
        v = Version(version)
        pieces = [v.major, v.minor, v.micro]

    except Exception:
        pass

    return pieces


major, minor, bugfix = split_version(version)

del split_version  # clean up namespace.

release = "dev" not in version
