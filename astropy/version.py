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

    # this branch may be executed in rare but still valid conditions. For instance,
    # if astropy is installed from a shallow clone, i.e., without most of the vcs
    # history, which is needed to compute the dev version number dynamically.

    # default to a virtually infinite version number, never to be reached, that should
    # - alert the reader if seen in production
    # - satisfy any real-world lower bound on astropy's version
    version = "999.999.999"


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
