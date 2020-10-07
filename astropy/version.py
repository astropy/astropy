# NOTE: First try _dev.scm_version if it exists and setuptools_scm is installed
# This file is not included in astropy wheels/tarballs, so otherwise it will
# fall back on the generated _version module.
try:
    try:
        from ._dev.scm_version import version
    except ImportError:
        from ._version import version
except Exception:
    import warnings
    warnings.warn(
        f'could not determine {__name__.split(".")[0]} package version; '
        f'this indicates a broken installation')
    del warnings

    version = '0.0.0'


# We use LooseVersion to define major, minor, micro, but ignore any suffixes.
def split_version(version):
    pieces = [0, 0, 0]

    try:
        from distutils.version import LooseVersion

        for j, piece in enumerate(LooseVersion(version).version[:3]):
            pieces[j] = int(piece)

    except Exception:
        pass

    return pieces


major, minor, bugfix = split_version(version)

del split_version  # clean up namespace.

release = 'dev' not in version
