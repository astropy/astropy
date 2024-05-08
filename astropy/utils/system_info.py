r"""
A helper script to automate system info in bug reports.

This can be run as::

    python -c "import astropy; astropy.system_info()"

(note that the interpreter might be called python3 instead of python, depending
on your system's details)

This script should not require anything else than the standard library and is
subject to change without notice.
"""

__all__ = ["system_info"]

import platform
from importlib.metadata import version
from importlib.util import find_spec


def _header(name: str) -> list[str]:
    return [name, "-" * len(name)]


def _report_platform() -> list[str]:
    return [
        f"{platform.platform() = }",
        f"{platform.version() = }",
        f"{platform.python_version() = }",
    ]


def _report_packages() -> list[str]:
    pkgs = [
        "astropy",
        "numpy",
        "scipy",
        "matplotlib",
        "pandas",
    ]
    lines = [f"{p:<20} {'--' if find_spec(p) is None else version(p)}" for p in pkgs]

    # pyerfa may not export its package metadata (as of 2.0.1.4)
    try:
        import erfa
    except ImportError:
        v = "--"
    else:
        v = erfa.__version__
    p = "pyerfa"
    lines.append(f"{p:<20} {v}")
    return lines


def system_info() -> None:
    """Print relevant system information for astropy bug reports.

    Examples
    --------
    >>> import astropy
    >>> astropy.system_info()  # doctest: +ELLIPSIS
    platform
    --------
    platform.platform() = ...
    platform.version() = ...
    platform.python_version() = ...
    <BLANKLINE>
    packages
    --------
    astropy              ...
    numpy                ...
    scipy                ...
    matplotlib           ...
    pandas               ...
    pyerfa               ...

    """
    report_lines = (
        *_header("platform"),
        *_report_platform(),
        "",
        *_header("packages"),
        *_report_packages(),
    )
    print("\n".join(report_lines))
