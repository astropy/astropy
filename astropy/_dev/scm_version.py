# Try to use setuptools_scm to get the current version; this is only used
# in development installations from the git repository.
import os.path as pth

try:
    from setuptools_scm import get_version

    version = get_version(root=pth.join("..", ".."), relative_to=__file__)
except Exception:
    raise ImportError("setuptools_scm broken or not installed")
