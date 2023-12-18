__all__ = ["get_include"]

def get_include():
    """
    Get the path to astropy.wcs's C header files.
    """
    import os

    return os.path.join(os.path.dirname(__file__), "include")
