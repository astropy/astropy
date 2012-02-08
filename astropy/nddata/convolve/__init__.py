from astropy import setup_helpers

if setup_helpers.is_in_build_mode():
    pass
else:
    from .convolve import convolve
