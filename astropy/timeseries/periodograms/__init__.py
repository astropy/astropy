# Licensed under a 3-clause BSD style license - see LICENSE.rst


from lazy_loader import attach_stub

__getattr__, __dir__, __all__ = attach_stub(__name__, __file__)
