try:
    from inspect import signature, Parameter, Signature, BoundArguments
except ImportError:
    from ._funcsigs import signature, Parameter, Signature, BoundArguments

__all__ = ['BoundArguments', 'Parameter', 'Signature', 'signature']
