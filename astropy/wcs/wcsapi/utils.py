import importlib


def deserialize_class(tpl, construct=True):
    """
    Deserialize classes recursively.
    """

    if not isinstance(tpl, tuple) or len(tpl) != 3:
        raise ValueError("Expected a tuple of three values")

    module, klass = tpl[0].rsplit('.', 1)
    module = importlib.import_module(module)
    klass = getattr(module, klass)

    args = tuple([deserialize_class(arg) if isinstance(arg, tuple) else arg for arg in tpl[1]])

    kwargs = dict((key, deserialize_class(val)) if isinstance(val, tuple) else (key, val) for (key, val) in tpl[2].items())

    if construct:
        return klass(*args, **kwargs)
    else:
        return klass, args, kwargs
