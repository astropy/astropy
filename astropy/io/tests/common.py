import asciitable

if asciitable.has_numpy:
    numpy_cases = (True, False)
else:
    numpy_cases = (False,)

def has_numpy_and_not_has_numpy(func):
    """Perform tests that should work for has_numpy==True and has_numpy==False"""
    def wrap():
        for numpy_case in numpy_cases:
            has_numpy = asciitable.has_numpy
            asciitable.has_numpy = numpy_case
            asciitable.core.has_numpy = numpy_case
            try:
                func(numpy=numpy_case)
            finally:
                asciitable.has_numpy = has_numpy
                asciitable.core.has_numpy = has_numpy
    wrap.__name__ = func.__name__
    return wrap

def has_numpy(func):
    """Tests that will only succeed if has_numpy == True"""
    def wrap():
        for numpy_case in numpy_cases:
            if numpy_case is True:
                func(numpy=numpy_case)
    wrap.__name__ = func.__name__
    return wrap

