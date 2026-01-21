import warnings


def issue_warning():
    warnings.warn("hello", UserWarning)


def test_warning_capture_safe():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        issue_warning()
        assert any(issubclass(x.category, UserWarning) for x in w)
