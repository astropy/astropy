import pytest
import warnings


# autouse makes this an all-coordinates-tests fixture
# this can be eliminated if/when warnings in pytest are all turned to errors (gh issue #7928)
@pytest.fixture(autouse=True)
def representation_deprecation_to_error():
    warnings.filterwarnings('error', 'The `representation` keyword/property name is deprecated in favor of `representation_type`')
    filt = warnings.filters[0]
    yield
    try:
        warnings.filters.remove(filt)
    except ValueError:
        pass  # already removed
