import pytest

ignore_pformat_default_pending_depr_warning = pytest.mark.filterwarnings(
    r"ignore:The default value for the (max_lines|max_width) argument "
    r"in \w+\.pformat:PendingDeprecationWarning"
)
