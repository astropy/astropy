from contextlib import contextmanager

_inf_dist_enabled = False  # Global flag

def get_inf_dist():
    return _inf_dist_enabled  # Return the current setting

@contextmanager
def set_inf_dist(enable=True):
    global _inf_dist_enabled
    old_value = _inf_dist_enabled
    _inf_dist_enabled = enable
    try:
        yield  # This allows use in a `with` statement
    finally:
        _inf_dist_enabled = old_value  # Reset after use
