class AstropyWarning(Warning):
    pass


class AstropyUserWarning(UserWarning):
    pass


class AstropyDeprecationWarning(DeprecationWarning, AstropyWarning):
    pass


class AstropyPendingDeprecationWarning(PendingDeprecationWarning, AstropyWarning):
    pass
