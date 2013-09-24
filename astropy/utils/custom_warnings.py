class AstropyWarning(Warning):
    pass


class AstropyUserWarning(UserWarning, AstropyWarning):
    pass


class AstropyDeprecationWarning(DeprecationWarning, AstropyWarning):
    pass


class AstropyPendingDeprecationWarning(PendingDeprecationWarning, AstropyWarning):
    pass
