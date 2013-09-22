class AstropyWarning(Warning):
    pass


class AstropyDeprecationWarning(DeprecationWarning, AstropyWarning):
    pass


class AstropyPendingDeprecationWarning(PendingDeprecationWarning, AstropyWarning):
    pass
