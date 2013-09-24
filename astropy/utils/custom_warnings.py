class AstropyWarning(Warning):
    pass


class AstropyUserWarning(UserWarning):
    pass


class AstropyDeprecationWarning(DeprecationWarning, AstropyUserWarning):
    pass


class AstropyPendingDeprecationWarning(PendingDeprecationWarning, AstropyUserWarning):
    pass
