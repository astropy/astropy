cdef extern from "skysub.h":
    double day_of_week(double jd)


def dow(jd):
    return day_of_week(jd)
