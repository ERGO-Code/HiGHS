# distutils: language=c++
# cython: language_level=3

cdef extern from "HConst.h" nogil:
    cdef enum HighsPrintMessageLevel:
        ML_MIN = 0
        ML_NONE = ML_MIN
        ML_VERBOSE = 1
        ML_DETAILED = 2
        ML_MINIMAL = 4
        ML_ALWAYS = ML_VERBOSE | ML_DETAILED | ML_MINIMAL
        ML_MAX = ML_ALWAYS

    cdef enum HighsBasisStatus:
        LOWER "HighsBasisStatus::LOWER" = 0, # (slack) variable is at its lower bound [including fixed variables]
        BASIC "HighsBasisStatus::BASIC" # (slack) variable is basic
        UPPER "HighsBasisStatus::UPPER" # (slack) variable is at its upper bound
        ZERO "HighsBasisStatus::ZERO" # free variable is non-basic and set to zero
        NONBASIC "HighsBasisStatus::NONBASIC" # nonbasic with no specific bound information - useful for users and postsolve
        SUPER "HighsBasisStatus::SUPER"
