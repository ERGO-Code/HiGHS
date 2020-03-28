# distutils: language=c++
# cython: language_level=3

from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "HighsLp.h" nogil:
    cdef cppclass HighsLp:
        int numCol_
        int numRow_
        int numInt_
        vector[double] Avalue_
        string model_name_

    cdef enum HighsModelStatus:
        HighsModelStatusNOTSET "HighsModelStatus::NOTSET"
        HighsModelStatusLOAD_ERROR "HighsModelStatus::LOAD_ERROR"
        HighsModelStatusMODEL_ERROR "HighsModelStatus::MODEL_ERROR"
        HighsModelStatusMODEL_EMPTY "HighsModelStatus::MODEL_EMPTY"
        HighsModelStatusPRESOLVE_ERROR "HighsModelStatus::PRESOLVE_ERROR"
        HighsModelStatusSOLVE_ERROR "HighsModelStatus::SOLVE_ERROR"
        HighsModelStatusPOSTSOLVE_ERROR "HighsModelStatus::POSTSOLVE_ERROR"
        HighsModelStatusPRIMAL_INFEASIBLE "HighsModelStatus::PRIMAL_INFEASIBLE"
        HighsModelStatusPRIMAL_UNBOUNDED "HighsModelStatus::PRIMAL_UNBOUNDED"
        HighsModelStatusOPTIMAL "HighsModelStatus::OPTIMAL"
        HighsModelStatusREACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND "HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND"
        HighsModelStatusREACHED_TIME_LIMIT "HighsModelStatus::REACHED_TIME_LIMIT"
        HighsModelStatusREACHED_ITERATION_LIMIT "HighsModelStatus::REACHED_ITERATION_LIMIT"
