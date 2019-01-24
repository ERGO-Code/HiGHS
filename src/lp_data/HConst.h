/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HConst.h
 * @brief Constants for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HCONST_H_
#define LP_DATA_HCONST_H_

const int HIGHS_CONST_I_INF = 32767;
const double HIGHS_CONST_INF = 1e200;
const double HIGHS_CONST_TINY = 1e-14;
const double HIGHS_CONST_ZERO = 1e-50;

constexpr double kBoundTolerance = 1e-8;

enum ModelLogLevel {
  ML_NONE = 0,
    ML_VERBOSE = 1,
    ML_DETAILED = 2,
    ML_MINIMAL = 4,
    ML_DEFAULT = ML_VERBOSE
    };
  

enum class ParallelOption {
  OFF = 0,
  ON,
  DEFAULT = OFF
};

enum class PresolveOption {
  OFF = 0,
  ON,
  DEFAULT = OFF
};

enum CrashOption {
  OFF = 0,
  ON,
  DEFAULT = OFF
};
  
enum class SimplexOption {
  OFF = 0,
  ON,
  DEFAULT = OFF
};

const double HIGHS_RUN_TIME_LIMIT_DEFAULT = HIGHS_CONST_INF;
const double PRIMAL_FEASIBILITY_TOLERANCE_DEFAULT = 1e-7;
const double DUAL_FEASIBILITY_TOLERANCE_DEFAULT = 1e-7;

#endif /* LP_DATA_HCONST_H_ */
