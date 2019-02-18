/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file interfaces/OsiHiGHSInterface.cpp
 * @brief Osi/HiGHS interface implementation
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "OsiHiGHSSolverInterface.hpp"

OsiHiGHSSolverInterface::OsiHiGHSSolverInterface() {
  HighsOptions options;
  this->highs = new Highs(options);
}

OsiHiGHSSolverInterface::~OsiHiGHSSolverInterface() {
   delete this->highs;
}