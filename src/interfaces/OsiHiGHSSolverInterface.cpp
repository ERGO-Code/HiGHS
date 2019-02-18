// TODO license
#include "OsiHiGHSSolverInterface.hpp"

OsiHiGHSSolverInterface::OsiHiGHSSolverInterface() {
  HighsOptions options;
  this->highs = new Highs(options);
}

OsiHiGHSSolverInterface::~OsiHiGHSSolverInterface() {
   delete this->highs;
}