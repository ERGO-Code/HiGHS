#ifndef __SRC_LIB_SETTINGS_HPP__
#define __SRC_LIB_SETTINGS_HPP__

#include "eventhandler.hpp"
#include "statistics.hpp"

enum class RatiotestStrategy { TwoPass, Textbook };

enum class PricingStrategy { SteepestEdge, DantzigWolfe, Devex };

enum class OutputLevel { LIGHT, MEDIUM, HEAVY };

enum class Phase1Strategy { HIGHS, QUASS };

struct Settings {
  RatiotestStrategy ratiotest = RatiotestStrategy::TwoPass;
  double ratiotest_t = 1E-9;
  double ratiotest_d = 1E-8;

  PricingStrategy pricing = PricingStrategy::Devex;

  double pnorm_zero_threshold = 10E-12; // if ||p|| < this threshold, p is determined to not be an improving search direction
  double improvement_zero_threshold = 10E-5; // if p^t gradient < this threshold, p is determined to not be an improving search direction
  double d_zero_threshold = 10E-13; // minimal value for pivot, will declare degeneracy if no larger pivot is found
  double lambda_zero_threshold = 10E-10; // used for pricing / optimality checking
  double pQp_zero_threshold = 10E-8; // if p'Qp < this, p is determined to not have curvature, a simplex-like iteration is performed.

  bool hessianregularization = false; // if true, a small multiple of the identity matrix will be added to the Hessian
  double hessianregularizationfactor = 1E-7; // multiple of identity matrix added to hessian in case of regularization

  Phase1Strategy phase1strategy = Phase1Strategy::HIGHS;
  bool phase1movefreevarsbasic = false;
  bool phase1boundfreevars = false;

  OutputLevel outputlevel = OutputLevel::LIGHT;
  HighsInt reportingfequency = 1;
  Eventhandler<Statistics&> endofiterationevent;

  HighsInt reinvertfrequency = 100;
  HighsInt gradientrecomputefrequency = 100;
  HighsInt reducedgradientrecomputefrequency =
      std::numeric_limits<HighsInt>::infinity();
  HighsInt reducedhessianrecomputefrequency =
      std::numeric_limits<HighsInt>::infinity();

  HighsInt iterationlimit = std::numeric_limits<HighsInt>::infinity();
  double timelimit = std::numeric_limits<double>::infinity();

  bool rowscaling = true;
  bool varscaling = true;

  bool perturbation = false;
};

#endif
