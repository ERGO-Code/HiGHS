#ifndef __SRC_LIB_SETTINGS_HPP__
#define __SRC_LIB_SETTINGS_HPP__

enum class RatiotestStrategy { TwoPass, Textbook };

enum class PricingStrategy { DantzigWolfe, Devex };

enum class OutputLevel { LIGHT, MEDIUM, HEAVY };

struct Settings {
  RatiotestStrategy ratiotest = RatiotestStrategy::TwoPass;
  double ratiotest_t = 1E-9;
  double ratiotest_d = 1E-8;

  PricingStrategy pricing = PricingStrategy::Devex;

  double pnorm_zero_threshold = 10E-12;
  double d_zero_threshold = 10E-13;
  double lambda_zero_threshold = 10E-10;

  double semidefiniteregularization = 0.0;

  bool phase1movefreevarsbasic = false;

  OutputLevel outputlevel = OutputLevel::LIGHT;
  HighsInt reportingfequency = 100;

  HighsInt reinvertfrequency = 100;
  HighsInt gradientrecomputefrequency = 1;
  HighsInt reducedgradientrecomputefrequency =
      std::numeric_limits<HighsInt>::infinity();
  HighsInt reducedhessianrecomputefrequency =
      std::numeric_limits<HighsInt>::infinity();

  HighsInt iterationlimit = std::numeric_limits<HighsInt>::infinity();
  double timelimit = std::numeric_limits<double>::infinity();

  bool rowscaling = true;
  bool varscaling = true;
};

#endif
