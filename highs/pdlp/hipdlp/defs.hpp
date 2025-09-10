#ifndef HIGHS_PDLP_DEFS_HPP
#define HIGHS_PDLP_DEFS_HPP

#include <cmath>
#include <vector>
#include "Highs.h"

enum class Device { CPU, GPU };

enum class ScalingMethod { NONE, RUIZ, POCK_CHAMBOLLE, L2_NORM, COMBINED };

enum class RestartStrategy { NO_RESTART, FIXED_RESTART, ADAPTIVE_RESTART };

enum class StepSizeStrategy { FIXED, ADAPTIVE, MALITSKY_POCK };

struct MalitskyPockParams {
  double step_size_interpolation = 0.5;  // Between 0 and 1
  double step_size_downscaling_factor = 0.7;
  double linesearch_contraction_factor = 0.99;
  void initialise();
};

struct AdaptiveLinesearchParams {
  double step_size_reduction_exponent = 0.3;
  double step_size_growth_exponent = 0.6;
  void initialise();
};

struct PrimalDualParams {
  double eta;
  double omega;
  double tolerance;
  size_t max_iterations;
  Device device_type;
  double time_limit = 3600.0;

  // Restart parameters
  RestartStrategy restart_strategy;
  int fixed_restart_interval;

  bool use_halpern_restart = false;

  // Scaling parameters
  ScalingMethod scaling_method = ScalingMethod::NONE;
  bool use_ruiz_scaling = false;
  bool use_pc_scaling = false;
  bool use_l2_scaling = false;

  // Ruiz scaling parameters
  int ruiz_iterations = 10;
  double ruiz_norm = INFINITY;

  // Pock-Chambolle scaling parameters
  double pc_alpha = 1.0;

  // Step sizes strategy
  StepSizeStrategy step_size_strategy = StepSizeStrategy::FIXED;

  MalitskyPockParams malitsky_pock_params;
  AdaptiveLinesearchParams adaptive_linesearch_params;
  void initialise();
};

struct PdlpIterate {
  // Primary variables 
  std::vector<double> x;
  std::vector<double> y;

  // Cached matrix-vector products
  mutable std::vector<double> Ax;
  mutable std::vector<double> Aty;
  mutable bool Ax_valid = false;
  mutable bool Aty_valid = false;

  // Constructors
  PdlpIterate() = default;
  PdlpIterate(int num_cols, int num_rows)
      : x(num_cols, 0.0), y(num_rows, 0.0), Ax(num_rows, 0.0), Aty(num_cols, 0.0) {}
  PdlpIterate(const std::vector<double>& x_init, const std::vector<double>& y_init)
      : x(x_init), y(y_init), Ax(y_init.size(), 0.0), Aty(x_init.size(), 0.0) {}

  // Arithmetic operations 
  // z = alpha * this + beta * other
  void LinearCombination(const PdlpIterate& other, double alpha, double beta);

  // this = this + alpha * other
  void AddScaled(const PdlpIterate& other, double alpha);

  // this = alpha * this
  void Scale(double alpha);

  // Copy operations
  void CopyFrom(const PdlpIterate& other);
  PdlpIterate& operator=(const PdlpIterate& other);

  // Norms and metrics
  double PrimalNorm() const;  // ||x||_2
  double DualNorm() const;  // ||y||_2
  double WeightedNorm(double omega) const;  // sqrt(||x||_2^2 + omega^2 * ||y||_2^2)

  // Distance metrics
  double Distance(const PdlpIterate& other, double omega = 1.0) const;

  // Matrix-vector product management
  void ComputeAx(const HighsLp& lp) const;
  void ComputeATy(const HighsLp& lp) const;
  void InvalidateProducts();  // Call when x or y change
  
  const std::vector<double>& GetAx(const HighsLp& lp) const;
  const std::vector<double>& GetATy(const HighsLp& lp) const;

  // For block-structured problems
  struct BlockStructure {
      std::vector<int> x_block_sizes;
      std::vector<int> y_block_sizes;
      // std::vector<std::vector<double>> x_blocks;  // Future: for block problems
      // std::vector<std::vector<double>> y_blocks;
  };

  // Optional: block structure for future extensions
  std::unique_ptr<BlockStructure> block_structure = nullptr;

private:
  void EnsureAxComputed(const HighsLp& lp) const;
  void EnsureATyComputed(const HighsLp& lp) const;
};

namespace pdlp_iterate_ops {
  // Compute z_new = z_old - 
};

#endif