#ifndef __SRC_LIB_PRICING_DEVEXHARRISPRICING_HPP__
#define __SRC_LIB_PRICING_DEVEXHARRISPRICING_HPP__

#include "qpsolver/basis.hpp"
#include "qpsolver/pricing.hpp"
#include "qpsolver/runtime.hpp"

// 44839, 78965849.088174, 559, 138.681866, 0.000671, 7998

class DevexHarrisPricing : public Pricing {
 private:
  Runtime& runtime;
  Basis& basis;
  ReducedCosts& redcosts;
  std::vector<double> weights;

  HighsInt chooseconstrainttodrop(const QpVector& lambda) {
    auto activeconstraintidx = basis.getactive();
    auto constraintindexinbasisfactor = basis.getindexinfactor();

    HighsInt minidx = -1;
    double maxabslambda = 0.0;
    for (size_t i = 0; i < activeconstraintidx.size(); i++) {
      HighsInt indexinbasis =
          constraintindexinbasisfactor[activeconstraintidx[i]];
      if (indexinbasis == -1) {
        printf("error\n");
      }
      assert(indexinbasis != -1);

      double val = lambda.value[indexinbasis] * lambda.value[indexinbasis] /
                   weights[indexinbasis];
      if (val > maxabslambda && fabs(lambda.value[indexinbasis]) >
                                    runtime.settings.lambda_zero_threshold) {
        if (basis.getstatus(activeconstraintidx[i]) ==
                BasisStatus::ActiveAtLower &&
            -lambda.value[indexinbasis] > 0) {
          minidx = activeconstraintidx[i];
          maxabslambda = val;
        } else if (basis.getstatus(activeconstraintidx[i]) ==
                       BasisStatus::ActiveAtUpper &&
                   lambda.value[indexinbasis] > 0) {
          minidx = activeconstraintidx[i];
          maxabslambda = val;
        } else {
          // TODO
        }
      }
    }

    return minidx;
  }

 public:
  DevexHarrisPricing(Runtime& rt, Basis& bas, ReducedCosts& rc)
      : runtime(rt),
        basis(bas), 
        redcosts(rc),
        weights(std::vector<double>(rt.instance.num_var, 1.0)) {};

  HighsInt price(const QpVector& x, const QpVector& gradient) {
    HighsInt minidx = chooseconstrainttodrop(redcosts.getReducedCosts());
    return minidx;
  }

  void recompute() {
    // do nothing
  }

  void update_weights(const QpVector& aq, const QpVector& ep, HighsInt p,
                      HighsInt q) {
    HighsInt rowindex_p = basis.getindexinfactor()[p];
    double weight_p = weights[rowindex_p];
    for (HighsInt i = 0; i < runtime.instance.num_var; i++) {
      if (i == rowindex_p) {
        weights[i] =
            1 / (aq.value[rowindex_p] * aq.value[rowindex_p]) * weight_p;
      } else {
        weights[i] =
            max(weights[i], (aq.value[i] * aq.value[i]) /
                                (aq.value[rowindex_p] * aq.value[rowindex_p]) *
                                weight_p * weight_p);
      }
      if (weights[i] > 10E6) {
        weights[i] = 1.0;
      }
    }
  }
};

#endif
