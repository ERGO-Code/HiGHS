#ifndef __SRC_LIB_NULLSPACE_HPP__
#define __SRC_LIB_NULLSPACE_HPP__

#include "basis.hpp"
#include "matrix.hpp"
#include "runtime.hpp"

struct NullspaceReductionResult {
  HighsInt maxabsd;
  HighsInt constrainttodrop;
  Vector& d;

  bool p_in_v;

  NullspaceReductionResult(bool pinv, HighsInt mad, Vector& d_, HighsInt ctd)
      : maxabsd(mad), constrainttodrop(ctd), d(d_), p_in_v(pinv) {
    ;
  }
};

class Nullspace {
  bool uptodateZ = false;

  Basis& basis;
  Runtime& runtime;
  Matrix bufferZ;

  Vector temp_unit;
  Vector buffer_d;

  Vector buffer_unit;

  void recompute() {
    HighsInt nvar = basis.getnumactive() + basis.getnuminactive();
    Matrix Z(nvar, 0);

    const std::vector<HighsInt>& nonactive = basis.getinactive();
    const std::vector<HighsInt>& indexinfactor = basis.getindexinfactor();

    for (HighsInt i = 0; i < nonactive.size(); i++) {
      HighsInt unit = indexinfactor[nonactive[i]];

      Vector::unit(nvar, unit, temp_unit);

      basis.btran(temp_unit, temp_unit);

      Z.append(temp_unit);
    }
    bufferZ = Z;
    uptodateZ = true;
  }

 public:
  Nullspace(Runtime& rt, Basis& bas)
      : basis(bas),
        runtime(rt),
        bufferZ(Matrix(rt.instance.num_var, 0)),
        temp_unit(rt.instance.num_var),
        buffer_d(rt.instance.num_var),
        buffer_unit(rt.instance.num_var) {
    if (bas.getnuminactive() > 0) {
      recompute();
    }
    uptodateZ = true;
  }

  Vector& expand_computenewcol(HighsInt conid, Vector& target) {
    HighsInt unit = basis.getindexinfactor()[conid];
    Vector::unit(runtime.instance.num_var, unit, target);

    basis.btran(target, target);

    return target;
  }

  NullspaceReductionResult reduce(Runtime& rt, HighsInt newactivecon) {
    uptodateZ = false;

    HighsInt idx = indexof(basis.getinactive(), newactivecon);
    if (idx != -1) {
      buffer_unit.dim = basis.getinactive().size() + 1;
      Vector::unit(basis.getinactive().size(), idx, buffer_unit);
      return NullspaceReductionResult(true, idx, buffer_unit, newactivecon);
    }

    // TODO: this operation is inefficient.
    Vector aq = rt.instance.A.t().extractcol(newactivecon);
    basis.Ztprod(aq, buffer_d, true, newactivecon);

    HighsInt maxabs = 0;
    for (HighsInt i = 0; i < buffer_d.num_nz; i++) {
      if (fabs(buffer_d.value[buffer_d.index[i]]) >
          fabs(buffer_d.value[maxabs])) {
        maxabs = buffer_d.index[i];
      }
    }

    if (fabs(buffer_d.value[maxabs]) < rt.settings.d_zero_threshold) {
      printf(
          "degeneracy? not possible to find non-active constraint to "
          "leave basis. max: log(d[%" HIGHSINT_FORMAT "]) = %lf\n",
          maxabs, log10(fabs(buffer_d.value[maxabs])));
      exit(1);
    }
    return NullspaceReductionResult(false, maxabs, buffer_d,
                                    basis.getinactive()[maxabs]);
  }

  Matrix& getNullspace() {
    if (!uptodateZ) {
      recompute();
    }
    return bufferZ;
  }
};

#endif
