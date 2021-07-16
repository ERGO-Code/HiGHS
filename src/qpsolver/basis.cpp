#include "basis.hpp"

#include <cassert>
#include <memory>

Basis::Basis(Runtime& rt, std::vector<HighsInt> active,
             std::vector<BasisStatus> lower, std::vector<HighsInt> inactive)
    : runtime(rt),
      buffer_column_aq(rt.instance.num_var),
      buffer_row_ep(rt.instance.num_var) {
  for (HighsInt i = 0; i < active.size(); i++) {
    activeconstraintidx.push_back(active[i]);
    basisstatus[activeconstraintidx[i]] = lower[i];
  }
  for (HighsInt i : inactive) {
    nonactiveconstraintsidx.push_back(i);
  }

  Atran = rt.instance.A.t();

  build();
}

void Basis::build() {
  updatessinceinvert = 0;

  baseindex =
      new HighsInt[activeconstraintidx.size() + nonactiveconstraintsidx.size()];
  constraintindexinbasisfactor.clear();

  basisfactor = HFactor();

  constraintindexinbasisfactor.assign(Atran.num_row + Atran.num_col, -1);
  assert(nonactiveconstraintsidx.size() + activeconstraintidx.size() ==
         Atran.num_row);

  HighsInt counter = 0;
  for (HighsInt i : nonactiveconstraintsidx) {
    baseindex[counter++] = i;
  }
  for (HighsInt i : activeconstraintidx) {
    baseindex[counter++] = i;
  }

  const bool empty_matrix = (int)Atran.index.size() == 0;
  if (empty_matrix) {
    // The index/value vectors have size zero if the matrix has no
    // columns. However, in the Windows build, referring to index 0 of a
    // vector of size zero causes a failure, so resize to 1 to prevent
    // this.
    assert(Atran.num_col == 0);
    Atran.index.resize(1);
    Atran.value.resize(1);
  }
  basisfactor.setup(Atran.num_col, Atran.num_row, (HighsInt*)&Atran.start[0],
                    (HighsInt*)&Atran.index[0], (const double*)&Atran.value[0],
                    baseindex);
  basisfactor.build();

  for (size_t i = 0;
       i < activeconstraintidx.size() + nonactiveconstraintsidx.size(); i++) {
    constraintindexinbasisfactor[baseindex[i]] = i;
  }
}

void Basis::rebuild() {
  updatessinceinvert = 0;
  constraintindexinbasisfactor.clear();

  constraintindexinbasisfactor.assign(Atran.num_row + Atran.num_col, -1);
  assert(nonactiveconstraintsidx.size() + activeconstraintidx.size() ==
         Atran.num_row);

  basisfactor.build();

  for (size_t i = 0;
       i < activeconstraintidx.size() + nonactiveconstraintsidx.size(); i++) {
    constraintindexinbasisfactor[baseindex[i]] = i;
  }
}

void Basis::report() {
  printf("basis: ");
  for (HighsInt a_ : activeconstraintidx) {
    printf("%" HIGHSINT_FORMAT " ", a_);
  }
  printf(" - ");
  for (HighsInt n_ : nonactiveconstraintsidx) {
    printf("%" HIGHSINT_FORMAT " ", n_);
  }
  printf("\n");
}

// move that constraint into V section basis (will correspond to Nullspace
// from now on)
void Basis::deactivate(HighsInt conid) {
  // printf("deact %" HIGHSINT_FORMAT "\n", conid);
  assert(contains(activeconstraintidx, conid));
  basisstatus.erase(conid);
  remove(activeconstraintidx, conid);
  nonactiveconstraintsidx.push_back(conid);
}

void Basis::activate(Runtime& rt, HighsInt conid, BasisStatus atlower,
                     HighsInt nonactivetoremove, Pricing* pricing) {
  // printf("activ %" HIGHSINT_FORMAT "\n", conid);
  if (!contains(activeconstraintidx, (HighsInt)conid)) {
    basisstatus[conid] = atlower;
    activeconstraintidx.push_back(conid);
  } else {
    printf("Degeneracy? constraint %" HIGHSINT_FORMAT
           " already in basis\n",
           conid);
    exit(1);
  }

  // printf("drop %d\n", nonactivetoremove);
  // remove non-active row from basis
  HighsInt rowtoremove = constraintindexinbasisfactor[nonactivetoremove];

  baseindex[rowtoremove] = conid;
  remove(nonactiveconstraintsidx, nonactivetoremove);
  updatebasis(rt, conid, nonactivetoremove, pricing);

  if (updatessinceinvert != 0) {
    constraintindexinbasisfactor[nonactivetoremove] = -1;
    constraintindexinbasisfactor[conid] = rowtoremove;
  }
}

void Basis::updatebasis(Runtime& rt, HighsInt newactivecon, HighsInt droppedcon,
                        Pricing* pricing) {
  if (newactivecon == droppedcon) {
    return;
  }

  HighsInt droppedcon_rowindex = constraintindexinbasisfactor[droppedcon];
  Atran.extractcol(newactivecon, buffer_column_aq);
  // column.report("col_pre_ftran");

  HVector column_aq_hvec = vec2hvec(buffer_column_aq);
  basisfactor.ftran(column_aq_hvec, 1.0);
  // column.report("col_post_ftran");

  Vector::unit(rt.instance.A.mat.num_col, droppedcon_rowindex, buffer_row_ep);
  // row_ep.report("rowep_pre_btran");

  HVector row_ep_hvec = vec2hvec(buffer_row_ep);
  basisfactor.btran(row_ep_hvec, 1.0);
  // row_ep.report("rowep_post_btran");

  pricing->update_weights(hvec2vec(column_aq_hvec), hvec2vec(row_ep_hvec),
                          droppedcon, newactivecon);

  HighsInt hint = 99999;
  HighsInt row_out = droppedcon_rowindex;

  updatessinceinvert++;
  basisfactor.update(&column_aq_hvec, &row_ep_hvec, &row_out, &hint);
  if (updatessinceinvert >= rt.settings.reinvertfrequency || hint != 99999) {
    rebuild();
    // printf("Hint: %d\n", hint);
    // printf("reinvert\n");
  }
}

Vector& Basis::btran(const Vector& rhs, Vector& target) const {
  HVector rhs_hvec = vec2hvec(rhs);
  basisfactor.btran(rhs_hvec, 1.0);
  return hvec2vec(rhs_hvec, target);
}

Vector Basis::btran(const Vector& rhs) const {
  HVector rhs_hvec = vec2hvec(rhs);
  basisfactor.btran(rhs_hvec, 1.0);

  return hvec2vec(rhs_hvec);
}

Vector& Basis::ftran(const Vector& rhs, Vector& target) const {
  HVector rhs_hvec = vec2hvec(rhs);
  basisfactor.ftran(rhs_hvec, 1.0);

  return hvec2vec(rhs_hvec, target);
}

Vector Basis::ftran(const Vector& rhs) const {
  HVector rhs_hvec = vec2hvec(rhs);
  basisfactor.ftran(rhs_hvec, 1.0);

  return hvec2vec(rhs_hvec);
}

Vector Basis::recomputex(const Instance& inst) {
  assert(activeconstraintidx.size() == inst.num_var);
  Vector rhs(inst.num_var);

  for (HighsInt i = 0; i < inst.num_var; i++) {
    HighsInt con = activeconstraintidx[i];
    if (constraintindexinbasisfactor[con] == -1) {
      printf("error\n");
    }
    if (basisstatus[con] == BasisStatus::ActiveAtLower) {
      if (con < inst.num_con) {
        rhs.value[constraintindexinbasisfactor[con]] = inst.con_lo[con];
      } else {
        rhs.value[constraintindexinbasisfactor[con]] =
            inst.var_lo[con - inst.num_con];
      }
    } else {
      if (con < inst.num_con) {
        rhs.value[constraintindexinbasisfactor[con]] = inst.con_up[con];
        // rhs.value[i] = inst.con_up[con];
      } else {
        rhs.value[constraintindexinbasisfactor[con]] =
            inst.var_up[con - inst.num_con];
        // rhs.value[i] = inst.var_up[con - inst.num_con];
      }
    }

    rhs.index[i] = i;
    rhs.num_nz++;
  }
  HVector rhs_hvec = vec2hvec(rhs);
  basisfactor.btran(rhs_hvec, 1.0);

  return hvec2vec(rhs_hvec);
}

// void Basis::write(std::string filename) {
//    FILE* file = fopen(filename.c_str(), "w");

//    fprintf(file, "%lu %lu\n", activeconstraintidx.size(),
//    nonactiveconstraintsidx.size()); for (HighsInt i=0;
//    i<activeconstraintidx.size(); i++) {
//       fprintf(file, "%" HIGHSINT_FORMAT " %" HIGHSINT_FORMAT "\n",
//       activeconstraintidx[i], (HighsInt)rowstatus[i]);
//    }
//    for (HighsInt i=0; i<nonactiveconstraintsidx.size(); i++) {
//       fprintf(file, "%" HIGHSINT_FORMAT " %" HIGHSINT_FORMAT "\n",
//       nonactiveconstraintsidx[i], (HighsInt)rowstatus[i]);
//    }
//    // TODO

//    fclose(file);
// }
