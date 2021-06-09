#include "basis.hpp"

#include <cassert>
#include <memory>

Basis::Basis(Runtime& rt, std::vector<int> active, std::vector<BasisStatus> lower, std::vector<int> inactive) : runtime(rt), buffer_column_aq(rt.instance.num_var), buffer_row_ep(rt.instance.num_var) {
   for (int i=0; i<active.size(); i++) {
      activeconstraintidx.push_back(active[i]);
      basisstatus[activeconstraintidx[i]] = lower[i];
   }
   for (int i : inactive) {
      nonactiveconstraintsidx.push_back(i);
   }

   Atran = rt.instance.A.t();

   build();
}

void Basis::build() {
   updatessinceinvert = 0;

   baseindex = new int[activeconstraintidx.size() + nonactiveconstraintsidx.size()];
   constraintindexinbasisfactor.clear();

   basisfactor = HFactor();
   
   constraintindexinbasisfactor.assign(Atran.num_row + Atran.num_col, -1);
   assert(nonactiveconstraintsidx.size() + activeconstraintidx.size() == Atran.num_row);
   
   int counter = 0;
   for (int i : nonactiveconstraintsidx) {
      baseindex[counter++] = i;
   }
   for (int i : activeconstraintidx) {
      baseindex[counter++] = i;
   }

   basisfactor.setup(Atran.num_col, Atran.num_row, (int*)&Atran.start[0],
                     (int*)&Atran.index[0], (const double*)&Atran.value[0],
                     baseindex);
   basisfactor.build();

   for (size_t i = 0; i < activeconstraintidx.size() + nonactiveconstraintsidx.size(); i++) {
      constraintindexinbasisfactor[baseindex[i]] = i;
   }
}

void Basis::rebuild() {
   updatessinceinvert = 0;
   constraintindexinbasisfactor.clear();
   
   constraintindexinbasisfactor.assign(Atran.num_row + Atran.num_col, -1);
   assert(nonactiveconstraintsidx.size() + activeconstraintidx.size() == Atran.num_row);

   basisfactor.build();

   for (size_t i = 0; i < activeconstraintidx.size() + nonactiveconstraintsidx.size(); i++) {
      constraintindexinbasisfactor[baseindex[i]] = i;
   }
}

void Basis::report() {
   printf("basis: ");
   for (int a_ : activeconstraintidx) {
      printf("%u ", a_);
   }
   printf(" - ");
   for (int n_ : nonactiveconstraintsidx) {
      printf("%u ", n_);
   }
   printf("\n");
}

   // move that constraint into V section basis (will correspond to Nullspace from
   // now on)
void Basis::deactivate(int conid) {
   // printf("deact %u\n", conid);
   assert(contains(activeconstraintidx, conid));
   basisstatus.erase(conid);
   remove(activeconstraintidx, conid);
   nonactiveconstraintsidx.push_back(conid);
}

void Basis::activate(Runtime& rt, int conid,
                        BasisStatus atlower, int nonactivetoremove, Pricing* pricing) {
   // printf("activ %u\n", conid);
   if (!contains(activeconstraintidx, (int)conid)) {
      basisstatus[conid] = atlower;
      activeconstraintidx.push_back(conid);
   } else {
      printf("Degeneracy? constraint %u already in basis\n", conid);
      exit(1);
   }

   // printf("drop %d\n", nonactivetoremove);
   // remove non-active row from basis
   int rowtoremove = constraintindexinbasisfactor[nonactivetoremove];
      
   baseindex[rowtoremove] = conid;
   remove(nonactiveconstraintsidx, nonactivetoremove);
   updatebasis(rt, conid, nonactivetoremove, pricing);

   if (updatessinceinvert != 0) {
      constraintindexinbasisfactor[nonactivetoremove] = -1;
      constraintindexinbasisfactor[conid] = rowtoremove;
   } 
}

void Basis::updatebasis(Runtime& rt, int newactivecon, int droppedcon, Pricing* pricing) {
   if (newactivecon == droppedcon) {
      return;
   }

   int droppedcon_rowindex = constraintindexinbasisfactor[droppedcon];
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

   pricing->update_weights(hvec2vec(column_aq_hvec), hvec2vec(row_ep_hvec), droppedcon, newactivecon);

   int hint = 99999;
   int row_out = droppedcon_rowindex;

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

   for (int i=0; i<inst.num_var; i++) {
      int con = activeconstraintidx[i];
      if (constraintindexinbasisfactor[con] == -1) {
         printf("error\n");
      }
      if (basisstatus[con] == BasisStatus::ActiveAtLower) {
         if (con < inst.num_con) {
            rhs.value[constraintindexinbasisfactor[con]] = inst.con_lo[con];
         } else {
            rhs.value[constraintindexinbasisfactor[con]] = inst.var_lo[con - inst.num_con];
         } 
      } else {
         if (con < inst.num_con) {
            rhs.value[constraintindexinbasisfactor[con]] = inst.con_up[con];
            // rhs.value[i] = inst.con_up[con];
         } else {
            rhs.value[constraintindexinbasisfactor[con]] = inst.var_up[con - inst.num_con];
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

   //    fprintf(file, "%lu %lu\n", activeconstraintidx.size(), nonactiveconstraintsidx.size());
   //    for (int i=0; i<activeconstraintidx.size(); i++) {
   //       fprintf(file, "%u %u\n", activeconstraintidx[i], (int)rowstatus[i]);
   //    }
   //    for (int i=0; i<nonactiveconstraintsidx.size(); i++) {
   //       fprintf(file, "%u %u\n", nonactiveconstraintsidx[i], (int)rowstatus[i]);
   //    }
   //    // TODO
      
   //    fclose(file);
   // }
