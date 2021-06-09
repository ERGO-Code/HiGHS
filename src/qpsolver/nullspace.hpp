#ifndef __SRC_LIB_NULLSPACE_HPP__
#define __SRC_LIB_NULLSPACE_HPP__

#include "matrix.hpp"
#include "basis.hpp"

#include "runtime.hpp"

struct NullspaceReductionResult {
   int maxabsd;
   int constrainttodrop;
   Vector& d;

   bool p_in_v;

   NullspaceReductionResult(bool pinv, int mad, Vector& d_, int ctd) : maxabsd(mad), constrainttodrop(ctd), d(d_), p_in_v(pinv) {;
   }
};

class Nullspace {
   bool uptodateZ = false;

   const Basis& basis;
   Runtime& runtime;
   Matrix bufferZ;

   Vector temp_unit;
   Vector buffer_d;
   Vector buffer_col;

   Vector buffer_aq;
   Vector buffer_col_p;

   Vector buffer_unit;

   void recompute() {
      int nvar = basis.getnumactive() + basis.getnuminactive();
      Matrix Z(nvar, 0);

      const std::vector<int>& nonactive = basis.getinactive(); 
      const std::vector<int>& indexinfactor = basis.getindexinfactor();
      
      for (int i = 0; i < nonactive.size(); i++) {
         int unit =
            indexinfactor[nonactive[i]];

         Vector::unit(nvar, unit, temp_unit);

         basis.btran(temp_unit, temp_unit);

         Z.append(temp_unit);
      }
      bufferZ = Z;
      uptodateZ = true;
   }

   Vector& aq_Z_prod(Runtime& rt, int q, Vector& target) {
      target.reset();
      Matrix& Z = getNullspace();

      if (q >= rt.instance.num_con) {
         // Vector aq = rt.instance.A.t().extractcol(q);
         // return Z.vec_mat(aq);
         int ep = q - rt.instance.num_con;
         // 
         for (int col=0; col<Z.mat.num_col; col++) {
            double dot = 0.0;
            for (int idx=Z.mat.start[col]; idx<Z.mat.start[col+1]; idx++) {
               if (Z.mat.index[idx] == ep) {
                  dot += Z.mat.value[idx];
                  break;
               }
            }

            if (dot != 0.0) {
               target.value[col] = dot;
               target.index[target.num_nz] = col;
               target.num_nz++;
            }
         }
         return target;
      } else {
         rt.instance.A.t().extractcol(q, buffer_aq);
         return Z.vec_mat(buffer_aq, target);
         // MatrixBase& Atran = rt.instance.A.t();
         // return Z.vec_mat(&Atran.index[Atran.start[q]], &Atran.value[Atran.start[q]], Atran.start[q+1] - Atran.start[q]);
         // // return res2;
      }
   }

public:
   Nullspace(Runtime& rt, Basis& bas) : basis(bas), runtime(rt), bufferZ(Matrix(rt.instance.num_var,0)), temp_unit(rt.instance.num_var), buffer_d(rt.instance.num_var), buffer_col(rt.instance.num_var), buffer_aq(rt.instance.num_var), buffer_col_p(rt.instance.num_var), buffer_unit(rt.instance.num_var) {
      if (bas.getnumactive() > 0) {
         recompute();
      }
      uptodateZ = true;
   }

   Vector& expand_computenewcol(int conid, Vector& target) {
      if (uptodateZ) {
         int unit = basis.getindexinfactor()[conid];
         Vector::unit(bufferZ.mat.num_row, unit, target);

         basis.btran(target, target);

         return target;
      } 
      exit(1);
   }

   void expand_appendnewcol(Vector& newcol) {
      if (uptodateZ) {
         bufferZ.append(newcol);
      } 
   }

   NullspaceReductionResult reduce(Runtime& rt, int newactivecon) { 
      Matrix& Z = getNullspace();

      int idx = indexof(basis.getinactive(), newactivecon);
      if (idx != -1) {
         bufferZ.dropcol(idx);
         buffer_unit.dim = Z.mat.num_col+1;
         Vector::unit(Z.mat.num_col, idx, buffer_unit);
         return NullspaceReductionResult(true, idx, buffer_unit, newactivecon);
      }

      // TODO: this operation is inefficient.
      // Vector aq = rt.instance.A.t().extractcol(newactivecon);
      // Vector d = Z.vec_mat(aq);
      aq_Z_prod(rt, newactivecon, buffer_d);

      int maxabs = 0;
      for (int i = 0; i < buffer_d.num_nz; i++) {
         if (fabs(buffer_d.value[buffer_d.index[i]]) > fabs(buffer_d.value[maxabs])) {
            maxabs = buffer_d.index[i];
         }
      }

      if (fabs(buffer_d.value[maxabs]) < rt.settings.d_zero_threshold) {
         printf("degeneracy? not possible to find non-active constraint to leave basis. max: log(d[%u]) = %lf\n", maxabs, log10(fabs(buffer_d.value[maxabs])));
         exit(1);
      }

      Matrix Z_(Z.mat.num_row, 0);

      Z.mat.extractcol(maxabs, buffer_col_p);
      // assert(col_p == row_ep);
      for (int i=0; i<Z.mat.num_col; i++) {
         if (i == maxabs) {
            continue;
         }
         if (buffer_d.value[i] != 0.0) {
            // TODO can this be done without creating too manu additional vectors?
            Z.mat.extractcol(i, buffer_col);
            buffer_col.saxpy(-buffer_d.value[i] / buffer_d.value[maxabs], buffer_col_p);
            Z_.append(buffer_col);
         } else {
            Z_.append(&Z.mat.index[Z.mat.start[i]], &Z.mat.value[Z.mat.start[i]], Z.mat.start[i+1] - Z.mat.start[i]);
         }
      }

      bufferZ = Z_;
      return NullspaceReductionResult(false, maxabs, buffer_d, basis.getinactive()[maxabs]);
   }

   Matrix& getNullspace() {
      if (!uptodateZ) {
         recompute();
      }
      return bufferZ;
   }

   // TODO: leads to cycling
   Vector& Ztprod(const Vector& rhs, Vector& target) {
      bool test = false;
      Matrix& Z = getNullspace();
      if (test) {
         Vector res_ = basis.ftran(rhs);
         
         target.reset();
         for (int i=0; i<Z.mat.num_col; i++) {
            int nonactive = basis.getinactive()[i];
            int idx = basis.getindexinfactor()[nonactive];
            target.index[i] = i;
            target.value[i] = res_.value[idx];
         }
         target.resparsify();
         return target;

         // x' = b * Z
         // x = Z' b
         // x = B^-1 b
         // B x = b
      } else {
         return Z.vec_mat(rhs, target);
      }
   }

   // TODO: leads to cycling
   Vector& Zprod(const Vector& rhs, Vector& target) {
      Matrix& Z = getNullspace();
      bool test = false;
      if (test) {
         for (int i=0; i<rhs.num_nz; i++) {
            int nonactive = basis.getinactive()[i];
            int idx = basis.getindexinfactor()[nonactive];
            target.index[i] = idx;
            target.value[idx] = rhs.value[i];
         }
         target.num_nz = rhs.num_nz;
         return basis.btran(target, target);
      } else {
         return Z.mat_vec(rhs, target);
      }
   }

   double density() {
      if (getNullspace().mat.value.size() > 0) {
         return (double)getNullspace().mat.value.size()/(getNullspace().mat.num_col * getNullspace().mat.num_row);
      } else {
         return 0.0;
      }
      
   }

};

#endif
