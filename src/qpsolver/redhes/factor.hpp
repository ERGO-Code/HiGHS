#ifndef __SRC_LIB_NEWFACTOR_HPP__
#define __SRC_LIB_NEWFACTOR_HPP__

#include <vector>
#include <cassert>
#include <cmath>

#include "../matrix.hpp"
#include "../nullspace.hpp"
#include "../runtime.hpp"

class NewCholeskyFactor {
private:
   bool uptodate = false;

   Runtime& runtime;
   Nullspace& nullspace;

   std::vector<std::vector<double>> orig;

   int current_k = 0;
   int current_k_max;
   std::vector<double> L;

   void recompute() {
      Matrix& Z = nullspace.getNullspace();
      // M = (Q' * Z)' * Z
      // Matrix m = Z.tran_mat(runtime.instance.Q).mat_mat(Z);
      Matrix m = Matrix(Matrix(Z.mat.tran_mat_(runtime.instance.Q.mat), false).t().tran_mat_(Z.mat), false);

      orig.assign(m.mat.num_col, std::vector<double>(m.mat.num_col, 0.0));
 
      for (int i=0; i<m.mat.num_col; i++) {
         for (int j=m.mat.start[i]; j<m.mat.start[i+1]; j++) {
            int row = m.mat.index[j];
            orig[row][i] = m.mat.value[j];
         }
      }

      for (size_t col = 0; col < orig.size(); col++) { 
        for (size_t row = 0; row <= col; row++) { 
            double sum = 0; 
            if (row == col) { 
               for (size_t k = 0; k < row; k++) 
                    sum += L[k * current_k_max + row] * L[k * current_k_max + row];
                L[row * current_k_max + row] = sqrt(orig[row][row] - sum); 
            } else { 
               for (size_t k = 0; k < row; k++) 
                  sum += (L[k * current_k_max + col] * L[k * current_k_max + row]); 
               L[row * current_k_max + col] = (orig[col][row] - sum) / L[row * current_k_max + row]; 
            } 
         } 
      } 
      current_k = Z.mat.num_col;
      uptodate = true;
   }

void resize() {
   std::vector<double> L_old = L;
   L.clear();
   L.resize( (current_k_max * 2) *  (current_k_max * 2));
   for (int i=0; i<current_k_max; i++) {
      for (int j=0; j<current_k_max; j++) {
         L[i * (current_k_max * 2) + j] = L_old[i * current_k_max + j];
      }
   }
   current_k_max *= 2;
}

public:
   NewCholeskyFactor(Runtime& rt, Basis& basis, Nullspace& ns) : runtime(rt), nullspace(ns) {
      uptodate = false;
      printf("computed size: %u, needed size: %u\n", min((int)ceil(rt.instance.num_var / 16.0), 1000), basis.getnuminactive());
      current_k_max = max(min((int)ceil(rt.instance.num_var / 16.0), 1000), basis.getnuminactive());
      L.resize(current_k_max * current_k_max);
   }

   void expand(Vector& yp, Vector& gyp, Vector& l) {
      if (!uptodate) {
         return;
      }
      double mu = gyp * yp;
      l.resparsify();
      double lambda = mu - l.norm2();

      assert(lambda > 0);
      
      if (current_k_max == current_k+1) {
         resize();
      }

      for (int i=0; i<current_k; i++) {
         L[i *current_k_max + current_k] = l.value[i];
      }
      L[current_k * current_k_max + current_k] = sqrt(lambda);

      current_k++;
   }


   void solveL(Vector& rhs) {
      if (!uptodate) {
         recompute();
      }

      for (int r=0; r<rhs.dim; r++) {

         for (int j=0; j<r; j++) {
            rhs.value[r] -= rhs.value[j]*L[j * current_k_max + r];
         }

         rhs.value[r] /=  L[r * current_k_max + r];
      }
   }

   // solve L' u = v
   void solveLT(Vector& rhs) {
      for (int i=rhs.dim-1; i>=0; i--) {
         double sum = 0.0;
         for (int j=rhs.dim-1; j>i; j--) {
            sum += rhs.value[j] * L[i* current_k_max + j];
         } 
         rhs.value[i] = (rhs.value[i] - sum) / L[i* current_k_max + i];
      }
   }


public:
   void solve(Vector& rhs) {
      if (!uptodate) {
         recompute();
      }
      solveL(rhs);
      solveLT(rhs);

      rhs.resparsify();
   }

   void eliminate(std::vector<double>& m, int i, int j, int kmax, int currentk) {
      // i = col, j = row
      if (m[j * kmax +i] == 0.0) {
         return;
      }
      double z = sqrt(m[i * kmax +i] *  m[i * kmax +i] +  m[j * kmax +i] *  m[j * kmax +i]);
      double cos_, sin_;
      if (z == 0) {
         cos_ = 1.0;
         sin_ = 0.0;
      } else {
         cos_ = m[i * kmax +i] / z;
         sin_ = - m[j * kmax +i] / z;
      }
      
      if (sin_ == 0.0) {
         if (cos_ > 0.0) {
            // nothing
         } else {
            for (int k=0; k<current_k; k++) {
               // update entry i and j of column k
               double a_ik = m[i*kmax + k];
               // entry i
               m[i*kmax+k] = -a_ik;
               m[j*kmax+k] = -m[j*kmax+k];
            }
         }
      } else if (cos_ == 0.0) {
         if (sin_ > 0.0) {
            for (int k=0; k<current_k; k++) {
               // update entry i and j of column k
               double a_ik = m[i*kmax + k];
               // entry i
               m[i*kmax+k] = - m[j*kmax+k];
               m[j*kmax+k] = a_ik;
            }
         } else {
            for (int k=0; k<current_k; k++) {
               // update entry i and j of column k
               double a_ik = m[i*kmax + k];
               // entry i
               m[i*kmax+k] = m[j*kmax+k];
               m[j*kmax+k] = -a_ik;
            }
         }
      } else {
         // #pragma omp parallel for
         for (int k=0; k<current_k; k++) {
            // update entry i and j of column k
            double a_ik = m[i*kmax + k];
            // entry i
            m[i*kmax+k] = cos_ * a_ik - sin_ * m[j*kmax+k];
            m[j*kmax+k] = sin_ * a_ik + cos_ * m[j*kmax+k];
         }
      }
      m[j * kmax + i] = 0.0;
   }

   void reduce(NullspaceReductionResult& nrr) {
      if (current_k == 0) {
         return;
      }

      unsigned p = nrr.maxabsd; // col we push to the right and remove

      // start situation: p=3, current_k = 5
      // |1 x  | |x    |       |1   | |xxxxx|  
      // | 1x  | |xx   |  ===  | 1  | | xxxx| 
      // |  x1 | |xxx  |       |xxxx| |  xxx|
      // |  x 1| |xxxx |       |  1 | |   xx|
      //         |xxxxx|       |   1| |    x|
      // next step: move row/col p to the bottom/right

      //> save row p
      std::vector<double> row_p(current_k, 0.0);
      for (int i=0; i<current_k; i++) {
         row_p[i] = L[p*current_k_max + i];
      }

      //> move all rows > p up by one row 
      for (int row=p; row<current_k-1; row++) {
         for (int i=0; i<current_k; i++) {
            L[row*current_k_max + i] = L[(row+1) * current_k_max + i];
         }
      }

      //> load row p
      for (int i=0; i<current_k; i++) {
         L[(current_k-1) * current_k_max + i] = row_p[i];
      }

      //> now move col p to the right in each row
      for (int row=0; row<current_k; row++) {
         double p_entry = L[row * current_k_max + p];
         for (int col=p; col<current_k-1; col++) {
            L[row*current_k_max+col] = L[row*current_k_max+col+1];
         }
         L[row * current_k_max + current_k-1] = p_entry;
      }

      if (current_k == 1) {
         current_k--;
         return;
      }
      
      if (!nrr.p_in_v) {
         // situation now:
         // |1   x| |x    |       |1   | |xxxxx|  
         // | 1  x| |xx   |  ===  | 1  | | xxxx| 
         // |  1 x| |xxx x|       |  1 | |  xx |
         // |   1x| |xxxxx|       |   1| |   x |
         //         |xx  x|       |xxxx| |  xxx|
         // next: remove nonzero entries in last column except for diagonal element
         for (int r=p-1; r>=0; r--) { // to current_k-1
            eliminate(L, current_k-1, r, current_k_max, current_k);
         }

         // situation now:
         // |1   x| |x   x|        |xxxx | |1   | 
         // | 1  x| |xx  x|  ===   | xxx | | 1  |
         // |  1 x| |xxx x|        |  xx | |  1 |
         // |   1x| |xxxxx|        |   x | |   1|
         //         |    x|        |xxxxx| |xxxx|
         // next: multiply product
         // new last row: old last row (first current_k-1 elements) + r * R_current_k_current_k

         for (int i=0; i<nrr.d.num_nz; i++) {
            int idx = nrr.d.index[i];
            if (idx == nrr.maxabsd) {
               continue;
            }
            if (idx < nrr.maxabsd) {
               L[(current_k-1) * current_k_max + idx] += -nrr.d.value[idx] / nrr.d.value[nrr.maxabsd] * L[(current_k-1) * current_k_max + current_k-1];
            } else {
               L[(current_k-1) * current_k_max + idx-1] += -nrr.d.value[idx] / nrr.d.value[nrr.maxabsd] * L[(current_k-1) * current_k_max + current_k-1];
            }
         }
         // situation now: as above, but no more product
      }
      // next: eliminate last row
      for (int i=0; i<current_k-1; i++) {
         eliminate(L, i, current_k-1, current_k_max, current_k);
      }
      current_k--;
   }

   void report(std::string name = "") {
      printf("%s\n", name.c_str());
      for (int i=0; i<current_k; i++) {
         for (int j=0; j<current_k; j++) {
            printf("%lf ", L[i * current_k_max + j]);
         }
         printf("\n");
      }
   }

   double density() {
      if (current_k == 0) {
         return 0.0;
      }
      
      int num_nz = 0;
      for (int i=0; i<current_k; i++) {
         for (int j=0; j<current_k; j++) {
            if (fabs(L[i * current_k_max + j]) > 10e-8) {
               num_nz++;
            }
         }
      }
      return (double)num_nz/(current_k * (current_k+1) / 2.0);
   }

};

#endif
