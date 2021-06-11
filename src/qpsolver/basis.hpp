#pragma once
#ifndef __SRC_LIB_BASIS_HPP__
#define __SRC_LIB_BASIS_HPP__

#include <map>
#include <vector>

#include "simplex/HFactor.h"

#include "pricing/pricing.hpp"
#include "runtime.hpp"
#include "instance.hpp"
#include "snippets.hpp"

#include "simplex/HVector.h"

inline HVector vec2hvec(const Vector& vec) {
   HVector hvec;
   hvec.setup(vec.dim);
   for (HighsInt i = 0; i < vec.num_nz; i++) {
      hvec.index[i] = vec.index[i];
      hvec.array[vec.index[i]] = vec.value[vec.index[i]];
   }
   hvec.count = vec.num_nz;
   hvec.packFlag = true;
   return hvec;
}

inline Vector& hvec2vec(const HVector& hvec, Vector& target) {
   target.reset();
   for (HighsInt i = 0; i < hvec.count; i++) {
      target.index[i] = hvec.index[i];
      target.value[target.index[i]] = hvec.array[hvec.index[i]];
   }
   target.num_nz = hvec.count;
   return target;
}

inline Vector hvec2vec(const HVector& hvec) {
   Vector vec(hvec.size);

   return hvec2vec(hvec, vec);
}

enum class BasisStatus {
   Default,
   ActiveAtLower = 1,
   ActiveAtUpper,
   ActiveAtZero,
   Inactive
};

class Basis {
   Runtime& runtime;
   HFactor basisfactor;
   HighsInt updatessinceinvert = 0;

   MatrixBase Atran;

   // indices of active constraints in basis
	std::vector<int> activeconstraintidx;

	// ids of constraints that are in the basis but not active
	// I need to extract those columns to get Z
	std::vector<int> nonactiveconstraintsidx;

	// ids of constraints that are in the basis 
	// std::vector<int> baseindex;
   int* baseindex;

   std::map<int, BasisStatus> basisstatus;

   // index i: -1 if constraHighsInt not in basis, [0, num_var] if constraHighsInt in basis (active or not)
	std::vector<int> constraintindexinbasisfactor;

   void build();
   void rebuild();

   // buffer to avoid recreating vectors
   Vector buffer_column_aq;
   Vector buffer_row_ep;

public:
   Basis(Runtime& rt, std::vector<int> active, std::vector<BasisStatus> atlower, std::vector<int> inactive);

   HighsInt getnupdatessinceinvert() {
      return updatessinceinvert;
   }

   HighsInt getnumactive() const { return activeconstraintidx.size(); };

   HighsInt getnuminactive() const { return nonactiveconstraintsidx.size(); };

   const std::vector<int>& getactive() const { return activeconstraintidx; };

   const std::vector<int>& getinactive() const { return nonactiveconstraintsidx; };

   const std::vector<int>& getindexinfactor() const { return constraintindexinbasisfactor; };

   BasisStatus getstatus(HighsInt conid) { return basisstatus[conid]; };

   void report();

   // move that constraHighsInt into V section basis (will correspond to Nullspace from now on)
   void deactivate(HighsInt conid);

   void activate(Runtime& rt, HighsInt conid, BasisStatus atlower, HighsInt nonactivetoremove, Pricing* pricing);

   void updatebasis(Runtime& rt, HighsInt newactivecon, HighsInt droppedcon, Pricing* pricing);

   Vector btran(const Vector& rhs) const;

   Vector ftran(const Vector& rhs) const;

   Vector& btran(const Vector& rhs, Vector& target) const;

   Vector& ftran(const Vector& rhs, Vector& target) const;

   Vector recomputex(const Instance& inst);

   void write(std::string filename);
};

#endif
