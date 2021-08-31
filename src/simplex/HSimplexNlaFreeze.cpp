/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HSimplexNla.cpp
 *
 * @brief Interface to HFactor allowing non-HFactor updates, NLA-only
 * scaling and shifting of NLA analysis below simplex level.
 */
#include "simplex/HSimplexNla.h"

#include <stdio.h>

void FrozenBasis::clear() {
  this->valid_ = false;
  this->clearAllButBasis();
  this->basis_.clear();
}

void FrozenBasis::clearAllButBasis() {
  this->prev_ = kNoLink;
  this->next_ = kNoLink;
  this->update_.clear();
}

void HSimplexNla::frozenBasisClearAllData() {
  this->first_frozen_basis_id_ = kNoLink;
  this->last_frozen_basis_id_ = kNoLink;
  this->frozen_basis_.clear();
  this->update_.clear();
}

void HSimplexNla::frozenBasisClearAllButBasis() {
  this->first_frozen_basis_id_ = kNoLink;
  this->last_frozen_basis_id_ = kNoLink;
  for (HighsInt frozen_basis_id = 0; frozen_basis_id < this->frozen_basis_.size(); frozen_basis_id++) 
    this->frozen_basis_[frozen_basis_id].clearAllButBasis();
  this->update_.clear();
}

bool HSimplexNla::frozenBasisIdValid(const HighsInt frozen_basis_id) const {
  bool valid_id = 0<=frozen_basis_id && frozen_basis_id < frozen_basis_.size();
  if (valid_id) valid_id = frozen_basis_[frozen_basis_id].valid_;
  return valid_id;
}

bool HSimplexNla::frozenBasisHasInvert(const HighsInt frozen_basis_id) const {
  return this->last_frozen_basis_id_ != kNoLink;
}

HighsInt HSimplexNla::freeze(const SimplexBasis& basis, const double col_aq_density) {
  frozen_basis_.push_back(FrozenBasis());
  HighsInt this_frozen_basis_id = frozen_basis_.size() - 1;
  FrozenBasis& frozen_basis = frozen_basis_[this_frozen_basis_id];
  frozen_basis.valid_ = false;
  frozen_basis.prev_ = last_frozen_basis_id_;
  frozen_basis.next_ = kNoLink;
  frozen_basis.update_.clear();
  frozen_basis.basis_ = basis;
  if (last_frozen_basis_id_ == kNoLink) {
    // There is thisly no frozen basis, so record this as the first
    first_frozen_basis_id_ = this_frozen_basis_id;
  } else {
    // Update the forward link from the previous last frozen basis
    FrozenBasis& frozen_basis = frozen_basis_[last_frozen_basis_id_];
    frozen_basis.next_ = this_frozen_basis_id;
    // The PF updates held in simplex NLA now become the updates to
    // apply in order to go from the previous frozen basis to the
    // latest one
    frozen_basis.update_ = std::move(update_);
  }
  last_frozen_basis_id_ = this_frozen_basis_id;
  update_.setup(lp_->num_row_, col_aq_density);
  return this_frozen_basis_id;
}

void HSimplexNla::unfreeze(const HighsInt unfreeze_basis_id, SimplexBasis& basis) {
  assert(frozenBasisIdValid(unfreeze_basis_id));
  FrozenBasis& frozen_basis = frozen_basis_[unfreeze_basis_id];
  // Move the frozen basis into the return basis
  basis = std::move(frozen_basis.basis_);
  // This frozen basis, and any linked forward from it, must be
  // cleared. Any link to it must also be cleared.
  HighsInt prev_frozen_basis_id = frozen_basis.prev_;
  HighsInt frozen_basis_id = unfreeze_basis_id;
  if (prev_frozen_basis_id == kNoLink) {
    // There is no previous frozen basis linking to this one, so all
    // frozen basis data can be cleared
    frozenBasisClearAllData();
  } else {
    // The previous frozen basis is now the last and has no link
    // forward
    this->last_frozen_basis_id_ = prev_frozen_basis_id;
    frozen_basis_[prev_frozen_basis_id].next_ = kNoLink;
    for(;;) {
      assert(frozen_basis_id != kNoLink);
      frozen_basis_id = frozen_basis_[frozen_basis_id].next_;
      if (frozen_basis_id == kNoLink) break;
      frozen_basis_[frozen_basis_id].clear();
    }
    // Clear any recent PF updates
    this->update_.clear();
  }
}

