void assignBasis() {
  const int numTot = model.lp_.numCol_ + model.lp_.numRow_;
  model.basis_.basicIndex_.resize(model.lp_.numRow_);
  model.basis_.nonbasicFlag_.assign(numTot, 0);
  model.basis_.nonbasicMove_.resize(numTot);
}

