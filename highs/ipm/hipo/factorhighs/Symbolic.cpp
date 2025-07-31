#include "Symbolic.h"

#include <iostream>

#include "FactorHiGHSSettings.h"
#include "ipm/hipo/auxiliary/Log.h"

namespace hipo {

Symbolic::Symbolic() : block_size_{kBlockSize} {}

void Symbolic::setParallel(bool par_tree, bool par_node) {
  parallel_tree_ = par_tree;
  parallel_node_ = par_node;
}

int64_t Symbolic::nz() const { return nz_; }
double Symbolic::flops() const { return flops_; }
double Symbolic::spops() const { return spops_; }
double Symbolic::critops() const { return critops_; }
Int Symbolic::blockSize() const { return block_size_; }
Int Symbolic::size() const { return n_; }
Int Symbolic::sn() const { return sn_; }
Int Symbolic::rows(Int i) const { return rows_[i]; }
Int Symbolic::ptr(Int i) const { return ptr_[i]; }
Int Symbolic::snStart(Int i) const { return sn_start_[i]; }
Int Symbolic::snParent(Int i) const { return sn_parent_[i]; }
Int Symbolic::relindCols(Int i) const { return relind_cols_[i]; }
Int Symbolic::relindClique(Int i, Int j) const { return relind_clique_[i][j]; }
Int Symbolic::consecutiveSums(Int i, Int j) const {
  return consecutive_sums_[i][j];
}
Int Symbolic::cliqueBlockStart(Int sn, Int bl) const {
  return clique_block_start_[sn][bl];
}
Int Symbolic::cliqueSize(Int sn) const {
  return clique_block_start_[sn].back();
}
bool Symbolic::parTree() const { return parallel_tree_; }
bool Symbolic::parNode() const { return parallel_node_; }

const std::vector<Int>& Symbolic::ptr() const { return ptr_; }
const std::vector<Int>& Symbolic::iperm() const { return iperm_; }
const std::vector<Int>& Symbolic::snParent() const { return sn_parent_; }
const std::vector<Int>& Symbolic::snStart() const { return sn_start_; }
const std::vector<Int>& Symbolic::pivotSign() const { return pivot_sign_; }

std::string memoryString(double mem) {
  std::stringstream ss;

  if (mem < 1024)
    ss << sci(mem, 0, 1) << " B";
  else if (mem < 1024 * 1024)
    ss << sci(mem / 1024, 0, 1) << " KB";
  else if (mem < 1024 * 1024 * 1024)
    ss << sci(mem / 1024 / 1024, 0, 1) << " MB";
  else
    ss << sci(mem / 1024 / 1024 / 1024, 0, 1) << " GB";

  return ss.str();
}

void Symbolic::print(const Log* log, bool verbose) const {
  std::stringstream log_stream;
  log_stream << "\nFactorisation statistics\n";
  log_stream << textline("Size:") << sci(n_, 0, 2) << '\n';
  log_stream << textline("Nnz:") << sci(nz_, 0, 2) << '\n';
  log_stream << textline("Fill-in:") << fix(fillin_, 0, 2) << '\n';
  log_stream << textline("Serial memory:") << memoryString(serial_storage_)
             << '\n';
  log_stream << textline("Flops:") << sci(flops_, 0, 1) << '\n';
  if (verbose) {
    log_stream << textline("Sparse ops:") << sci(spops_, 0, 1) << '\n';
    log_stream << textline("critical ops:") << sci(critops_, 0, 1) << '\n';
    log_stream << textline("Max tree speedup:") << fix(flops_ / critops_, 0, 2)
               << '\n';
    log_stream << textline("Artificial nz:") << sci(artificial_nz_, 0, 1)
               << '\n';
    log_stream << textline("Artificial ops:") << sci(artificial_ops_, 0, 1)
               << '\n';
    log_stream << textline("Largest front:") << integer(largest_front_, 0)
               << '\n';
    log_stream << textline("Largest supernode:") << integer(largest_sn_, 0)
               << '\n';
    log_stream << textline("Supernodes:") << integer(sn_, 0) << '\n';
    log_stream << textline("Sn size <= 1:") << integer(sn_size_1_, 0) << '\n';
    log_stream << textline("Sn size <= 10:") << integer(sn_size_10_, 0) << '\n';
    log_stream << textline("Sn size <= 100:") << integer(sn_size_100_, 0)
               << '\n';
    log_stream << textline("Sn avg size:") << sci(n_, 0, 1) << '\n';
  }

  log_stream << '\n';
  if (log) log->print(log_stream);
}

}  // namespace hipo