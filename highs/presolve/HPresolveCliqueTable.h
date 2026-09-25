/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HPRESOLVE_CLIQUE_TABLE_H_
#define HPRESOLVE_CLIQUE_TABLE_H_

#include <vector>

#include "mip/HighsCliqueTable.h"

class HPresolveCliqueTable {
  HighsCliqueTable* table = nullptr;

  enum class ColState : int8_t { kActive, kFixedZero, kFixedOne, kEliminated };

  std::vector<ColState> colStates;
  std::vector<HighsCliqueTable::CliqueVar> fixingQueue;
  std::vector<HighsInt> incidentCliques;
  std::vector<HighsCliqueTable::CliqueVar> shortenedClique;

  void checkCompactClique(HighsInt cliqueId, HighsInt threshold,
                          HighsInt activeSize, HighsInt actualSize,
                          bool equality, HighsInt origin);

 public:
  void setCliqueTable(HighsCliqueTable* table) { this->table = table; }

  void rebuild(HighsCliqueTable& table);

  bool fixCol(HighsInt col, bool val,
              std::vector<HighsCliqueTable::CliqueVar>& impliedFixings);

  void eliminateCol(HighsInt col);

  bool substituteCol(HighsInt substCol, HighsCliqueTable::CliqueVar replacement,
                     std::vector<HighsCliqueTable::CliqueVar>& impliedFixings);
};

#endif
