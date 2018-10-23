#ifndef HDUALRHS_H_
#define HDUALRHS_H_

#include <vector>
using namespace std;

#include "HModel.h"
#include "HVector.h"

/**
 * This class deal with optimality test
 * and some update primal/weight tasks
 */

class HDualRHS {
 public:
  void setup(HModel *model);
  void setup_partition(const char *filename);

  void choose_normal(int *chIndex);
  void choose_multi_global(int *chIndex, int *chCount, int chLimit);
  void choose_multi_HGauto(int *chIndex, int *chCount, int chLimit);
  void choose_multi_HGpart(int *chIndex, int *chCount, int chLimit);

  void update_primal(HVector *column, double theta);
  void update_weight(HVector *column, double devex, double Kai, double *dse);
  void update_weight_Dvx(HVector *column, double dvx_wt_o_rowOut);
  void update_pivots(int iRow, double value);

  void update_infeasList(HVector *column);
  void create_infeasList(double columnDensity);
  void create_infeasArray();

  HModel *workModel;
  double workCutoff;
  int workCount;
  vector<char> workMark;
  vector<int> workIndex;
  vector<double> workArray;
  vector<double> workEdWt;
  vector<double> workEdWtFull;

  int partNum;
  int partNumRow;
  int partNumCol;
  int partNumCut;
  int partSwitch;
  vector<int> workPartition;
};

#endif /* HDUALRHS_H_ */
