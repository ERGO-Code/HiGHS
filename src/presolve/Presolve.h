/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/HPresolve.h
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef PRESOLVE_PRESOLVE_H_
#define PRESOLVE_PRESOLVE_H_

#include <list>
#include <stack>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "HFactor.h"
#include "HPreData.h"
#include "HTimerPre.h"
#include "KktChStep.h"
#include "HighsLp.h"

using namespace std;

enum class HighsPostsolveStatus {
  ReducedSolutionEmpty,
  ReducedSolutionDimenionsError,
  SolutionRecovered,
  LpOrPresolveObjectMissing,
  NoPostsolve
};

enum class HighsPresolveStatus {
    NotReduced,
    Infeasible,
    Unbounded,
    Empty,
    Reduced,
    ReducedToEmpty,
    NullError
};


class Presolve : public HPreData {
 public:
  Presolve();
  HighsPresolveStatus presolve();
  HighsPostsolveStatus postsolve();

  void load(const HighsLp& lp);
 // todo: clear the public from below. 
  string modelName;

 private: 

  int iPrint;
  int iKKTcheck;
  int presolve(int print);


  double objShift;
  void initializeVectors();
  void setProblemStatus(const int s);
  void reportTimes();

  friend class HModel;
  // new bounds on primal variables for implied free detection
  vector<double> implColLower;
  vector<double> implColUpper;
  vector<int> implColLowerRowIndex;
  vector<int> implColUpperRowIndex;

  vector<int> implRowDualLowerSingColRowIndex;
  vector<int> implRowDualUpperSingColRowIndex;

  // new bounds on row duals y_i
  vector<double> implRowDualLower;
  vector<double> implRowDualUpper;

  vector<double> implColDualLower;
  vector<double> implColDualUpper;
  vector<double> implRowValueLower;
  vector<double> implRowValueUpper;

  HTimerPre timer;  // holds enum for main presolve rules

  enum stat {
    Unset = 0,
    Infeasible = 1,
    Unbounded = 2,
    Empty = 3,
    Optimal = 4,
    Reduced = 5,
  };

 private:
  bool hasChange;
  int status = 0;  // 0 is unassigned, see enum stat
  friend class PresolveInfo;

  list<int> singRow;  // singleton rows
  list<int> singCol;  // singleton columns

  // original data
 public:
  vector<double> colCostOriginal;

 private:
  vector<double> rowLowerOriginal;
  vector<double> rowUpperOriginal;
  vector<double> colLowerOriginal;
  vector<double> colUpperOriginal;

  // functions
  void setPrimalValue(int j, double value);
  void checkForChanges(int iteration);
  void resizeProblem();
  void resizeImpliedBounds();

  // easy transformations
  void removeIfFixed(int j);
  void removeEmptyRow(int i);
  void removeEmptyColumn(int j);
  void removeRow(int i);

  // singleton rows
  void removeRowSingletons();
  int getSingRowElementIndexInAR(int i);
  int getSingColElementIndexInA(int j);

  // forcing constraints
  void removeForcingConstraints(int mainIter);
  pair<double, double> getImpliedRowBounds(int row);
  void setVariablesToBoundForForcingRow(const int row, const bool isLower);
  void dominatedConstraintProcedure(const int i, const double g,
                                    const double h);

  // doubleton equations
  void removeDoubletonEquations();
  pair<int, int> getXYDoubletonEquations(const int row);
  void processRowDoubletonEquation(const int row, const int x, const int y,
                                   const double akx, const double aky,
                                   const double b);
  pair<double, double> getNewBoundsDoubletonConstraint(int row, int col, int j,
                                                       double aik, double aij);
  void UpdateMatrixCoeffDoubletonEquationXzero(const int i, const int x,
                                               const int y, const double aiy,
                                               const double akx,
                                               const double aky);
  void UpdateMatrixCoeffDoubletonEquationXnonZero(const int i, const int x,
                                                  const int y, const double aiy,
                                                  const double akx,
                                                  const double aky);

  // column singletons
  void removeColumnSingletons();
  bool removeIfImpliedFree(int col, int i, int k);
  void removeFreeColumnSingleton(const int col, const int row, const int k);
  bool removeColumnSingletonInDoubletonInequality(const int col, const int i,
                                                  const int k);
  void removeSecondColumnSingletonInDoubletonRow(const int j, const int i);
  pair<double, double> getBoundsImpliedFree(double lowInit, double uppInit,
                                            const int col, const int i,
                                            const int k);
  void removeImpliedFreeColumn(const int col, const int i, const int k);

  // dominated columns
  void removeDominatedColumns();
  void rowDualBoundsDominatedColumns();
  pair<double, double> getImpliedColumnBounds(int j);
  void removeIfWeaklyDominated(const int j, const double d, const double e);

  //    void findDuplicateRows();
  //    void findDuplicateColumns();
  //    void removeDuplicateRows(int i, int k, double v);
  //    int makeCheckForDuplicateRows(int k, int i, vector<double>& coeff,
  //    vector<int>& colIndex, double v, int whichIsFirst); void
  //    removeDuplicateColumns(int j,int k, double v); bool
  //    checkDuplicateRows(int i, int k) ;
  //	  bool checkDuplicateColumns(int i, int k) ;

  // old or test
  // void updateRemovedColRow(int dim);
  // void updateRowsByNZ();
  void testAnAR(int post);

  vector<int> countRemovedRows;
  vector<int> countRemovedCols;
  double tol;

  // postsolve
  bool noPostSolve;

  void addChange(int type, int row, int col);
  void fillStackRowBounds(int col);
  void setKKTcheckerData();

  void getBoundOnLByZj(int row, int j, double *lo, double *up, double colLow,
                       double colUpp);
  double getRowDualPost(int row, int col);
  double getColumnDualPost(int col);
  string getDualsForcingRow(int row, vector<int> &fRjs);
  void getDualsSingletonRow(int row, int col);
  void getDualsDoubletonEquation(int row, int col);
  void recordCounts(const string fileName);
  void trimA();

  void setBasisElement(change c);

  // test basis matrix singularity
  //
  // public:
  //	vector<int> nbffull;
  //	vector<int> bindfull;
  //	void cmpNBF(int row, int col);
  //	void setNBFfullproblem(vector<int>& nbfFull, vector<int>& bnFull);
  //	int testBasisMatrixSingularity();
  //

  string countsFile;
};

// comment out whole class and see what the issue is.

// Class for easy communication between Presolve and Highs. A single
// instance of PresolveInfo handles a single presolve execution on one
// LP.
class PresolveInfo {
 public:
  PresolveInfo() {}
  // option_presolve : 0 means don't presolve.
  PresolveInfo(int option_presolve, const HighsLp& lp) {
    if (option_presolve) {
      lp_ = &lp;
      presolve_.push_back(Presolve());
     }
   }

  HighsLp& getReducedProblem();
  HighsPresolveStatus presolve_status_;
  HighsPostsolveStatus postsolve_status_;
 public: 
  // Original problem is lp_.
  const  HighsLp * lp_;
  std::vector<Presolve> presolve_;
  HighsLp reduced_lp_;

  HighsSolution reduced_solution_;
  HighsSolution recovered_solution_;
};
#endif /* PRESOLVE_HPRESOLVE_H_ */
