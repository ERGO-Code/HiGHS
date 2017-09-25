#ifndef HPRESOLVE_H_
#define HPRESOLVE_H_

#include <string>
#include <vector>
#include <list>
#include <stack>
#include <utility>

#include "KktChStep.h"
#include "HTimerPre.h"
#include "HFactor.h"
#include <stdexcept>

#include "HPreData.h"

using namespace std;

class HPresolve : public HPreData {
public:
	HPresolve();
	string modelName;

	int iPrint;
	int iKKTcheck;
	int  presolve(int print);
	int  presolve();

	void postsolve();

 	double objShift;
 	void initializeVectors();
	void setProblemStatus(const int s);
 	void reportTimes();
 	
 	//new bounds on primal variables for implied free detection
	vector<double> implColLower;
    vector<double> implColUpper;
	vector<int> implColLowerRowIndex;
	vector<int> implColUpperRowIndex;

	vector<int> implRowDualLowerSingColRowIndex;
	vector<int> implRowDualUpperSingColRowIndex;

	//new bounds on row duals y_i
	vector<double> implRowDualLower;
    vector<double> implRowDualUpper;



    vector<double> implColDualLower;
    vector<double> implColDualUpper;
    vector<double> implRowValueLower;
    vector<double> implRowValueUpper;

    HTimerPre timer; //holds enum for main presolve rules

    enum stat {
    	Infeasible = 1,
    	Unbounded = 2,
    	Empty = 3
    };

private: 

    bool hasChange;
    int debug;

    int status=0;  //0 is unassigned, see enum stat

    list<int> singRow;    		//singleton rows 
	list<int> singCol;    		//singleton columns

	//original data
public:
	vector<double> colCostOriginal;

private:
	vector<double> rowLowerOriginal;
	vector<double> rowUpperOriginal;
	vector<double> colLowerOriginal;
	vector<double> colUpperOriginal;


	//functions
    void checkForChanges(int iteration);
    void setPrimalValue(int j, double value);
    void resizeProblem();
    void removeRowSingletons();
    int  getSingRowElementIndexInAR(int i);
    int  getSingColElementIndexInA(int j);
    void removeForcingConstraints(int mainIter);
    void removeRow(int i);
    void removeColumnSingletons();
    void addChange(int type, int row, int col);

    void removeDominatedColumns();
    bool removeIfImpliedFree(int col, int i, int k);
    pair<double, double> getNewBoundsDoubletonConstraint(int row, int col, int j, double aik, double aij);

	void removeIfFixed(int j);
	void removeEmptyRow(int i);
	void removeEmptyColumn(int j);
	void resizeImpliedBounds();
	void removeDoubletonEquations();


//    void findDuplicateRows();
//    void findDuplicateColumns();
//    void removeDuplicateRows(int i, int k, double v);
//    int makeCheckForDuplicateRows(int k, int i, vector<double>& coeff, vector<int>& colIndex, double v, int whichIsFirst);
//    void removeDuplicateColumns(int j,int k, double v);
//    bool checkDuplicateRows(int i, int k) ;
//	  bool checkDuplicateColumns(int i, int k) ;


	//old or test
	//void updateRemovedColRow(int dim);
	//void updateRowsByNZ();
	int  testSingRows(int i);
    int  testSingCols(int i);
    void testAnAR(int post);

    vector<int> countRemovedRows;
    vector<int> countRemovedCols;
    double tol;
    
    //postsolve
    bool noPostSolve;
    
    void fillStackRowBounds(int col);
	void setKKTcheckerData();

	void getBoundOnLByZj(int row, int j, double* lo, double* up, double colLow, double colUpp);
	double getRowDualPost(int row, int col);
	double getColumnDualPost(int col);
	string getDualsForcingRow( int row, vector<int>& fRjs);
	void getDualsSingletonRow( int row, int col );
	void getDualsDoubletonEquation(int row, int col);
	void recordCounts(string fileName);
	void trimA();

	void setBasisElement(change c);

	//test basis matrix singularity
//
//public:
//	vector<int> nbffull;
//	vector<int> bindfull;
//	void cmpNBF(int row, int col);
//	void setNBFfullproblem(vector<int>& nbfFull, vector<int>& bnFull);
//	int testBasisMatrixSingularity();
//

	string countsFile;
};

#endif /* HPRESOLVE_H_ */
