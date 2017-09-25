#include "HPresolve.h"
#include "HConst.h"

#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <queue>
#include <sstream>
#include <algorithm>

using namespace std;

int HPresolve::presolve(int print) {

	iPrint = print;
	iKKTcheck = 0;

	chk.print = 3; // 3 for experiments mode
	if (chk.print==3) {
		iPrint = 0;
		if (iKKTcheck) {
			iKKTcheck = 2;
			countsFile = "../experiments/t2";
		}
	}
	if (iPrint || iKKTcheck)
		debug = 1;
		
	//counter for the different types of reductions
	countRemovedCols.resize(HTICK_ITEMS_COUNT_PRE, 0);
	countRemovedRows.resize(HTICK_ITEMS_COUNT_PRE, 0);
	
	if (iPrint > 0) {
		cout<<"Presolve started ..."<<endl;
		cout<<"Original problem ... N="<<numCol<<"  M="<<numRow<<endl;
	}

	int i,j,k;

	initializeVectors();
	if (status) return status;

	int iter = 1;
	//print(0);

	timer.recordStart(HTICK_PRE_FIXED);
	for (int j=0;j<numCol;j++)
		if (flagCol.at(j)) {
			removeIfFixed(j);
			if (status) return status;
		}
	timer.recordFinish(HTICK_PRE_FIXED);


	while (hasChange) {

		hasChange = false;
		if (iPrint > 0)
			cout<<"PR: main loop "<<iter<<":"<<endl;
		//***************** main loop ******************

		removeRowSingletons();
		if (status) return status;
		removeForcingConstraints(iter);
		if (status) return status;

		removeRowSingletons();
		if (status) return status;
		removeDoubletonEquations();
		if (status) return status;

		removeRowSingletons();
		if (status) return status;
		removeColumnSingletons();
		if (status) return status;

		removeDominatedColumns();
		if (status) return status;


		//***************** main loop ******************
		iter++;
	}

	checkForChanges(iter);
	
	if (countsFile.length() > 0) {
		recordCounts(countsFile);
	}
	return status;
}

int HPresolve::presolve() {
	return presolve(0);
}

void HPresolve::removeDoubletonEquations() {
	if (flagCol.size() == numCol)
		flagCol.push_back(0);
	double col1, col2, x, y, b, low, upp;
	int iter = 0;

	for (int row=0;row<numRow;row++)
		if (flagRow[row])
			if (nzRow[row]==2 && abs(rowLower[row]-rowUpper[row])< tol    ) {
				//row is of form akx_x + aky_y = b, where k=row and y is present in fewer constraints

				double b = rowLower[row];
				col1 = col2 = -1;
				int kk = ARstart[row];
				while (kk<ARstart[row+1]) {
					if (flagCol[ARindex[kk]]) {
						if (col1==-1)
							col1 = ARindex[kk];
						else if (col2==-1)
							col2 = ARindex[kk];
						else {
							cout<<"ERROR: doubleton eq row"<< row << " has more than two variables. \n";
							col2 = -2;
							break;
						}
						kk++;
					}
					else
						kk++;
				}
				if (col2==-1)
					cout<<"ERROR: doubleton eq row"<< row << " has less than two variables. \n";
				if (col2 < 0)
					continue;

				if (nzCol[col1] <= nzCol[col2]) {
					y = col1;
					x = col2;
				}
				else {
					x = col1;
					y = col2;
				}

				double akx = getaij(row, x);
				double aky = getaij(row, y);


				/***********************************************
				//different cases:

				//singletons skip this
				for (int k = Astart[y]; k<Aend[y]; k++)
					if (flagRow[Aindex.at(k)]) {
						int i = Aindex.at(k);
						if (i==row)
							continue;
						double aiy = Avalue.at(k);
						if (isZeroA(i,x)) {
							//col2 = -1;   break;	 // no X
						}
						else {
							double xNew;
							int ind;
							for (ind = ARstart.at(i); ind<ARstart[i+1]; ind++)
								if (ARindex[ind] == x)
									break;

							xNew = ARvalue[ind] - (aiy*akx)/aky;
							//if (abs(xNew)>tol) { col2 = -1;		break;	}    // new x != 0
							//if (abs(xNew)<tol) { col2 = -1;		break;	}    // new x == 0 only nonsingleton y anyway, singletons separately
						}
					}

				//singleton case
				//if (nzCol[y] == 1)  col2 = -1;                 // nonsingleton y
				if (col2 < 0)
					continue;  */
				//***********************************************

				if (nzCol[y] == 1 && nzCol[x] == 1 )            //two singletons case handled elsewhere
					continue;

				postValue.push(akx);
				postValue.push(aky);
				postValue.push(b);

				//modify bounds on variable x (j), variable y (col,k) is substituted out
				//double aik = Avalue.at(k);
				//double aij = Avalue[kk];
				pair<double, double> p = getNewBoundsDoubletonConstraint(row, y, x, aky, akx);
				low = get<0>(p);
				upp = get<1>(p);

				//add old bounds of x to checker and for postsolve
				if (iKKTcheck == 1) {
					vector<pair<int, double>> bndsL, bndsU, costS;
					bndsL.push_back( make_pair( x, colLower[x]));
					bndsU.push_back( make_pair( x, colUpper[x]));
					costS.push_back( make_pair( x, colCost[x]));
					chk.cLowers.push(bndsL);
					chk.cUppers.push(bndsU);
					chk.costs.push(costS);
				}

				vector<double> bnds, bnds2, bnds3;
				bnds.push_back(colLower[y]);
				bnds.push_back(colUpper[y]);
				bnds.push_back(colCost[y]);
				oldBounds.push(make_pair( y, bnds));

				bnds2.push_back(colLower[x]);
				bnds2.push_back(colUpper[x]);
				bnds2.push_back(colCost[x]);
				oldBounds.push(make_pair( x, bnds2));


				if (low > colLower[x])
					colLower[x] = low;
				if (upp < colUpper[x])
					colUpper[x] = upp;

				//modify cost of xj
				colCost[x] = colCost[x] - colCost[y]*akx/aky;

				//for postsolve: need the new bounds too

				bnds3.push_back(colLower[x]);
				bnds3.push_back(colUpper[x]);
				bnds3.push_back(colCost[x]);
				oldBounds.push(make_pair( x, bnds3));

				addChange(17, row, y);

				//remove y (col) and the row
				if (iPrint > 0)
					//cout<<"PR: Doubleton equation removed. Row "<<row<<", column "<<y<<", column left is "<<x<<endl;
					cout<<"PR: Doubleton equation removed. Row "<<row<<", column "<<y<<", column left is "<<x<<"    nzy="<<nzCol[y]<<endl;
				flagRow[row] = 0;
				nzCol[x]--;

				countRemovedRows[HTICK_PRE_DOUBLETON_EQUATION]++;
				countRemovedCols[HTICK_PRE_DOUBLETON_EQUATION]++;

				//----------------------------
				flagCol[y] = 0;
				if (!hasChange)
					hasChange = true;


				vector<pair<int, double>> bndsL, bndsU;

				for (int k = Astart[y]; k<Aend[y]; k++)
					if (flagRow[Aindex.at(k)] && Aindex.at(k) != row) {
						int i = Aindex.at(k);
						double aiy = Avalue.at(k);

						//update row bounds
						if (iKKTcheck == 1) {
							bndsL.push_back( make_pair( i, rowLower.at(i)));
							bndsU.push_back( make_pair( i, rowUpper.at(i)));
							chk.rLowers.push(bndsL);
							chk.rUppers.push(bndsU);
							addChange(171, i, y);
						}

						if (rowLower.at(i) > -HSOL_CONST_INF)
							rowLower.at(i) -= b*aiy/aky;
						if (rowUpper.at(i) < HSOL_CONST_INF)
							rowUpper.at(i) -= b*aiy/aky;

						if (implRowValueLower.at(i) > -HSOL_CONST_INF)
							implRowValueLower.at(i) -= b*aiy/aky;
						if (implRowValueUpper.at(i) < HSOL_CONST_INF )
							implRowValueUpper.at(i) -= b*aiy/aky;


						//update matrix coefficients
						if (isZeroA(i,x)) {
							//case x is zero initially
							// row nonzero count doesn't change here
							//cout<<"case: x not present "<<i<<" "<<endl;

							//update AR
							int ind;
							for (ind = ARstart.at(i); ind<ARstart[i+1]; ind++)
								if (ARindex[ind] == y) {
									break;
								}
							postValue.push(ARvalue[ind]);
							postValue.push(y);
							addChange(173, i, x);

							ARindex[ind] = x;
							ARvalue[ind] = -aiy*akx/aky;

							//just row rep in checker
							if (iKKTcheck == 1) {
								bndsL.push_back( make_pair( i, rowLower.at(i)));
								bndsU.push_back( make_pair( i, rowUpper.at(i)));
								chk.ARvalue[ind] = ARvalue[ind];
								chk.ARindex[ind] = ARindex[ind];
							}

							//update A: append X column to end of array
							int st = Avalue.size();
							for (int ind = Astart[x]; ind < Aend[x]; ind++) {
								Avalue.push_back(Avalue[ind]);
								Aindex.push_back(Aindex[ind]);
							}
							Avalue.push_back(-aiy*akx/aky);
							Aindex.push_back(i);
							Astart[x] = st;
							Aend[x] = Avalue.size();

							nzCol[x]++;
							//nzRow does not change here.
							if (nzCol[x] == 2)
								singCol.remove(x);

						}
						else {
							int ind;

							//update nonzeros: for removal of
							nzRow.at(i)--;
							if (nzRow.at(i)==1)
								singRow.push_back(i);
							if (nzRow.at(i)==0) {
								singRow.remove(i);
								removeEmptyRow(i);
								countRemovedRows[HTICK_PRE_DOUBLETON_EQUATION]++;
							}

							double xNew;
							for (ind = ARstart.at(i); ind<ARstart[i+1]; ind++)
								if (ARindex[ind] == x)
									break;

							xNew = ARvalue[ind] - (aiy*akx)/aky;
							if (abs(xNew)>tol) {
								//case new x != 0
								//cout<<"case: x still there row "<<i<<" "<<endl;

								postValue.push(ARvalue[ind]);
								addChange(172, i, x);
								ARvalue[ind] = xNew;

								if (iKKTcheck == 1)
									chk.ARvalue[ind] = xNew;

								//update A:
								for (ind = Astart[x]; ind<Aend[x]; ind++)
									if (Aindex[ind] ==i) {
										break;
									}
								Avalue[ind] = xNew;

							}
							else if (xNew<tol) {
								//case new x == 0
								//cout<<"case: x also disappears from row "<<i<<" "<<endl;
								//update nz row
								nzRow.at(i)--;
								//update singleton row list
								if (nzRow.at(i)==1)
									singRow.push_back(i);
								if (nzRow.at(i)==0) {
									singRow.remove(i);
									removeEmptyRow(i);
									countRemovedRows[HTICK_PRE_DOUBLETON_EQUATION]++;
								}



								if (nzRow.at(i) > 0) {
									// AR update
									//set ARindex of element for x to numCol
									//flagCol[numCol] = false
									//mind when resizing: should be OK
									postValue.push(ARvalue[ind]);

									ARindex[ind] = numCol;
									if (iKKTcheck == 1) {
										chk.ARindex[ind] = ARindex[ind];
										chk.ARvalue[ind] = ARvalue[ind];
									}


									addChange(174, i, x);
								}

								if (nzCol[x] > 0) {
									// A update for case when x is zero: move x entry to end and set
									// Aend to be Aend - 1;
									int indi;
									for (indi=Astart[x];indi<Aend[x];indi++)
											if (Aindex[indi]==i)
												break;

									postValue.push(Avalue[indi]);

									//if indi is not Aend-1 swap elements indi and Aend-1
									if (indi != Aend[x]-1) {
										double tmp = Avalue[Aend[x]-1];
										int   tmpi = Aindex[Aend[x]-1];
										Avalue[Aend[x]-1] = Avalue[indi];
										Aindex[Aend[x]-1] = Aindex[indi];
										Avalue[indi] = tmp;
										Aindex[indi] = tmpi;

									}
									Aend[x]--;
									addChange(175, i, x);
								}

								//update nz col
								nzCol[x]--;
								//update singleton col list
								if (nzCol[x]==1)
									singCol.push_back(x);
								if (nzCol[x]==0) {
									nzRow.at(i)++;  //need this because below we decrease it by 1 too
									removeEmptyColumn(x);
								}
							}
						}
					}
				if (Avalue.size() > 40000000) {
					trimA();
				}

				iter++;
			}
}

void HPresolve::trimA() {
	int cntEl=0;
	for (int j=0;j<numCol;j++)
		if (flagCol.at(j))
			cntEl+=nzCol.at(j);

	vector<pair<int,size_t> > vp;
	vp.reserve(numCol);

	for (size_t i = 0 ; i != numCol ; i++) {
		vp.push_back(make_pair(Astart.at(i), i));
	}

	// Sorting will put lower values ahead of larger ones,
	// resolving ties using the original index
	sort(vp.begin(), vp.end());

//	for (size_t i = 0 ; i != vp.size() ; i++) {
//		cout << vp.at(i).first << " " << vp.at(i).second << endl;
//	}

	vector<int> Aendtmp;
	Aendtmp = Aend;

	int iPut=0;
	for (size_t i = 0 ; i != vp.size() ; i++) {
		int col = vp.at(i).second;
		if (flagCol[col]) {
			int k = vp.at(i).first; // = Astarttmp[col]
			Astart[col] = iPut;
			while(k < Aendtmp[col]) {
				if (flagRow[Aindex.at(k)]) {
					Avalue[iPut] = Avalue.at(k);
					Aindex[iPut] = Aindex.at(k);

					iPut++;
				}
				k++;
			}
			Aend[col] = iPut;
		}
	}
	Avalue.resize(iPut);
	Aindex.resize(iPut);

}


void HPresolve::resizeProblem() {

	int i, j, k;

	int nz = 0;
	int nR = 0;
	int nC = 0;

	//arrays to keep track of indices
	rIndex.assign(numRow, -1);
	cIndex.assign(numCol, -1);

	for (i=0;i<numRow;i++)
		if (flagRow.at(i)) {
			nz += nzRow.at(i);
			rIndex.at(i) = nR;
			nR++;
			}

	for (i=0;i<numCol;i++)
		if (flagCol.at(i)) {
			cIndex.at(i) = nC;
			nC++;
		}

	//counts
	numRowOriginal = numRow;
	numColOriginal = numCol;
	numRow = nR;
	numCol = nC;
	numTot = nR + nC;

	if (1) {
    //if (iPrint == -1) {
    	cout<<"Presolve m="<<setw(2)<<numRow<<"(-"<<setw(2)<< numRowOriginal - numRow <<") ";
    	cout<<"n="<<setw(2)<< numCol <<"(-"<<setw(2)<< numColOriginal - numCol <<") ";
    	cout<<"nz="<<setw(4)<< nz << "(-"<<setw(4)<< ARindex.size() - nz  <<") ";
		cout<<endl;
    }

	if (nR + nC == 0) {
    	status = Empty;
		return;
	}

	//matrix
    vector<int> iwork(numCol, 0);
    Astart.assign(numCol + 1, 0);
    Aend.assign(numCol + 1, 0);
    Aindex.resize(nz);
    Avalue.resize(nz);


    for (i = 0;i<numRowOriginal; i++)
    	if (flagRow.at(i))
    	    for (int k = ARstart.at(i); k < ARstart[i+1]; k++) {
    	    	j = ARindex.at(k);
    	    	if (flagCol.at(j))
        			iwork[cIndex.at(j)]++;
        		}
    for (i = 1; i <= numCol; i++)
        Astart.at(i) = Astart[i - 1] + iwork[i - 1];
    for (i = 0; i < numCol; i++)
        iwork.at(i) = Aend.at(i) = Astart.at(i);
    for (i = 0; i < numRowOriginal; i++) {
    	if (flagRow.at(i)) {
			int iRow = rIndex.at(i);
		    for (k = ARstart.at(i); k < ARstart[i + 1]; k++) {
		        j = ARindex.at(k);
		        if (flagCol.at(j)) {
		        	int iCol = cIndex.at(j);
				    int iPut = iwork[iCol]++;
				    Aindex[iPut] = iRow;
				    Avalue[iPut] = ARvalue.at(k);
				}
		    }
		}
    }
     //For KKT checker: pass vectors before you trim them
    if (iKKTcheck == 1) {
		chk.setFlags(flagRow, flagCol);
		chk.setBoundsCostRHS(colUpper, colLower, colCost, rowLower, rowUpper);
	}

    //also call before trimming
    resizeImpliedBounds();

    //cost, bounds
    colCostAtEl = colCost;
    vector<double> tempCost = colCost;
    vector<double> temp = colLower;
    vector<double> teup = colUpper;

    colCost.resize(numCol);
    colLower.resize(numCol);
    colUpper.resize(numCol);

    k=0;
    for (i=0;i<numColOriginal;i++)
    	if (flagCol.at(i)) {
    		colCost.at(k)  = tempCost.at(i);
    		colLower.at(k) = temp.at(i);
    		colUpper.at(k) = teup.at(i);
    		k++;
	    }

    //RHS and bounds
    rowLowerAtEl = rowLower;
    rowUpperAtEl = rowUpper;
    temp = rowLower;
    teup = rowUpper;
    rowLower.resize(numRow);
    rowUpper.resize(numRow);
    k=0;
    for (i=0;i<numRowOriginal;i++)
    	if (flagRow.at(i)) {
    		rowLower.at(k) = temp.at(i);
    		rowUpper.at(k) = teup.at(i);
    		k++;
	    }

    if (chk.print == 3) {
		ofstream myfile;
  		myfile.open ("../experiments/out", ios::app );
		myfile << " eliminated rows "<<(numRowOriginal - numRow)<<" cols "<<(numColOriginal - numCol);
		myfile.close();

		myfile.open ("../experiments/t3", ios::app );
		myfile<<(numRowOriginal )<<"  &  "<<(numColOriginal ) << "  & ";
		myfile<<(numRowOriginal - numRow)<<"  &  "<<(numColOriginal - numCol) << "  & "<<endl;


		myfile.close();

	}
}


void HPresolve::initializeVectors() {

	//copy original bounds
	colCostOriginal = colCost;
	rowUpperOriginal = rowUpper;
	rowLowerOriginal = rowLower;
	colUpperOriginal = colUpper;
	colLowerOriginal = colLower;

	makeARCopy();

	valueRowDual.resize(numRow);
	valuePrimal.resize(numCol);
	valueColDual.resize(numCol);

	flagCol.assign(numCol, 1);
	flagRow.assign(numRow, 1);

	if (iKKTcheck)
		setKKTcheckerData();

	nzCol.assign(numCol,0);
	nzRow.assign(numRow,0);

	for (int i = 0; i < numRow; i++) {
		nzRow.at(i) = ARstart[i+1]-ARstart.at(i);
		if (nzRow.at(i) == 1)
			singRow.push_back(i);
		if (nzRow.at(i) == 0) {
			timer.recordStart(HTICK_PRE_EMPTY_ROW);
			removeEmptyRow(i);
			countRemovedRows[HTICK_PRE_EMPTY_ROW]++;
			timer.recordFinish(HTICK_PRE_EMPTY_ROW);
		}
	}

	Aend.resize(numCol+1);
	for (int i = 0; i < numCol; i++) {
		Aend.at(i)  = Astart[i+1];
		nzCol.at(i) = Aend.at(i)-Astart.at(i);
		if (nzCol.at(i) == 1)
			singCol.push_back(i);
	}
	objShift = 0;

	implColUpper = colUpper;	//working copies of primal variable bounds
	implColLower = colLower;
	implColLowerRowIndex.assign(numCol,-1);
	implColUpperRowIndex.assign(numCol,-1);

	implRowDualLowerSingColRowIndex.assign(numRow,-1);
	implRowDualUpperSingColRowIndex.assign(numRow,-1);
	implRowDualLower.assign(numRow, -HSOL_CONST_INF);
	implRowDualUpper.assign(numRow, HSOL_CONST_INF);

	implColDualLower.assign(numCol, -HSOL_CONST_INF);
    implColDualUpper.assign(numCol, HSOL_CONST_INF);
    implRowValueLower = rowLower;
    implRowValueUpper = rowUpper;

    for (int i = 0; i < numRow; i++) {
		if (rowLower.at(i) == -HSOL_CONST_INF)
			implRowDualUpper.at(i) = 0;
		if (rowUpper.at(i) == HSOL_CONST_INF)
			implRowDualLower.at(i) = 0;
	}

    for (int i = 0; i < numCol; i++) {
    	if (colLower.at(i) == -HSOL_CONST_INF)
    		implColDualUpper.at(i) = 0;
    	if (colUpper.at(i) == HSOL_CONST_INF)
    		implColDualLower.at(i) = 0;
    }

	colCostAtEl = colCost;
	rowLowerAtEl = rowLower;
    rowUpperAtEl = rowUpper;
}

HPresolve::HPresolve() {
	tol = 0.0000001;
	noPostSolve = false;
	objShift = 0;
	hasChange = true;
	iKKTcheck = 0;
	iPrint = 0;
	debug = 0;
	countsFile = "";
}

void HPresolve::removeIfFixed(int j) {

		if (colLower.at(j) == colUpper.at(j)) {

			setPrimalValue(j, colUpper.at(j));
			addChange(7,0,j);
			if (iPrint > 0)
				cout<<"PR: Fixed variable "<<j<<" = "<<colUpper.at(j)<<". Column eliminated."<< endl;
			countRemovedCols[HTICK_PRE_FIXED]++;

			for (int k = Astart.at(j); k < Aend.at(j); k++) {
				if (flagRow[Aindex.at(k)]) {
					int i = Aindex.at(k);

					if (nzRow.at(i) == 0) {
						removeEmptyRow(i);
						countRemovedRows[HTICK_PRE_FIXED]++;
					}
				}
			}
		}
}

void HPresolve::removeEmptyRow(int i) {
	if (rowLower.at(i) <= tol && rowUpper.at(i) >= -tol) {
		if (iPrint > 0)
			cout<<"PR: Empty row "<<i<<" removed. "<<endl;
		flagRow.at(i) = 0;
		valueRowDual.at(i) = 0;
		addChange(0, i, 0);
	}
	else {
		if (iPrint > 0)
			cout<<"PR: Problem infeasible."<<endl;
		status = Infeasible;
		return;
	}
}

void HPresolve::removeEmptyColumn(int j) {
	flagCol.at(j) = 0;
	singCol.remove(j);
	double value;
	if ((colCost.at(j) < 0 && colUpper.at(j) == HSOL_CONST_INF) || (colCost.at(j) > 0 && colLower.at(j) == -HSOL_CONST_INF) ) {
		if (iPrint > 0)
			cout<<"PR: Problem unbounded."<<endl;
		status = Unbounded;
		return;
	}

	if (colCost.at(j) > 0)
		value = colLower.at(j);
	else if (colCost.at(j) < 0)
		value = colUpper.at(j);
	else if (colUpper.at(j) >= 0 && colLower.at(j) <=0)
		value = 0;
	else if (colUpper.at(j) < 0)
		value = colUpper.at(j);
	else
		value = colLower.at(j);

	valuePrimal.at(j) = value;
	setPrimalValue(j, value);
	valueColDual.at(j) = colCost.at(j);

	addChange(6, 0, j);

	if (iPrint > 0)
		cout<<"PR: Column: "<<j<<" eliminated: all nonzero rows have been removed. Cost = "<< colCost.at(j) <<", value = "<<value<<endl;
	countRemovedCols[HTICK_PRE_EMPTY_COL]++;
}


void HPresolve::removeDominatedColumns() {
	int col, i, k;
	
	//for each row calc yihat and yibar and store in implRowDualLower and implRowDualUpper
	for (list<int>::iterator it = singCol.begin(); it != singCol.end(); ++it) 
		if (flagCol[*it]) {
			col = *it;
			k = getSingColElementIndexInA(col);
			i = Aindex.at(k);
			if (!flagRow.at(i)) {
				cout<<"ERROR: column singleton "<<col<<" is in row "<<i<<" which is already mapped off\n";
				exit(18);
			}

			if (colLower[col] == -HSOL_CONST_INF || colUpper[col] == HSOL_CONST_INF) {

				if (colLower[col] > -HSOL_CONST_INF && colUpper[col] == HSOL_CONST_INF)  {
					if (Avalue.at(k)>0)
						if ((colCost[col]/Avalue.at(k)) < implRowDualUpper.at(i))
							implRowDualUpper.at(i) = colCost[col]/Avalue.at(k);
					if (Avalue.at(k)<0)
						if ((colCost[col]/Avalue.at(k)) > implRowDualLower.at(i))
							implRowDualLower.at(i) = colCost[col]/Avalue.at(k);
				}
				else if (colLower[col] == -HSOL_CONST_INF && colUpper[col] < HSOL_CONST_INF) {
					if (Avalue.at(k)>0)
						if ((colCost[col]/Avalue.at(k)) > implRowDualLower.at(i))
							implRowDualUpper.at(i) = -colCost[col]/Avalue.at(k);
					if (Avalue.at(k)<0)
						if ((colCost[col]/Avalue.at(k)) < implRowDualUpper.at(i))
							implRowDualUpper.at(i) = colCost[col]/Avalue.at(k);
				}
				else if (colLower[col] == -HSOL_CONST_INF && colUpper[col] == HSOL_CONST_INF) {
					//all should be removed earlier but use them
					if ((colCost[col]/Avalue.at(k)) > implRowDualLower.at(i))
							implRowDualLower.at(i) = colCost[col]/Avalue.at(k);
					if ((colCost[col]/Avalue.at(k)) < implRowDualUpper.at(i))
							implRowDualUpper.at(i) = colCost[col]/Avalue.at(k);
				}
				
				if (implRowDualLower.at(i) > implRowDualUpper.at(i)) {
					cout<<"Error: inconstistent bounds for Lagrange multiplier for row "<< i<< " detected after column singleton "<<col<<". In presolve::dominatedColumns"<<endl;
					exit(0); 
				}	
			}		
		}
	
	
	//for each column j calculate e and d and check:
	double e,d;
	for(int j=0;j<numCol;j++) 
		if (flagCol.at(j)) {
			timer.recordStart(HTICK_PRE_DOMINATED_COLS);
			e=0;
			d=0;

			for (k=Astart.at(j); k<Aend.at(j); k++) {
				i = Aindex.at(k);
				if (flagRow.at(i))
					if (Avalue.at(k) < 0) {
						if (implRowDualUpper.at(i) < HSOL_CONST_INF)
							e+= Avalue.at(k)*implRowDualUpper.at(i);
						else {
							e = -HSOL_CONST_INF;
							break;
						}
						
					}
					else {
						if (implRowDualLower.at(i) > -HSOL_CONST_INF)
							e+= Avalue.at(k)*implRowDualLower.at(i);
						else {
							e = -HSOL_CONST_INF;
							break;
						}
					}	
			}
			for (k=Astart.at(j); k<Aend.at(j); k++) {
				i = Aindex.at(k);
				if (flagRow.at(i))
					if (Avalue.at(k) < 0) {
						if (implRowDualLower.at(i) > -HSOL_CONST_INF)
							d+= Avalue.at(k)*implRowDualLower.at(i);
						else {
							d = HSOL_CONST_INF;
							break;
						}
						
					}
					else {
						if (implRowDualUpper.at(i) < HSOL_CONST_INF)
							d+= Avalue.at(k)*implRowDualUpper.at(i);
						else {
							d = HSOL_CONST_INF;
							break;
						}
					}	
			}
	
			if (e>d) {
				cout<<"Error: inconstistent bounds for Lagrange multipliers for column "<<j<<": e>d. In presolve::dominatedColumns"<<endl;
				exit(-2); 
			}	

			//check if it is dominated 
			if (colCost.at(j) - d > tol) {
				if (colLower.at(j) == -HSOL_CONST_INF) {
					if (iPrint > 0)
						cout<<"PR: Problem unbounded."<<endl;
					status = Unbounded;
					return;
				}
				setPrimalValue(j, colLower.at(j));
				addChange(9, 0, j);
				if (iPrint > 0)
					cout<<"PR: Dominated column "<<j<<" removed. Value := "<<valuePrimal.at(j)<<endl;
				timer.recordFinish(HTICK_PRE_DOMINATED_COLS);
				countRemovedCols[HTICK_PRE_DOMINATED_COLS]++;
			}
			else if (colCost.at(j) - e < -tol) {
				if (colUpper.at(j) == HSOL_CONST_INF) {
					if (iPrint > 0)
						cout<<"PR: Problem unbounded."<<endl;
					status = Unbounded;
					return;
				}
				setPrimalValue(j, colUpper.at(j));
				addChange(9, 0, j);
				if (iPrint > 0)
					cout<<"PR: Dominated column "<<j<<" removed. Value := "<<valuePrimal.at(j)<<endl;
				timer.recordFinish(HTICK_PRE_DOMINATED_COLS);
				countRemovedCols[HTICK_PRE_DOMINATED_COLS]++;
			}
			else {
				if (implColDualLower.at(j) < (colCost.at(j) - d ))
					implColDualLower.at(j) = colCost.at(j) - d;
				if (implColDualUpper.at(j) > (colCost.at(j) - e ))
					implColDualUpper.at(j) = colCost.at(j) - e;
				if (implColDualLower.at(j) > implColDualUpper.at(j) )
					cout<<"INCONSISTENT\n";

				timer.recordFinish(HTICK_PRE_DOMINATED_COLS);


				//check if it is weakly dominated: Excluding singletons!
				if (nzCol.at(j)>1) {
					if (abs(colCost.at(j) - d) < tol && colLower.at(j) > -HSOL_CONST_INF) {
						timer.recordStart(HTICK_PRE_WEAKLY_DOMINATED_COLS);
						setPrimalValue(j, colLower.at(j));
						addChange(10, 0, j);
						if (iPrint > 0)
							cout<<"PR: Weakly Dominated column "<<j<<" removed. Value := "<<valuePrimal.at(j)<<endl;
						countRemovedCols[HTICK_PRE_WEAKLY_DOMINATED_COLS]++;
						timer.recordFinish(HTICK_PRE_WEAKLY_DOMINATED_COLS);
					}
					else if (abs(colCost.at(j) - e) < tol && colUpper.at(j) < HSOL_CONST_INF) {
						timer.recordStart(HTICK_PRE_WEAKLY_DOMINATED_COLS);
						setPrimalValue(j, colUpper.at(j));
						addChange(10, 0, j);
						if (iPrint > 0)
							cout<<"PR: Weakly Dominated column "<<j<<" removed. Value := "<<valuePrimal.at(j)<<endl;
						countRemovedCols[HTICK_PRE_WEAKLY_DOMINATED_COLS]++;
						timer.recordFinish(HTICK_PRE_WEAKLY_DOMINATED_COLS);
					}
					else {

						timer.recordStart(HTICK_PRE_DOMINATED_COL_BOUNDS);
						double bnd;
						//calculate new bounds
						if (colLower.at(j) > -HSOL_CONST_INF || colUpper.at(j) == HSOL_CONST_INF)
							for (int kk=Astart.at(j); kk<Aend.at(j);kk++)
								if (flagRow[Aindex[kk]] && d < HSOL_CONST_INF ) {
									i = Aindex[kk];
									if (Avalue[kk]>0 && implRowDualLower.at(i) > -HSOL_CONST_INF) {
										bnd = -(colCost.at(j) + d)/Avalue[kk] + implRowDualLower.at(i);
										if (bnd < implRowDualUpper.at(i) && !(bnd < implRowDualLower.at(i)))
											implRowDualUpper.at(i) = bnd;

									}
									else if (Avalue[kk]<0 && implRowDualUpper.at(i) < HSOL_CONST_INF) {
											bnd = -(colCost.at(j) + d)/Avalue[kk] + implRowDualUpper.at(i);
										if (bnd > implRowDualLower.at(i) && !(bnd > implRowDualUpper.at(i)))
											implRowDualLower.at(i) = bnd;

									}
								}

						if (colLower.at(j) == -HSOL_CONST_INF || colUpper.at(j) < HSOL_CONST_INF)
							for (int kk=Astart.at(j); kk<Aend.at(j);kk++)
								if (flagRow[Aindex[kk]] && e > -HSOL_CONST_INF ) {
									i = Aindex[kk];
									if (Avalue[kk]>0 && implRowDualUpper.at(i) < HSOL_CONST_INF) {
										bnd = -(colCost.at(j) + e)/Avalue[kk] + implRowDualUpper.at(i);
										if (bnd > implRowDualLower.at(i) && !(bnd > implRowDualUpper.at(i)))
											implRowDualLower.at(i) = bnd;

									}
									else  if (Avalue[kk]<0  && implRowDualLower.at(i) > -HSOL_CONST_INF)     {
										bnd = -(colCost.at(j) + e)/Avalue[kk] + implRowDualLower.at(i);
										if (bnd < implRowDualUpper.at(i) && !(bnd < implRowDualLower.at(i)))
											implRowDualUpper.at(i) = bnd;

									}
								}
						timer.recordFinish(HTICK_PRE_DOMINATED_COL_BOUNDS);
					}
				}
			}
		}
}

void HPresolve::setProblemStatus(int s) {
	if (s==Infeasible)
		cout<<"NOT-OPT status = 1, returned from solver after presolve: Problem infeasible.\n";
	else if (s==Unbounded)
		cout<<"NOT-OPT status = 2, returned from solver after presolve: Problem unbounded.\n";
	else if (s==0)
		return;
	else
		cout<<"unknown problem status returned from solver after presolve: "<<s<<endl;
	status = s;

}

void HPresolve::setKKTcheckerData(){
	//after initializing equations.
	chk.setMatrixAR(numCol, numRow, ARstart, ARindex, ARvalue);
	chk.setFlags(flagRow, flagCol); 
	chk.setBoundsCostRHS(colUpper, colLower, colCost, rowLower, rowUpper);
}

pair<double, double> HPresolve::getNewBoundsDoubletonConstraint(int row, int col, int j, double aik, double aij){
	int i=row;

	double upp = HSOL_CONST_INF;
	double low = -HSOL_CONST_INF;

	if ( aij > 0 && aik > 0 ) {
		if (colLower[col] > -HSOL_CONST_INF)
			upp = (rowUpper.at(i) - aik*colLower[col]) / aij;
		if (colUpper[col] < HSOL_CONST_INF)
			low = (rowLower.at(i) - aik*colUpper[col]) / aij;
	}
	else if (aij > 0 && aik < 0)  {
		if (colLower[col] > -HSOL_CONST_INF)
			low = (rowLower.at(i) - aik*colLower[col]) / aij;
		if (colUpper[col] < HSOL_CONST_INF)
			upp = (rowUpper.at(i) - aik*colUpper[col]) / aij;
	}
	else if ( aij < 0 && aik > 0 ) {
		if (colLower[col] > -HSOL_CONST_INF)
			low = (rowUpper.at(i) - aik*colLower[col]) / aij;
		if (colUpper[col] < HSOL_CONST_INF)
			upp = (rowLower.at(i) - aik*colUpper[col]) / aij;
	}
	else {
		if (colLower[col] > -HSOL_CONST_INF)
			upp = (rowLower.at(i) - aik*colLower[col]) / aij;
		if (colUpper[col] < HSOL_CONST_INF)
			low = (rowUpper.at(i) - aik*colUpper[col]) / aij;
	}

	return make_pair(low, upp);
}

void HPresolve::removeColumnSingletons()  {
	int i,j,k, col;
	list<int>::iterator it = singCol.begin();
	int iter=1;
	while (it != singCol.end()) {
		if (flagCol[*it]) {
			col = *it;
			k = getSingColElementIndexInA(col);
			i = Aindex.at(k);

			if (debug) {
				if (i != testSingCols(col)) {
				cout<<"ERROR: column "<<col<<" not singleton\n";
				it++;
				continue;}
			}

			if (!flagRow.at(i)) {
				cout<<"ERROR: column singleton "<<col<<" is in row "<<i<<" which is already mapped off\n";
				exit(18);
			}

			//free
			if (colLower[col] == -HSOL_CONST_INF && colUpper[col] == HSOL_CONST_INF) {
				timer.recordStart(HTICK_PRE_FREE_SING_COL);
				if (iPrint > 0)
					cout<<"PR: Free column singleton "<<col<<" removed. Row "<<i<<" removed."<<endl;
					
				//modify costs
				vector<pair<int, double>> newCosts;
				for (int kk=ARstart.at(i); kk<ARstart[i+1]; kk++) {
					j = ARindex[kk];
					if (flagCol.at(j) && j!=col) {
						newCosts.push_back(make_pair(j, colCost.at(j)));
						colCost.at(j) = colCost.at(j) -  colCost[col]*ARvalue[kk]/Avalue.at(k);
					}
				}		
				if (iKKTcheck == 1)
					chk.costs.push(newCosts);
				
				flagCol[col] = 0;
				postValue.push(colCost[col]);
				fillStackRowBounds(i);

				valueColDual[col] = 0;
				valueRowDual.at(i) = -colCost[col]/Avalue.at(k);
				addChange(4, i, col);
				removeRow(i);
				it = singCol.erase(it);
				countRemovedCols[HTICK_PRE_FREE_SING_COL]++;
				countRemovedRows[HTICK_PRE_FREE_SING_COL]++;
				timer.recordFinish(HTICK_PRE_FREE_SING_COL);
				continue;
		}
			//singleton column in a doubleton equation 
			else if (nzRow.at(i)==2) {
				//count
				int kk = ARstart.at(i);
				while (kk<ARstart[i+1]) {
					j = ARindex[kk];
					if (flagCol.at(j) && j!=col)
						break;
					else
						kk++;
				}
				if (kk==ARstart[i+1]) 
					cout<<"ERROR: nzRow["<< i<< "]=2, but no second variable in row. \n";

				//only inequality case and case two singletons here,
				//others in doubleton equation
				if (abs(rowLower.at(i)-rowUpper.at(i)) < tol) {
					if (nzCol.at(j)>1) {
						it++;
						continue;
					}
				}

				timer.recordStart(HTICK_PRE_SING_COL_DOUBLETON_INEQ);
				// additional check if it is indeed implied free: need
				// low and upp to be tighter than original bounds for variable col
				// so it is indeed implied free and we can remove it
				pair<double, double> p = getNewBoundsDoubletonConstraint(i, j, col, ARvalue[kk], Avalue.at(k));
				double low, upp;
				low = get<0>(p);
				upp = get<1>(p);
				if (!(colLower[col] <= low && colUpper[col] >= upp)) {
					it++;
					timer.recordFinish(HTICK_PRE_SING_COL_DOUBLETON_INEQ);
					continue;
				}

				postValue.push(ARvalue[kk]);
				postValue.push(Avalue.at(k));

				//modify bounds on variable j, variable col (k) is substituted out
				//double aik = Avalue.at(k);
				//double aij = Avalue[kk];
				p = getNewBoundsDoubletonConstraint(i, col, j, Avalue.at(k), ARvalue[kk]);
				low = get<0>(p);
				upp = get<1>(p);
				
				//add old bounds of xj to checker and for postsolve
				if (iKKTcheck == 1) {
					vector<pair<int, double>> bndsL, bndsU, costS;
					bndsL.push_back( make_pair( j, colLower.at(j)));
					bndsU.push_back( make_pair( j, colUpper.at(j)));
					costS.push_back( make_pair( j, colCost.at(j)));
					chk.cLowers.push(bndsL);
					chk.cUppers.push(bndsU);
					chk.costs.push(costS);
				}

				vector<double> bnds;
				bnds.push_back(colLower[col]);
				bnds.push_back(colUpper[col]);
				bnds.push_back(colCost[col]);
				oldBounds.push(make_pair( col, bnds));
				bnds.clear();
				bnds.push_back(colLower.at(j));
				bnds.push_back(colUpper.at(j));
				bnds.push_back(colCost.at(j));
				oldBounds.push(make_pair( j, bnds));
				
				if (low > colLower.at(j))
					colLower.at(j) = low;
				if (upp < colUpper.at(j))
					colUpper.at(j) = upp;
				
				//modify cost of xj
				colCost.at(j) = colCost.at(j) - colCost[col]*ARvalue[kk]/Avalue.at(k);

				//for postsolve: need the new bounds too
				//oldBounds.push_back(colLower.at(j)); oldBounds.push_back(colUpper.at(j));
				bnds.clear();
				bnds.push_back(colLower.at(j));
				bnds.push_back(colUpper.at(j));
				bnds.push_back(colCost.at(j));
				oldBounds.push(make_pair( j, bnds));
				
						//remove col as free column singleton
				if (iPrint > 0)
					cout<<"PR: Column singleton "<<col<<" in a doubleton inequality constraint removed. Row "<<i<<" removed. variable left is "<<j<<endl;
				flagCol[col] = 0;
				fillStackRowBounds(i);
				countRemovedCols[HTICK_PRE_SING_COL_DOUBLETON_INEQ]++;
				countRemovedRows[HTICK_PRE_SING_COL_DOUBLETON_INEQ]++;

				valueColDual[col] = 0;
				valueRowDual.at(i) = -colCost[col]/Avalue.at(k); //may be changed later, depending on bounds.

				addChange(5, i, col);
				
				if (nzCol.at(j) > 1)
					removeRow(i);
				else if (nzCol.at(j)==1) {
					// case two singleton columns
					// when we get here bounds on xj are updated so we can choose low/upper one
					// depending on the cost of xj
					flagRow.at(i) = 0;
					double value;
					if (colCost.at(j) > 0) {
						if (colLower.at(j) == -HSOL_CONST_INF) {
							if (iPrint > 0)
								cout<<"PR: Problem unbounded."<<endl;
							status = Unbounded;
							return;
						}
						value = colLower.at(j);
					}
					else if (colCost.at(j) < 0) {
						if (colUpper.at(j) == HSOL_CONST_INF) {
							if (iPrint > 0)
								cout<<"PR: Problem unbounded."<<endl;
							status = Unbounded;
							return;
						}
						value = colUpper.at(j);
					}
					else { //(colCost.at(j) == 0)
						if (colUpper.at(j) >= 0 && colLower.at(j) <= 0)
							value = 0;
						else if ( abs(colUpper.at(j)) < abs(colLower.at(j)) )
							value = colUpper.at(j);
						else
							value = colLower.at(j);
					}
					setPrimalValue(j, value);
					addChange(19, 0, j);
					if (iPrint > 0)
						cout<<"PR: Second singleton column "<<j<<" in doubleton row "<<i<< " removed.\n";
					countRemovedCols[HTICK_PRE_SING_COL_DOUBLETON_INEQ]++;
					singCol.remove(j);

				}
				
				it = singCol.erase(it); 
				iter++;
				timer.recordFinish(HTICK_PRE_SING_COL_DOUBLETON_INEQ);
				continue;

			}
			//implied free
			else{
				timer.recordStart(HTICK_PRE_IMPLIED_FREE_SING_COL);
				bool res = removeIfImpliedFree(col, i, k);
				if (res) {
					it = singCol.erase(it);
					iter++;
				}
				timer.recordFinish(HTICK_PRE_IMPLIED_FREE_SING_COL);
			}
		
			it++;
		}
		else 
			it = singCol.erase(it);

	}
}

bool HPresolve::removeIfImpliedFree(int col, int i, int k) {
	//first find which bound is active for row i
	//A'y + c = z so yi = -ci/aij
	double aij = getaij(i,col);
	if (aij != Avalue.at(k))
		cout<<"ERROR during implied free";
	double yi = -colCost[col]/aij;
	double low, upp;

	if (yi > 0) {
		if (rowUpper.at(i) == HSOL_CONST_INF)
			return false;
		low = rowUpper.at(i);
		upp = rowUpper.at(i);
	}
	else if (yi < 0) {
		if (rowLower.at(i) == -HSOL_CONST_INF)
			return false;
		low = rowLower.at(i);
		upp = rowLower.at(i);
	}
	else  {
		low = rowLower.at(i);
		upp = rowUpper.at(i);
	}

	//use implied bounds with original bounds
	int kk = ARstart.at(i);
	int j;
	double l,u;
	// if at any stage low becomes  or upp becomes inf break loop  
	// can't use bounds for variables generated by the same row. 
	//low
	for (int kk = ARstart.at(i); kk<ARstart[i+1]; kk++) {
		j = ARindex[kk];
		if (flagCol.at(j) && j!=col) {
			// check if new bounds are precisely implied bounds from same row
			if (i != implColLowerRowIndex.at(j))
				l = max(colLower.at(j), implColLower.at(j));
			else 
				l = colLower.at(j);
			if (i != implColUpperRowIndex.at(j))
				u = min(colUpper.at(j), implColUpper.at(j));
			else 
				u = colUpper.at(j);
		
			if ((Avalue.at(k)<0 && ARvalue[kk]>0) || (Avalue.at(k)>0 && ARvalue[kk]<0))
				if (l==-HSOL_CONST_INF) {
					low = -HSOL_CONST_INF;
					break; }
				else 
					low -= ARvalue[kk]*l;	
			else 
				if (u==HSOL_CONST_INF) {
					low = -HSOL_CONST_INF;
					break; }
				else 
					low -= ARvalue[kk]*u;
		}
	}
	//upp
	for (int kk = ARstart.at(i); kk<ARstart[i+1]; kk++) {
		j = ARindex[kk];
		if (flagCol.at(j) && j!=col) {
			// check if new bounds are precisely implied bounds from same row
			if (i != implColLowerRowIndex.at(j))
				l = max(colLower.at(j), implColLower.at(j));
			else 
				l = colLower.at(j);
			if (i != implColUpperRowIndex.at(j))
				u = min(colUpper.at(j), implColUpper.at(j));
			else 
				u = colUpper.at(j);
			// if at any stage low becomes  or upp becomes inf it's not implied free 
			//low:: 
			if ((Avalue.at(k)<0 && ARvalue[kk]>0) || (Avalue.at(k)>0 && ARvalue[kk]<0))
				if (u==HSOL_CONST_INF) {
					upp = HSOL_CONST_INF;
					break; }
				else 
					upp -= ARvalue[kk]*u;	
			else 
				if (l==-HSOL_CONST_INF) {
					upp = HSOL_CONST_INF;
					break; }
				else 
					upp -= ARvalue[kk]*l;	
		}
	}

	if (low>-HSOL_CONST_INF)
		low = low/Avalue.at(k);
	if (upp<HSOL_CONST_INF)				
		upp = upp/Avalue.at(k);

	if (colLower[col] <= low && low <= upp && upp <= colUpper[col]) {
		if (iPrint > 0)
			cout<<"PR: Implied free column singleton "<<col<<" removed.  Row "<<i<<" removed."<<endl;
		
		countRemovedCols[HTICK_PRE_IMPLIED_FREE_SING_COL]++;
		countRemovedRows[HTICK_PRE_IMPLIED_FREE_SING_COL]++;

		//modify costs
		vector<pair<int, double>> newCosts;
		for (int kk=ARstart.at(i); kk<ARstart[i+1]; kk++) {
			j = ARindex[kk];
			if (flagCol.at(j) && j!=col) {
				newCosts.push_back(make_pair(j, colCost.at(j)));
				colCost.at(j) = colCost.at(j) -  colCost[col]*ARvalue[kk]/Avalue.at(k);
			}
		}
		if (iKKTcheck == 1)
			chk.costs.push(newCosts);

		flagCol[col] = 0;
		postValue.push(colCost[col]);
		fillStackRowBounds(i);

		valueColDual[col] = 0;
		valueRowDual.at(i) = -colCost[col]/Avalue.at(k);
		addChange(8, i, col);
		removeRow(i);
		return true;

	}

	//implied bounds
	else if (colLower[col] <= low && low <= upp) {
		if (implColLower[col] < low) {
			implColLower[col] = low;
			implColUpperRowIndex[col] = i;
		}
		implColDualUpper.at(i) = 0;
	}
	else if (low <= upp && upp <= colUpper[col]) {
		if (implColUpper[col] > upp) {
			implColUpper[col] = upp;
			implColUpperRowIndex[col] = i;
		}
		implColDualLower[col] = 0;
	}


	return false;
}
	
	


//used to remove column too, now possible to just modify bounds
void HPresolve::removeRow(int i) {
	hasChange = true;
	flagRow.at(i) = 0;
	for(int k=ARstart.at(i);k<ARstart[i+1];k++) {
		int j = ARindex.at(k);
		if (flagCol.at(j)) {
			nzCol.at(j)--;
			//if now singleton add to list
			if (nzCol.at(j)==1) {
				if (debug) {
					int singletonColIndex = testSingCols(j);
					if (singletonColIndex >= 0)
							singCol.push_back(j);
					else
						cout<<"Warning: Column "<<j<<" with 1 nz but not in singCol or? Row removing of "<<i<<". Ignored.\n";
				}
				else
					singCol.push_back(j);
			}
			//if it was a singleton column remove from list and problem
			if (nzCol.at(j)==0) {
				removeEmptyColumn(j);
			}	
		}
	}
}


void HPresolve::fillStackRowBounds(int row) {
	postValue.push(rowUpper[row]);
	postValue.push(rowLower[row]);
}


void HPresolve::removeForcingConstraints(int mainIter) {
	double val;
	int i,j,k;
	for(i=0;i<numRow;i++) 
		if (flagRow.at(i)) {
			if (nzRow.at(i)==0) {
				removeEmptyRow(i);
				countRemovedRows[HTICK_PRE_EMPTY_ROW]++;
				continue;
			}

			//removeRowSingletons will handle just after removeForcingConstraints
			if (nzRow.at(i)==1)
				continue;

			timer.recordStart(HTICK_PRE_FORCING_ROW);

			double g=0;
			double h=0;
			k = ARstart[i+1]-1;
			while (k>=ARstart.at(i)) {
				j = ARindex.at(k);
				if (flagCol.at(j)) {
					if (ARvalue.at(k) < 0) {
						if (colUpper.at(j) < HSOL_CONST_INF)
							g+= ARvalue.at(k)*colUpper.at(j);
						else {
							g = -HSOL_CONST_INF;
							break;
						}
						
					}
					else {
						if (colLower.at(j) > -HSOL_CONST_INF)
							g+= ARvalue.at(k)*colLower.at(j);
						else {
							g = -HSOL_CONST_INF;
							break;
						}
					}	
				}
				k--;
			}
			k = ARstart[i+1]-1;
			while (k>=ARstart.at(i)) {
				j = ARindex.at(k);
				if (flagCol.at(j)) {
					if (ARvalue.at(k) < 0) {
						if (colLower.at(j) > -HSOL_CONST_INF)
							h+= ARvalue.at(k)*colLower.at(j);
						else {
							h = HSOL_CONST_INF;
							break;
						}
						
					}
					else {
						if (colUpper.at(j) < HSOL_CONST_INF)
							h+= ARvalue.at(k)*colUpper.at(j);
						else {
							h = HSOL_CONST_INF;
							break;
						}
					}	
				}
				k--;
			}
	
			timer.recordFinish(HTICK_PRE_FORCING_ROW);

			if (g>rowUpper.at(i) || h<rowLower.at(i)) {
				if (iPrint > 0)
					cout<<"PR: Problem infeasible."<<endl;
				status = Infeasible;
				return;
			}
			else if (g == rowUpper.at(i)) {
				//set all variables to lower bound
				if (iPrint > 0)
					cout<<"PR: Forcing row "<<i<<" removed. Following variables too:   nzRow="<<nzRow.at(i)<<endl;
				flagRow.at(i) = 0;
	        	addChange(3, i, 0);
				k = ARstart.at(i);
				while (k<ARstart[i+1]) {
					j = ARindex.at(k);
					if (flagCol.at(j)) {
						double value;
						if (ARvalue.at(k) < 0)
							value = colUpper.at(j);
						else
							value = colLower.at(j);
						setPrimalValue(j, value);
						valueColDual.at(j) = colCost.at(j);
						vector<double> bnds;
						bnds.push_back(colLower.at(j));
						bnds.push_back(colUpper.at(j));
						oldBounds.push(make_pair( j, bnds));
						addChange(2, 0, j);
						
						if (iPrint > 0)
							cout<<"PR:      Variable  "<<j<<" := "<<value<<endl;
						countRemovedCols[HTICK_PRE_FORCING_ROW]++;
					} k++; 
				}
				if (nzRow.at(i)==1)
					singRow.remove(i);
				countRemovedRows[HTICK_PRE_FORCING_ROW]++;
			}
			else if (h == rowLower.at(i)) {
				//set all variables to upper bound 
				if (iPrint > 0)
					cout<<"PR: Forcing row "<<i<<" removed. Following variables too:"<<endl;
				flagRow.at(i) = 0;
	        	addChange(3, i, 0);
				k = ARstart.at(i);
				while (k<ARstart[i+1]) {
					j = ARindex.at(k);
					if (flagCol.at(j)) {
						double value;
						if (ARvalue.at(k) < 0)
							value = colLower.at(j);
						else
							value = colUpper.at(j);
						setPrimalValue(j, value);
						valueColDual.at(j) = colCost.at(j);
						vector<double> bnds;
						bnds.push_back(colLower.at(j));
						bnds.push_back(colUpper.at(j));
						oldBounds.push(make_pair( j, bnds));
						addChange(2, 0, j);
						if (iPrint > 0)
							cout<<"PR:      Variable  "<<j<<" := "<<value<<endl;
						countRemovedCols[HTICK_PRE_FORCING_ROW]++;
					} k++;
				}	
				if (nzRow.at(i)==1)
					singRow.remove(i);
				countRemovedRows[HTICK_PRE_FORCING_ROW]++;
			}
			//redundant row: for any assignment of variables
			//constraint is satisfied
			else if (g >= rowLower.at(i) && h <= rowUpper.at(i)) {
				removeRow(i);
				addChange(16, i, 0);
				if (iPrint > 0)
					cout<<"PR: Redundant row "<<i<<" removed."<<endl;
				countRemovedRows[HTICK_PRE_REDUNDANT_ROW]++;
			}

			//Dominated constraints:
			else if (h < HSOL_CONST_INF) { 
				//fill in implied bounds arrays
				if (h < implRowValueUpper.at(i)) {
					implRowValueUpper.at(i) = h;} //	 cout<<"NEW UB row "<<i<< "iter = "<<mainIter<<endl; }
				if (h <= rowUpper.at(i))
					implRowDualLower.at(i) = 0;

				//calculate implied bounds for discovering free column singletons
				timer.recordStart(HTICK_PRE_DOMINATED_ROW_BOUNDS);
				for (k=ARstart.at(i); k<ARstart[i+1]; k++) {
					j = ARindex.at(k);
					if (flagCol.at(j)) {
						if (ARvalue.at(k) < 0 && colLower.at(j)> -HSOL_CONST_INF) {
							val =  (rowLower.at(i) - h)/ARvalue.at(k) + colLower.at(j);
							if (val<implColUpper.at(j)) {
								implColUpper.at(j) = val;
								implColUpperRowIndex.at(j) = i;
							}
						}
						else if (ARvalue.at(k) > 0 && colUpper.at(j) < HSOL_CONST_INF)  {
							 val = (rowLower.at(i) - h)/ARvalue.at(k) + colUpper.at(j);
							 if (val>implColLower.at(j)) {
							 	implColLower.at(j) = val;
							 	implColLowerRowIndex.at(j) = i;
							 }
						}
					}
				}
				timer.recordFinish(HTICK_PRE_DOMINATED_ROW_BOUNDS);
			}
			else if (g > -HSOL_CONST_INF) {
				//fill in implied bounds arrays
				if (g > implRowValueLower.at(i)) {
					implRowValueLower.at(i) = g; } //cout<<"NEW LB row "<<i<< "iter = "<<mainIter<<endl; }
				if (g >= rowLower.at(i))
					implRowDualUpper.at(i) = 0;

				//calculate implied bounds for discovering free column singletons
				timer.recordStart(HTICK_PRE_DOMINATED_ROW_BOUNDS);
				for (k=ARstart.at(i); k<ARstart[i+1]; k++) {
					j = ARindex.at(k);
					if (flagCol.at(j)) {
						if (ARvalue.at(k) < 0 && colUpper.at(j) < HSOL_CONST_INF) {
							val = (rowUpper.at(i) - g)/ARvalue.at(k) + colUpper.at(j);
							if (val > implColLower.at(j)) {
								implColLower.at(j) = val;
								implColLowerRowIndex.at(j) = i;
							}
						}
						else if (ARvalue.at(k) > 0 && colLower.at(j)> -HSOL_CONST_INF) {
							val = (rowUpper.at(i) - g)/ARvalue.at(k) + colLower.at(j);
							if (val < implColUpper.at(j)) {
								implColUpper.at(j) = val;
								implColUpperRowIndex.at(j) = i;
							}
						}
					}
				}
				timer.recordFinish(HTICK_PRE_DOMINATED_ROW_BOUNDS);
			}
		}

}					


void HPresolve::removeRowSingletons() {
	timer.recordStart(HTICK_PRE_SING_ROW);
	int i;
	while (!(singRow.empty()) ) {
		i=singRow.front();
		singRow.pop_front();

		//if (i==25 || i==22) 			continue; //on TEST

		if (!flagRow.at(i))	{
			cout<<"Warning: Row "<<i<<" already flagged off but in singleton row list. Ignored.\n";
			continue;
		}

		if (debug) {
			int vvv =  testSingRows(i);
			if (vvv < 0) {
				cout<<"Warning: Row "<<i<<" in singleton row list but not actually singleton. Ignored.\n";
				continue;
			}
		}
		

		int k = getSingRowElementIndexInAR(i); 
		int j = ARindex.at(k);

		//add old bounds OF X to checker and for postsolve
		if (iKKTcheck == 1) {
			vector<pair<int, double>> bndsL, bndsU, costS;
			bndsL.push_back( make_pair( j, colLower.at(j)));
			bndsU.push_back( make_pair( j, colUpper.at(j)));
			chk.cLowers.push(bndsL);
			chk.cUppers.push(bndsU);
		}
		vector<double> bnds;
		bnds.push_back(colLower.at(j));
		bnds.push_back(colUpper.at(j));
		bnds.push_back(rowLower.at(i));
		bnds.push_back(rowUpper.at(i));
		oldBounds.push(make_pair( j, bnds));

		double aij = ARvalue.at(k);
/*		//before update bounds of x take it out of rows with implied row bounds
		for (int r = Astart.at(j); r<Aend.at(j); r++) {
			if (flagRow[Aindex[r]]) {
				int rr = Aindex[r];
				if (implRowValueLower[rr] > -HSOL_CONST_INF) {
					if (aij > 0)
						implRowValueLower[rr] = implRowValueLower[rr] - aij*colLower.at(j);
					else
						implRowValueLower[rr] = implRowValueLower[rr] - aij*colUpper.at(j);
				}
				if (implRowValueUpper[rr] < HSOL_CONST_INF) {
					if (aij > 0)
						implRowValueUpper[rr] = implRowValueUpper[rr] - aij*colUpper.at(j);
					else
						implRowValueUpper[rr] = implRowValueUpper[rr] - aij*colLower.at(j);
				}
			}
		}*/


		//update bounds of X
		if (aij > 0) {
			if (rowLower.at(i) != -HSOL_CONST_INF)
				colLower.at(j) = max( max(rowLower.at(i)/aij, -HSOL_CONST_INF), colLower.at(j));
			if (rowUpper.at(i) != HSOL_CONST_INF)
				colUpper.at(j) = min( min(rowUpper.at(i)/aij, HSOL_CONST_INF), colUpper.at(j));
		}
		else if (aij < 0) {
			if (rowLower.at(i) != -HSOL_CONST_INF)
				colUpper.at(j) = min( min(rowLower.at(i)/aij, HSOL_CONST_INF), colUpper.at(j));
			if (rowUpper.at(i) != HSOL_CONST_INF)
				colLower.at(j) = max( max(rowUpper.at(i)/aij, -HSOL_CONST_INF), colLower.at(j));
		}

/*		//after update bounds of x add to rows with implied row bounds
		for (int r = Astart.at(j); r<Aend.at(j); r++) {
			if (flagRow[r]) {
				int rr = Aindex[r];
				if (implRowValueLower[rr] > -HSOL_CONST_INF) {
					if (aij > 0)
						implRowValueLower[rr] = implRowValueLower[rr] + aij*colLower.at(j);
					else
						implRowValueLower[rr] = implRowValueLower[rr] + aij*colUpper.at(j);
				}
				if (implRowValueUpper[rr] < HSOL_CONST_INF) {
					if (aij > 0)
						implRowValueUpper[rr] = implRowValueUpper[rr] + aij*colUpper.at(j);
					else
						implRowValueUpper[rr] = implRowValueUpper[rr] + aij*colLower.at(j);
				}
			}
		}*/

		//check for feasibility
		if (colLower.at(j) > colUpper.at(j) + tol) {
			status = Infeasible;
			return;
		}

		if (iPrint > 0)
					cout<<"PR: Singleton row "<<i<<" removed. Bounds of variable  "<<j<<" modified: l= "<<colLower.at(j) <<" u="<< colUpper.at(j) << ", aij = "<<aij<<endl;
		addChange(1, i, j);
		postValue.push(colCost.at(j));
		removeRow(i);


		if (flagCol.at(j) && colLower.at(j) == colUpper.at(j))
			removeIfFixed(j);

		countRemovedRows[HTICK_PRE_SING_ROW]++;

		}
	timer.recordFinish(HTICK_PRE_SING_ROW);
}

void HPresolve::addChange(int type, int row, int col) {
	change ch;
	ch.type = type;
	ch.row = row;
	ch.col = col;
	chng.push(ch);
}


//when setting a value to a primal variable and eliminating row update b, singleton Rows linked list, number of nonzeros in rows
void HPresolve::setPrimalValue(int j, double value) {
    flagCol.at(j) = 0;
    if (!hasChange)
	    hasChange = true;
	valuePrimal.at(j) = value;
	
	//update nonzeros
	for(int k=Astart.at(j);k<Aend.at(j);k++) {
			int i = Aindex.at(k);
			if (flagRow.at(i)) {
				nzRow.at(i)--;

			//update singleton row list
			if (nzRow.at(i)==1)
				singRow.push_back(i);
			if (nzRow.at(i)==0)
				singRow.remove(i);

			}
		}
	
	//update values if necessary 
	if (fabs(value)>0) {
		//RHS
		vector<pair<int, double>> bndsL, bndsU;
		for(int k=Astart.at(j);k<Aend.at(j);k++)
			if (flagRow[Aindex.at(k)]) {
				if (iKKTcheck == 1) {
		        	bndsL.push_back( make_pair( Aindex.at(k), rowLower[Aindex.at(k)]));
					bndsU.push_back( make_pair( Aindex.at(k), rowUpper[Aindex.at(k)]));
				}
				if (rowLower[Aindex.at(k)] > -HSOL_CONST_INF)
					rowLower[Aindex.at(k)] -= Avalue.at(k)*value;
				if (rowUpper[Aindex.at(k)] < HSOL_CONST_INF)
					rowUpper[Aindex.at(k)] -= Avalue.at(k)*value;
				if (implRowValueLower[Aindex.at(k)] > -HSOL_CONST_INF)
					implRowValueLower[Aindex.at(k)] -= Avalue.at(k)*value;
				if (implRowValueUpper[Aindex.at(k)] < HSOL_CONST_INF )
					implRowValueUpper[Aindex.at(k)] -= Avalue.at(k)*value;
			}
		if (iKKTcheck == 1) {
			chk.rLowers.push(bndsL);
			chk.rUppers.push(bndsU);
		}
		//shift objective 
		if (colCost.at(j) != 0)
			objShift += colCost.at(j)*value;
	}
}




void HPresolve::checkForChanges(int iteration) {
	int i;
	if (iteration==2) {
		bool allPresent = true;
		for (i=0;i<numCol;i++)
			if (!flagCol.at(i)) {
				allPresent = false;
				break;
			}
		for (i=0;i<numRow;i++)
			if (!flagRow.at(i)) {
				allPresent = false;
				break;
			}
		if (allPresent) {
			if (iPrint > 0)
				cout<<"PR: No variables were eliminated at presolve."<<endl;
			noPostSolve = true;
			return;
			}
		else
			resizeProblem();
	}
	else
		resizeProblem();
}


bool HPresolve::checkIfRedundant(int r, int type, double bound) { //>= 1, <= 2
	double rval = 0;
	
	for (int k = ARstart[r]; k<ARstart[r+1];k++)
		if ( (((type==1) && ARvalue.at(k)>0)  || ((type==2) && ARvalue.at(k)<0)))
			if (colLower[ARindex.at(k)]>-HSOL_CONST_INF)
				rval += colLower[ARindex.at(k)]*ARvalue.at(k) ;
			else return false;
		else if ( (((type==1) && ARvalue.at(k)<0)  || ((type==2) && ARvalue.at(k)>0)))
			if (colUpper[ARindex.at(k)]<HSOL_CONST_INF)
				rval += colUpper[ARindex.at(k)]*ARvalue.at(k) ;
			else return false;
	if (type==1 && rval >= bound)
		return true;
	if (type==2 && rval <= bound)
		return true;
	return false;
}





void HPresolve::reportTimes() {
	int reportList[] = {
				HTICK_PRE_EMPTY_ROW,
				HTICK_PRE_FIXED,
				HTICK_PRE_SING_ROW,
				HTICK_PRE_DOUBLETON_EQUATION,
				HTICK_PRE_FORCING_ROW,
				HTICK_PRE_REDUNDANT_ROW,
				HTICK_PRE_FREE_SING_COL,
				HTICK_PRE_SING_COL_DOUBLETON_INEQ,
				HTICK_PRE_IMPLIED_FREE_SING_COL,
				HTICK_PRE_DOMINATED_COLS,
				HTICK_PRE_WEAKLY_DOMINATED_COLS };
	int reportCount = sizeof(reportList) / sizeof(int);

	double totalTick = timer.getTick();
	printf("Presolve rules ");
	for (int i = 0; i < reportCount; i++) {
		printf(" %s", timer.itemNames[reportList[i]].c_str());
		cout<<flush;
	}

	printf("\n");
	cout<<"Time spent     "<<flush;
	for (int i = 0; i < reportCount; i++) {
		//int percent = 1000.0 * itemTicks[itemList.at(i)] / totalTick;
		float f = (float) timer.itemTicks[reportList[i]];
		if (f<0.01)
			cout<<setw(4)<<" <.01 ";
		else
			printf(" %3.2f ",f);
	}
	printf("\n");

}


void HPresolve::recordCounts(string fileName) {

	ofstream myfile;
	myfile.open(fileName, ios::app);
	int reportList[] = { HTICK_PRE_EMPTY_ROW, HTICK_PRE_FIXED,
			HTICK_PRE_SING_ROW, HTICK_PRE_DOUBLETON_EQUATION,
			HTICK_PRE_FORCING_ROW, HTICK_PRE_REDUNDANT_ROW,
			HTICK_PRE_FREE_SING_COL, HTICK_PRE_SING_COL_DOUBLETON_INEQ,
			HTICK_PRE_IMPLIED_FREE_SING_COL, HTICK_PRE_DOMINATED_COLS,
			HTICK_PRE_WEAKLY_DOMINATED_COLS, HTICK_PRE_EMPTY_COL };
	int reportCount = sizeof(reportList) / sizeof(int);

	myfile << "Problem " << modelName << ":\n";

	myfile << "Rule   , removed rows , removed cols , time  \n";

	int cRows=0, cCols=0;
	for (int i = 0; i < reportCount; i++) {
		float f = (float) timer.itemTicks[reportList[i]];

		myfile << setw(7) << timer.itemNames[reportList[i]].c_str() << ", "
				<< setw(7) << countRemovedRows[reportList[i]] << ", " << setw(7)
				<< countRemovedCols[reportList[i]] << ", ";
		if (f < 0.001)
			myfile << setw(7) << " <.001 ";
		else
			myfile << setw(7) << setprecision(3) << f;
		myfile << endl;

		cRows += countRemovedRows[reportList[i]];
		cCols += countRemovedCols[reportList[i]];

	}

	if (!noPostSolve) {
		if (cRows != numRowOriginal - numRow)
			cout<<"Wrong row reduction count\n";
		if (cCols != numColOriginal - numCol)
			cout<<"Wrong col reduction count\n";

		myfile << setw(7) << "Total " << ", " << setw(7) << numRowOriginal - numRow
				<< ", " << setw(7) << numColOriginal - numCol;
	}
	else {
		myfile << setw(7) << "Total " << ", " << setw(7) << 0
				<< ", " << setw(7) << 0;
	}
	myfile << endl << " \\\\ " << endl;
	myfile.close();
}

void HPresolve::resizeImpliedBounds() {
    //implied bounds for crashes
	//row duals
    vector<double> temp = implRowDualLower;
    vector<double> teup = implRowDualUpper;
    implRowDualLower.resize(numRow);
    implRowDualUpper.resize(numRow);
    int k, i;
    k=0;
    for (i=0;i<numRowOriginal;i++)
    	if (flagRow.at(i)) {
    		implRowDualLower.at(k) = temp.at(i);
    		implRowDualUpper.at(k) = teup.at(i);
    		k++;
	    }

    //row value
    temp = implRowValueLower;
    teup = implRowValueUpper;
    implRowValueLower.resize(numRow);
    implRowValueUpper.resize(numRow);
    k=0;
    for (i=0;i<numRowOriginal;i++)
    	if (flagRow.at(i)) {
    		if (temp.at(i) < rowLower.at(i))
    			temp.at(i) = rowLower.at(i);
    		implRowValueLower.at(k) = temp.at(i);
    		if (teup.at(i) > rowUpper.at(i))
    			teup.at(i) = rowUpper.at(i);
    		implRowValueUpper.at(k) = teup.at(i);
    		k++;
	    }

    //column dual
    temp = implColDualLower;
    teup = implColDualUpper;
    implColDualLower.resize(numCol);
    implColDualUpper.resize(numCol);
    k=0;
    for (i=0;i<numColOriginal;i++)
    	if (flagCol.at(i)) {
    		implColDualLower.at(k) = temp.at(i);
    		implColDualUpper.at(k) = teup.at(i);
    		k++;
	    }

    //column value
    temp = implColLower;
    teup = implColUpper;
    implColLower.resize(numCol);
    implColUpper.resize(numCol);
    k=0;
    for (i=0;i<numColOriginal;i++)
    	if (flagCol.at(i)) {
    		if (temp.at(i) < colLower.at(i))
    			temp.at(i) = colLower.at(i);
    		implColLower.at(k) = temp.at(i);
    		if (teup.at(i) > colUpper.at(i))
    			teup.at(i) = colUpper.at(i);
    		implColUpper.at(k) = teup.at(i);
    		k++;
	    }
}

int HPresolve::getSingRowElementIndexInAR(int i) {
	int k=ARstart.at(i);
    while (!flagCol[ARindex.at(k)])
            k++;  
    return k;
}

int HPresolve::getSingColElementIndexInA(int j) {
	int k=Astart.at(j);
    while (!flagRow[Aindex.at(k)])
            k++;  
    return k;
}

int HPresolve::testSingRows(int i) {
	int out = -1;

	int k=ARstart.at(i);
    while (!flagCol[ARindex.at(k)])
            k++;
            
	char buff[5];
	if (k>=ARstart[i+1]) { 
		if (1) {
			cout<< "all flagCol are false. row "<<i<<endl;
		
			cout<< "Avalue: ";
			for (int j=ARstart.at(i); j< ARstart[i+1]; j++) {
				sprintf(buff, "%2.1f ", ARvalue.at(j));
				cout<<setw(5)<<buff; 
				
			}
			cout<<endl<<"Column: ";
			for (int j=ARstart.at(i); j< ARstart[i+1]; j++) {
				sprintf(buff, "%i ", ARindex.at(j));
				cout<<setw(5)<<buff; 
			}
			cout<<endl<<"isValid: ";
			for (int j=ARstart.at(i); j< ARstart[i+1]; j++) {
				//sprintf(buff, "%i ", flagCol.at(j));
				cout<<setw(5)<<flagCol.at(j);
			}
			cout<<endl;
		}
		return -3;	
	}
	else 
		out = k;
	k++;
	while (k<ARstart[i+1] && !flagCol[ARindex.at(k)] )
		k++;
	if (k<ARstart[i+1]) {
		cout<< "more than one flagCol is true."<<endl;
		
		cout<< "Avalue: ";
		for (int j=ARstart.at(i); j< ARstart[i+1]; j++) {
			sprintf(buff, "%2.1f ", ARvalue.at(j));
    		cout<<setw(5)<<buff; 
    		
		}
		cout<< "Avalue: ";
		for (int j=ARstart.at(i); j< ARstart[i+1]; j++) {
			sprintf(buff, "%2.1f ", ARvalue.at(j));
    		cout<<setw(5)<<buff; 
		}
		cout<<endl<<"Column: ";
		for (int j=ARstart.at(i); j< ARstart[i+1]; j++) {
			sprintf(buff, "%i ", ARindex.at(j));
    		cout<<setw(5)<<buff; 
		}
		cout<<endl<<"isValid: ";
		for (int j=ARstart.at(i); j< ARstart[i+1]; j++) {
			//sprintf(buff, "%i ", flagCol.at(j));
    		cout<<setw(5)<<flagCol.at(j);
		}
		cout<<endl;
		return -2;
	}   
	if (out>=0) 
		return ARindex[out];
	return out;
}


int HPresolve::testSingCols(int i) {
	int out = -1;
	int k=Astart.at(i);
    while (!flagRow[Aindex.at(k)])
            k++;
            
	if (k>=Aend.at(i)) {
		cout<< "all flagRow are false."<<endl;
		return -3;	
	}
	else 
		out = k;
	k++;
	while (k<Aend.at(i) && !flagRow[Aindex.at(k)] )
		k++;
	if (k<Aend.at(i)) {
		cout<< "more than one flagRow is true. j="<<i<<endl;
		return -2;
	}    
	if (out>=0) 
		return Aindex[out];
	return out;
}

void HPresolve::testAnAR(int post) {
	int rows = numRow;
	int cols = numCol;
	int i,j,k;

	double valueA, valueAR;
	bool hasValueA, hasValueAR;

	if (post) {
		rows = numRowOriginal;
		cols = numColOriginal;
	}

	// check that A = AR
	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++) {
			if (post==0)
				if (!flagRow.at(i) || !flagCol.at(j))
				continue;
			hasValueA = false;
			for (k = Astart.at(j); k<Aend.at(j); k++)
				if (Aindex.at(k) == i) {
					hasValueA = true;
					valueA = Avalue.at(k);
				}

			hasValueAR = false;
			for (k = ARstart.at(i); k<ARstart[i+1]; k++)
				if (ARindex.at(k) == j) {
					hasValueAR = true;
					valueAR = ARvalue.at(k);
				}

			if (hasValueA != hasValueAR)
				cout<<"    MATRIX is0   DIFF row="<<i<< " col="<<j<<"           ------------A: "<<hasValueA<<"  AR: "<<hasValueAR <<endl;
			else if (hasValueA && valueA != valueAR)
				cout<<"    MATRIX VAL  DIFF row="<<i<< " col="<<j<<"           ------------A: "<<valueA<<"  AR: "<<valueAR <<endl;
		}
	}

	if (post == 0) {
		//check nz
		int nz=0;
		for (i=0; i<rows; i++) {
			if (!flagRow.at(i))
				continue;
			nz=0;
			for (k = ARstart.at(i); k<ARstart[i+1]; k++)
				if (flagCol[ARindex.at(k)])
					nz++;
			if (nz != nzRow.at(i))
				cout<<"    NZ ROW      DIFF row="<<i<< " nzRow="<<nzRow.at(i)<<" actually "<<nz <<"------------"<<endl;
		}

		for (j=0; j<cols; j++) {
			if (!flagCol.at(j))
				continue;
			nz=0;
			for (k = Astart.at(j); k<Aend.at(j); k++)
				if (flagRow[Aindex.at(k)])
					nz++;
			if (nz != nzCol.at(j))
				cout<<"    NZ COL      DIFF col="<<j<< " nzCol="<<nzCol.at(j)<<" actually "<<nz <<"------------"<<endl;
		}
	}

}


void HPresolve::postsolve() {

		if (noPostSolve) {
			//set valuePrimal
			for (int i=0;i<numCol;i++) {
				valuePrimal.at(i) = colValue.at(i);
				valueColDual.at(i) = colDual.at(i);
			}
			for (int i=0;i<numRow;i++)
				valueRowDual.at(i) = rowDual.at(i);
			//For KKT check: first check solverz` results before we do any postsolve
			if (iKKTcheck == 1) {
				chk.passSolution(colValue, colDual,  rowDual);
				chk.makeKKTCheck();
			}
			//testBasisMatrixSingularity();
			return;
		}

		//For KKT check: first check solver results before we do any postsolve
		if (iKKTcheck == 1) {
			cout<<"----KKT check on hsol solution-----\n";

			chk.passSolution(colValue, colDual,  rowDual);
			chk.makeKKTCheck();
		}
		//So there have been changes definitely ->
		makeACopy(); // so we can efficiently calculate primal and dual values

		//	iKKTcheck = false;
		//set corresponding parts of solution vectors:
		int j=0;
		vector<int> eqIndexOfReduced(numCol, -1);
		vector<int> eqIndexOfReduROW(numRow, -1);
		for (int i=0;i<numColOriginal;i++)
			if (cIndex.at(i)>-1) {
				eqIndexOfReduced.at(j) = i;
				j++;
			}
		j=0;
		for (int i=0;i<numRowOriginal;i++)
			if (rIndex.at(i)>-1) {
				eqIndexOfReduROW.at(j) = i;
				j++;
			}

		vector<int> temp = nonbasicFlag;

		nonbasicFlag.assign(numColOriginal + numRowOriginal, 1);

		for (int i=0;i<numCol;i++) {
			valuePrimal[eqIndexOfReduced.at(i)] = colValue.at(i);
			valueColDual[eqIndexOfReduced.at(i)] = colDual.at(i);
			nonbasicFlag[eqIndexOfReduced.at(i)] = temp.at(i);
		}

		for (int i=0;i<numRow;i++) {
			valueRowDual[eqIndexOfReduROW.at(i)] = rowDual.at(i);
			nonbasicFlag[numColOriginal + eqIndexOfReduROW.at(i)] = temp[numCol + i];
		}

		//cmpNBF(-1, -1);

		int kk, jj;
		double y,z,x;
		vector<int> fRjs;
		while (!chng.empty()) {
			change c = chng.top();
			chng.pop();
			//cout<<"chng.pop:       "<<c.col<<"       "<<c.row << endl;

			setBasisElement(c);
			switch(c.type) {
				case 17: { //Doubleton equation row
					getDualsDoubletonEquation(c.row,c.col);

					if (iKKTcheck == 1) {
						if (chk.print == 1)
							cout<<"----KKT check after doubleton equation re-introduced. Row: "<< c.row <<", column "<< c.col <<" -----\n";
						chk.addChange(17, c.row, c.col, valuePrimal[c.col], valueColDual[c.col], valueRowDual[c.row]);
						chk.makeKKTCheck();
					}
					//exit(2);
					break;
				}
				case 171: { //new bounds from doubleton equation, retrieve old ones
					//just for KKT check, not called otherwise
					chk.addChange(171, c.row, c.col, 0, 0, 0);
					break;
				}
				case 172: { //matrix transformation from doubleton equation, case x still there
					//case new x is not 0
					//just change value of entry in row for x

					int indi;
					for (indi=ARstart[c.row];indi<ARstart[c.row+1];indi++)
							if (ARindex[indi]==c.col)
								break;
					ARvalue[indi] = postValue.top();
					for (indi=Astart[c.col];indi<Aend[c.col];indi++)
							if (Aindex[indi]==c.row)
								break;
					Avalue[indi] = postValue.top();

					if (iKKTcheck == 1)
						chk.addChange(172, c.row, c.col, postValue.top(), 0, 0);
					postValue.pop();

					break;
				}
			case 173: { //matrix transformation from doubleton equation, retrieve old value
				//case when row does not have x initially: entries for row i swap x and y cols

				int indi, yindex;
				yindex = (int) postValue.top();
				postValue.pop();

				//reverse AR for case when x is zero and y entry has moved
				for (indi=ARstart[c.row];indi<ARstart[c.row+1];indi++)
						if (ARindex[indi]==c.col)
							break;
				ARvalue[indi] = postValue.top();
				ARindex[indi] = yindex;

				//reverse A for case when x is zero and y entry has moved
				for (indi=Astart[c.col];indi<Aend[c.col];indi++)
						if (Aindex[indi]==c.row)
							break;

				//recover x: column decreases by 1
				//if indi is not Aend-1 swap elements indi and Aend-1
				if (indi != Aend[c.col]-1) {
					double tmp = Avalue[Aend[c.col]-1];
					int   tmpi = Aindex[Aend[c.col]-1];
					Avalue[Aend[c.col]-1] = Avalue[indi];
					Aindex[Aend[c.col]-1] = Aindex[indi];
					Avalue[indi] = tmp;
					Aindex[indi] = tmpi;
				}
				Aend[c.col]--;

				//recover y: column increases by 1
				//update A: append X column to end of array
				int st = Avalue.size();
				for (int ind = Astart[yindex]; ind < Aend[yindex]; ind++) {
					Avalue.push_back(Avalue[ind]);
					Aindex.push_back(Aindex[ind]);
				}
				Avalue.push_back(postValue.top());
				Aindex.push_back(c.row);
				Astart[yindex] = st;
				Aend[yindex] = Avalue.size();

				if (iKKTcheck == 1)
					chk.addChange(173, c.row, c.col, postValue.top(), (double) yindex, 0);
				postValue.pop();

				break;
			}
			case 174: { //sp case x disappears row representation change
				int indi;
				for (indi=ARstart[c.row];indi<ARstart[c.row+1];indi++)
						if (ARindex[indi]==numColOriginal)
							break;
				ARindex[indi] = c.col;
				ARvalue[indi] = postValue.top();

				if (iKKTcheck == 1) {
					chk.ARindex[indi] = c.col;
					chk.ARvalue[indi] = postValue.top();
				}

				postValue.pop();

				break;
			}
			case 175: { //sp case x disappears column representation change
				//here A is copied from AR array at end of presolve so need to expand x column
				//Aend[c.col]++; wouldn't do because old value is overriden
				double oldXvalue = postValue.top();postValue.pop();
				int x = c.col;

				//update A: append X column to end of array
				int st = Avalue.size();
				for (int ind = Astart[x]; ind < Aend[x]; ind++) {
					Avalue.push_back(Avalue[ind]);
					Aindex.push_back(Aindex[ind]);
				}
				Avalue.push_back(oldXvalue);
				Aindex.push_back(c.row);
				Astart[x] = st;
				Aend[x] = Avalue.size();

				break;
			}
				case 0: {//empty row
					valueRowDual[c.row] = 0;
					flagRow[c.row] = 1;
					if (iKKTcheck == 1) {
						if (chk.print == 1)
							cout<<"----KKT check after empty row "<<c.row  <<" re-introduced-----\n";
						chk.addChange(0, c.row, 0, 0, 0, 0);
						chk.makeKKTCheck();
						}
					break;
				}
				case 1: {//sing row
					//valuePrimal is already set for this one, colDual also, we need rowDual. AR copy keeps full matrix.
					//col dual maybe infeasible, we need to check.
					//recover old bounds and see
					getDualsSingletonRow(c.row, c.col);

					if (iKKTcheck == 1) {
						if (chk.print == 1)
							cout<<"----KKT check after singleton row "<<c.row <<" re-introduced. Variable: "<<c.col<<" -----\n";
						chk.addChange(1, c.row, c.col, valuePrimal[c.col], valueColDual[c.col], valueRowDual[c.row]);
						chk.makeKKTCheck();
						}
					break;
				}
				case 2: //variables set at a bound by forcing row
					fRjs.push_back(c.col);
					flagCol[c.col] = 1;
					if (iKKTcheck == 1 && valuePrimal[c.col] != 0)
						chk.addChange(22, c.row, c.col, 0, 0, 0);
					break;
				case 3: {//forcing row

					string str = getDualsForcingRow( c.row, fRjs);

					if (iKKTcheck == 1) {
						if (chk.print == 1)
							cout<<"----KKT check after forcing row "<<c.row<<" re-introduced. Variable(s): "<<str<<" -----\n";
						chk.addChange(3, c.row, 0, 0, 0, valueRowDual[c.row]);
						chk.makeKKTCheck();
						}
					fRjs.clear();
					break;
				}
				case 16: {//redundant row
					valueRowDual[c.row] = 0;

					flagRow[c.row] = 1;

					if (iKKTcheck == 1) {
						if (chk.print == 1)
							cout<<"----KKT check after redundant row "<<c.row<<" re-introduced.----------------\n";
						chk.addChange(0, c.row, 0, 0, 0, 0);
						chk.makeKKTCheck();
						}
					break;
				}
				case 4 : case 8: { //implied free
					//colDual rowDual already set.

					//calculate row value without xj
					double aij = getaij(c.row,c.col);
					double sum = 0;
					for (int k=ARstart[c.row]; k<ARstart[c.row+1]; k++)
						if (flagCol[ARindex.at(k)])
							sum += valuePrimal[ARindex.at(k)]*ARvalue.at(k);

					double rowlb = postValue.top(); postValue.pop();
					double rowub = postValue.top(); postValue.pop();

					//calculate xj
					if (valueRowDual[c.row] < 0) {
						//row is at lower bound
						valuePrimal[c.col] = (rowlb - sum )/aij;
					}
					else if (valueRowDual[c.row] > 0) {
						//row is at upper bound
						valuePrimal[c.col] = (rowub - sum )/aij;
					}
					else if (rowlb == rowub )
						valuePrimal[c.col] = (rowlb - sum )/aij;
					else if (colCostAtEl[c.col] > 0) {
						//we are interested in the lowest possible value of x:
						//max { l_j, bound implied by row i }
						double bndL;
						if (aij>0)
							bndL = (rowlb - sum )/aij;
						else
							bndL = (rowub - sum )/aij;
						valuePrimal[c.col] = max(colLower[c.col], bndL);
					}
					else if (colCostAtEl[c.col] < 0) {
						//we are interested in the highest possible value of x:
						//min { u_j, bound implied by row i }
						double bndU;
						if (aij<0)
							bndU = (rowlb - sum )/aij;
						else
							bndU = (rowub - sum )/aij;
						valuePrimal[c.col] = min(colUpper[c.col], bndU);
					}
					else { //cost is zero
						double bndL, bndU;
						if (aij>0) {
							bndL = (rowlb - sum )/aij;
							bndU = (rowub - sum )/aij;
						}
						else {
							bndL = (rowub - sum )/aij;
							bndU = (rowlb - sum )/aij;
						}
						double valuePrimalUB = min(colUpper[c.col], bndU);
						double valuePrimalLB = max(colLower[c.col], bndL);
						if (abs(valuePrimalLB) < abs(valuePrimalUB))
							valuePrimal[c.col] = valuePrimalLB;
						else
							valuePrimal[c.col] = valuePrimalUB;
					}
					sum = sum + valuePrimal[c.col] * aij;

					double costAtTimeOfElimination = postValue.top(); postValue.pop();
					objShift += (costAtTimeOfElimination* sum)/aij;

					flagRow[c.row] = 1;
					flagCol[c.col] = 1;
					//valueRowDual[c.row] = 0;

					if (iKKTcheck == 1) {
						chk.addCost(c.col, costAtTimeOfElimination);
						if (c.type == 4 && chk.print == 1)
							cout<<"----KKT check after free col singleton "<<c.col <<" re-introduced. Row: "<<c.row<<" -----\n";
						else if (c.type == 8 && chk.print == 1)
							cout<<"----KKT check after implied free col singleton "<<c.col <<" re-introduced. Row: "<<c.row<<" -----\n";
						chk.addChange(4, c.row, c.col, valuePrimal[c.col], valueColDual[c.col], valueRowDual[c.row]);
						chk.makeKKTCheck();
						}
					break;
				}
				case 5: {
					// column singleton in a doubleton equation.
					// colDual already set. need valuePrimal from stack. maybe change rowDual depending on bounds. old bounds kept in oldBounds.
					// variables j,k : we eliminated j and are left with changed bounds on k and no row.
					// c.col is column COL (K) - eliminated, j is with new bounds
					pair< int ,vector<double>> p = oldBounds.top(); oldBounds.pop();
					vector<double> v = get<1>(p);
					int j = get<0>(p);
					double ubNew = v[1];
					double lbNew = v[0];
					double cjNew = v[2];
					p = oldBounds.top(); oldBounds.pop();
					v = get<1>(p);
					double ubOld = v[1];
					double lbOld = v[0];
					double cjOld = v[2];
					p = oldBounds.top(); oldBounds.pop();
					v = get<1>(p);
					double ubCOL = v[1];
					double lbCOL = v[0];
					double ck = v[2];

					double rowlb = postValue.top(); postValue.pop();
					double rowub = postValue.top(); postValue.pop();
					double aik = postValue.top(); postValue.pop();
					double aij = postValue.top(); postValue.pop();
					double xj = valuePrimal.at(j);

					//calculate xk, depending on signs of coeff and cost
					double upp = HSOL_CONST_INF;
					double low = -HSOL_CONST_INF;

					if (( aij > 0 && aik > 0 ) || ( aij < 0 && aik < 0 )) {
						upp = (rowub - aij*xj) / aik;
						low = (rowlb - aij*xj) / aik;
					}
					else {
						upp = (rowub - aij*xj) / aik;
						low = (rowlb - aij*xj) / aik;
					}

					double xkValue;
					if (ck == 0) {
						if (low < 0 && upp > 0)
							xkValue = 0;
						else if (abs(low) < abs(upp))
							xkValue = low;
						else
							xkValue = upp;
					}

					else if ((ck > 0 && aik > 0) || (ck < 0 && aik < 0)) {
						if (low == -HSOL_CONST_INF)
							cout<<"ERROR UNBOUNDED? unnecessary check";
						xkValue = low;
					}
					else if ((ck > 0 && aik < 0) || (ck < 0 && aik > 0)) {
						if (upp == HSOL_CONST_INF)
							cout<<"ERROR UNBOUNDED? unnecessary check";
						xkValue = upp;
					}

					// primal value and objective shift
					valuePrimal[c.col] = xkValue;
					objShift += -cjNew*xj + cjOld*xj + ck*xkValue;

					//fix duals

					double rowVal = aij* xj + aik*xkValue;
					if (rowub - rowVal > tol && rowVal - rowlb > tol) {
						valueRowDual[c.row] = 0;
						flagRow[c.row] = 1;
						valueColDual.at(j) = getColumnDualPost(j);
						valueColDual[c.col] = getColumnDualPost(c.col);
					}
					else {
						double lo,up;
						if (abs(rowlb- rowub) < tol) {
							lo = -HSOL_CONST_INF;
							up = HSOL_CONST_INF;
						}
						else if (abs(rowub - rowVal) <= tol ) {
							lo = 0;
							up = HSOL_CONST_INF;
						}
						else if (abs(rowlb - rowVal) <= tol ) {
							lo = -HSOL_CONST_INF;;
							up = 0;
						}

						colCostAtEl.at(j) = cjOld; //revert cost before calculating duals
						getBoundOnLByZj(c.row, j, &lo, &up, lbOld, ubOld);
						getBoundOnLByZj(c.row, c.col, &lo, &up, lbCOL,    ubCOL);

							//calculate yi
						if (lo-up > tol)
							cout<<"PR: Error in postsolving doubleton inequality "<<c.row <<" : inconsistent bounds for it's dual value.\n";

						if (lo<=0 && up>=0) {
							valueRowDual[c.row] = 0;
						}
						else if (lo>0) {
							valueRowDual[c.row] = lo;
						}
						else if (up<0) {
							valueRowDual[c.row] = up;
						}

						flagRow[c.row] = 1;
						valueColDual.at(j) = getColumnDualPost(j);
						if (iKKTcheck == 1)
							chk.colDual.at(j) = valueColDual.at(j);

						valueColDual[c.col] = getColumnDualPost(c.col);

					}

					if (abs(valueColDual[c.col]) > tol) {
						nonbasicFlag[c.col] = 1;
						nonbasicFlag.at(j) = 0;
						basicIndex.pop_back();
						basicIndex.push_back(j);
					}


					flagRow[c.row] = 1;
					flagCol[c.col] = 1;

					if (iKKTcheck == 1) {
						if (chk.print == 1)
							cout<<"----KKT check after col singleton "<<c.col <<" in doubleton eq re-introduced. Row: "<<c.row<<" -----\n";

						chk.addChange(5, c.row, c.col, valuePrimal[c.col], valueColDual[c.col], valueRowDual[c.row]);

						chk.makeKKTCheck();
					}
					//exit(2);
					break;
				}
				case 6: case 9: case 10: {//empty column: got valuePrimal, need colDual, also dominated column and weakly dominated column
					if (c.type != 6) {
						z = colCostAtEl[c.col];
						for (int k=Astart[c.col]; k<Astart[c.col+1]; k++)
							if (flagRow[Aindex.at(k)])
								z = z + valueRowDual[Aindex.at(k)]*Avalue.at(k);
						valueColDual[c.col] = z;
					}


					flagCol[c.col] = 1;
					if (iKKTcheck == 1) {
						if (c.type == 6 && chk.print == 1)
							cout<<"----KKT check after empty column "<<c.col <<" re-introduced.-----------\n";
						else if (c.type == 9 && chk.print == 1)
							cout<<"----KKT check after dominated column "<<c.col <<" re-introduced.-----------\n";
						else if (c.type == 10 && chk.print == 1)
							cout<<"----KKT check after weakly dominated column "<<c.col <<" re-introduced.-----------\n";

						chk.addChange(6, 0, c.col, valuePrimal[c.col], valueColDual[c.col], 0);
						chk.makeKKTCheck();
					}
					break;
				}

				case 7: { //fixed variable: got valuePrimal, need colDual

					valueColDual[c.col] = getColumnDualPost(c.col);

					flagCol[c.col] = 1;
					if (iKKTcheck == 1) {
						if (chk.print == 1)
							cout<<"----KKT check after fixed variable "<<c.col <<" re-introduced.-----------\n";
						chk.addChange(7, 0, c.col, valuePrimal[c.col], valueColDual[c.col], 0);
						chk.makeKKTCheck();
					}
					break;
				}
				case 11: {//empty row from duplucate rows
					valueRowDual[c.row] = 0;
					flagRow[c.row] = 1;

					//check duals
					pair< int ,vector<double>> p = oldBounds.top(); oldBounds.pop();
					vector<double> vv = get<1>(p);
					double ubOld = vv[1];
					double lbOld = vv[0];
					double ubNew = vv[2];
					double lbNew = vv[3];

					//then the linear transformation of A
					double v = postValue.top(); postValue.pop();
					int i = (int) postValue.top(); postValue.pop();
					int k = (int) postValue.top(); postValue.pop();

					if (k!=c.row)
						cout<<"PR: Error in postsolving implied free row after duplicate row "<<i<<" added to row "<<k<<" "<<v<<" times. ";

					double rv = getRowValue(i);

					if ((valueRowDual.at(i) != 0 && rv == lbNew && lbNew > lbOld) ||
						(valueRowDual.at(i) != 0 && rv == ubNew && ubNew < ubOld)) {
						valueRowDual.at(k) = valueRowDual.at(i)/v;
						valueRowDual.at(i) = 0;
						if (iKKTcheck == 1)
							chk.addChange(21, i, 0, 0, 0, 0);
					}

					if (iKKTcheck == 1) {
						if (chk.print == 1)
							cout<<"----KKT check after empty row "<<c.row <<" from duplucate rows re-introduced-----\n";
						chk.addChange(11, c.row, 0, 0, 0, valueRowDual.at(k));
						chk.makeKKTCheck();
						}
					break;
				}

			}
			//cmpNBF(c.row, c.col);
		}

		//cmpNBF();
		//cout<<"Singularity check at end of postsolve: ";
		//testBasisMatrixSingularity();

		int cntVar = 0;
		basicIndex.resize(numRowOriginal);
		for (int i=0; i< numColOriginal + numRowOriginal; i++) {
			if (nonbasicFlag.at(i) == 0) {
				if (cntVar==numRowOriginal) {
					cout<<"Error during postsolve: wrong basic variables number: basic variables more than rows."<<endl;
					break;
				}
				basicIndex[cntVar] = i;
				cntVar++;
			}
		}
		if (cntVar<numRowOriginal)
			cout<<"Error during postsolve: wrong basic variables number: basic variables less than rows."<<endl;


		if (iKKTcheck == 2) {
			if (chk.print==3)
				chk.print=2;
			chk.passSolution(valuePrimal, valueColDual,  valueRowDual);
			chk.makeKKTCheck();
		}

		// now recover original model data to pass back to hsol
		// A is already recovered!
		// however, A is expressed in terms of Astart, Aend and columns are in different order so
		makeACopy();

		numRow = numRowOriginal;
		numCol = numColOriginal;
		numTot = numRow + numCol;

		rowUpper = rowUpperOriginal;
		rowLower = rowLowerOriginal;

		colUpper = colUpperOriginal;
		colLower = colLowerOriginal;

		colCost = colCostOriginal;

		nonbasicMove.resize(numTot, 0);
		for (int i=0; i< numColOriginal; i++) {
			if (colLower.at(i) != colUpper.at(i) && colLower.at(i) != -HSOL_CONST_INF )
				nonbasicMove.at(i) = 1;
			else if (colUpper.at(i) != HSOL_CONST_INF)
				nonbasicMove.at(i) = -1;
			else
				nonbasicMove.at(i) = 0;
		}

}




void HPresolve::setBasisElement(change c) {

        //nonbasicFlag starts off as [numCol + numRow] and is already
        //increased to [numColOriginal + numRowOriginal] so fill in gaps

        switch (c.type) {
                        case 0: {//empty row
                                nonbasicFlag[numColOriginal + c.row] = 0;
                                break;
                        }
                        case 1: {//sing row : elsewhere
                                break;
                        }
                        case 2: // variables set at a bound by forcing row fRjs.p
                        		// all nonbasic

                                break;
                        case 3: {//forcing row
                        		nonbasicFlag[numColOriginal + c.row] = 0;
                                break;
                        }
                        case 16: {//redundant row
                                nonbasicFlag[numColOriginal + c.row] = 0;
                                break;
                        }
                        case 4 : case 8: { //(implied) free
                                nonbasicFlag[c.col] = 0;
                                basicIndex.push_back(c.col);

                                nonbasicFlag[numColOriginal + c.row] = 1;
                                break;
                        }
                        case 5: {
                                // column singleton in a doubleton inequality.
                                nonbasicFlag[c.col] = 0;
                                basicIndex.push_back(c.col);

                                nonbasicFlag[numColOriginal + c.row] = 1;
                                break;
                        }
                        case 6: case 9: case 10: {//empty column, also dominated
                                nonbasicFlag[c.col] = 1;
                                break;
                        }
                        case 7: { //fixed variable:
                        		//check if it was NOT after singRow
                        		if (chng.size() > 0)
									if (chng.top().type!=1)
										nonbasicFlag[c.col] = 1;
                                break;
                        }

        }

}

/* testing and dev
int HPresolve::testBasisMatrixSingularity() {

	HFactor factor;

	//resize matrix in M so we can pass to factor
	int i, j, k;
	int nz = 0;
	int nR = 0;
	int nC = 0;

	numRowOriginal = rowLowerOriginal.size();
	numColOriginal = colLowerOriginal.size();
	//arrays to keep track of indices
	vector<int> rIndex_(numRowOriginal, -1);
	vector<int> cIndex_(numColOriginal, -1);

	for (i=0;i<numRowOriginal;i++)
		if (flagRow.at(i)) {
			for (j = ARstart.at(i); j<ARstart[i+1]; j++)
				if (flagCol[ARindex.at(j)])
					nz ++;
			rIndex_.at(i) = nR;
			nR++;
			}

	for (i=0;i<numColOriginal;i++)
		if (flagCol.at(i)) {
			cIndex_.at(i) = nC;
			nC++;
		}


	//matrix
	vector<int>    Mstart(nC + 1, 0);
	vector<int>    Mindex(nz);
	vector<double> Mvalue(nz);

    vector<int> iwork(nC, 0);

    for (i = 0;i<numRowOriginal; i++)
    	if (flagRow.at(i))
    	    for (int k = ARstart.at(i); k < ARstart[i+1]; k++) {
    	    	j = ARindex.at(k);
    	    	if (flagCol.at(j))
        			iwork[cIndex_.at(j)]++;
        		}
    for (i = 1; i <= nC; i++)
        Mstart.at(i) = Mstart[i - 1] + iwork[i - 1];
   for (i = 0; i < numColOriginal; i++)
        iwork.at(i) = Mstart.at(i);

   for (i = 0; i < numRowOriginal; i++) {
    	if (flagRow.at(i)) {
			int iRow = rIndex_.at(i);
		    for (k = ARstart.at(i); k < ARstart[i + 1]; k++) {
		        j = ARindex.at(k);
		        if (flagCol.at(j)) {
		        	int iCol = cIndex_.at(j);
				    int iPut = iwork[iCol]++;
				    Mindex[iPut] = iRow;
				    Mvalue[iPut] = ARvalue.at(k);
				}
		    }
		}
    }

    vector<int>  bindex(nR);
    int countBasic=0;

     for (int i=0; i< nonbasicFlag.size();i++) {
    	 if (nonbasicFlag.at(i) == 0)
    			 countBasic++;
     }

     if (countBasic != nR)
    	 cout<<" Wrong count of basic variables: != numRow"<<endl;

     int c=0;
     for (int i=0; i< nonbasicFlag.size();i++) {
    	 if (nonbasicFlag.at(i) == 0) {
			if (i < numColOriginal)
				bindex[c] = cIndex_.at(i);
			else
				bindex[c] = nC + rIndex_[i - numColOriginal];
			c++;
    	 }
    }

	factor.setup(nC, nR, &Mstart[0], &Mindex[0], &Mvalue[0],  &bindex[0]);
/*	if (1) // for this check both A and M are the full matrix again
	{
		if (nC - numColOriginal != 0)
			cout<<"columns\n";
		if (nR - numRowOriginal != 0)
			cout<<"rows\n";
		for (int i=0; i< Mstart.size();i++)
			if (Mstart.at(i) - Astart.at(i) != 0)
				cout<<"Mstart "<<i<<"\n";
		for (int i=0; i< Mindex.size();i++)
			if (Mindex.at(i) - Aindex.at(i) != 0)
				cout<<"Mindex "<<i<<"\n";
		for (int i=0; i< Mvalue.size();i++)
			if (Mvalue.at(i) - Avalue.at(i) != 0)
				cout<<"Mvalue "<<i<<"\n";
		for (int i=0; i< bindex.size();i++)
			if (nonbasicFlag.at(i) - nbffull.at(i) != 0)
				cout<<"nbf "<<i<<"\n";
	} * /

	try {
        factor.build();
    } catch (runtime_error& error) {
        cout << error.what() << endl;
        cout << "Postsolve: could not factorize basis matrix." << endl;
        return 0;
    }
    cout << "Postsolve: basis matrix successfully factorized." << endl;

    return 1;
}*/




/***
 * lo and up refer to the place storing the current bounds on y_row
 *
 */
void HPresolve::getBoundOnLByZj(int row, int j, double* lo, double* up, double colLow, double colUpp) {

	double cost = colCostAtEl.at(j);//valueColDual.at(j);
	double x = -cost;

	double sum = 0;
	for (int kk=Astart.at(j); kk<Aend.at(j);kk++)
		if (flagRow[Aindex[kk]]) {
			sum = sum + Avalue[kk]*valueRowDual[Aindex[kk]];
		}
	x = x - sum;

	double aij=getaij(row,j);
	x = x/aij;


	if (abs(colLow-colUpp) < tol)
		return; //here there is no restriction on zj so no bound on y

	if ((valuePrimal.at(j) - colLow) > tol && (colUpp - valuePrimal.at(j)) > tol) {
		//set both bounds
		if (x<*up)
			*up = x;
		if (x>*lo)
			*lo = x;
	}

	else if ((valuePrimal.at(j) == colLow && aij<0) || (valuePrimal.at(j) == colUpp && aij>0)) {
		if (x<*up)
			*up = x;

	}
	else if ((valuePrimal.at(j) == colLow && aij>0) || (valuePrimal.at(j) == colUpp && aij<0)) {
		if (x>*lo)
			*lo = x;
	}

}



/**
 * returns z_col
 * z = A'y + c
 */
double HPresolve::getColumnDualPost(int col) {
	int row;
	double z;
	double sum = 0;
	for (int cnt=Astart[col]; cnt<Aend[col]; cnt++)
		if (flagRow[Aindex[cnt]]) {
			row = Aindex[cnt];
			sum = sum + valueRowDual[row]*Avalue[cnt];
		}
	z = sum + colCostAtEl[col] ;
	return z;
}


/***
 * A'y + c = z
 *
 * returns y_row = -(A'y      +   c   - z )/a_rowcol
 *               (except row)  (at el)
 */
double HPresolve::getRowDualPost(int row, int col) {
	double x = 0;

	for (int kk=Astart[col]; kk<Aend[col];kk++)
		if (flagRow[Aindex[kk]] && Aindex[kk] != row)
			x = x + Avalue[kk]*valueRowDual[Aindex[kk]];

	x = x + colCostAtEl[col] - valueColDual[col];

	double y=getaij(row,col);
	return -x/y;
}




string HPresolve::getDualsForcingRow( int row, vector<int>& fRjs) {
	double x,aij,z;
	stringstream ss;
	int j;

	double lo = -HSOL_CONST_INF;
	double up = HSOL_CONST_INF;

	double cost, sum;

	for (int jj=0;jj<fRjs.size();jj++) {
			j = fRjs[jj];

			pair<int, vector<double>> p = oldBounds.top();
			vector<double> v = get<1>(p);
			oldBounds.pop();
			double colLow = v[0];
			double colUpp = v[1];


			//calculate bound x imposed by zj
			getBoundOnLByZj(row, j, &lo, &up, colLow, colUpp);
		}

	//calculate yi
	if (lo>up)
		cout<<"PR: Error in postsolving forcing row "<<row <<" : inconsistent bounds for it's dual value.\n";

	if (lo<=0 && up>=0) {
		valueRowDual[row] = 0;
	}
	else if (lo>0) {
		valueRowDual[row] = lo;
	}
	else if (up<0) {
		valueRowDual[row] = up;
	}

	flagRow[row] = 1;


	for (int jj=0;jj<fRjs.size();jj++) {
			j = fRjs[jj];
			cost = valueColDual.at(j);
			sum = 0;
			for (int k=Astart.at(j); k<Aend.at(j); k++)
				if (flagRow[Aindex.at(k)]) {
					sum = sum + valueRowDual[Aindex.at(k)]*Avalue.at(k);
					//cout<<" row "<<Aindex.at(k)<<" dual "<<valueRowDual[Aindex.at(k)]<<" a_"<<Aindex.at(k)<<"_"<<j<<"\n";
				}
			z = cost + sum;

			valueColDual.at(j) = z;

			if (iKKTcheck == 1) {
				ss<<j;
				ss<<" ";
				chk.addChange(2, 0, j, valuePrimal.at(j), valueColDual.at(j), cost);
			}
		}

	return ss.str();
}

void HPresolve::getDualsSingletonRow( int row, int col ) {

	pair< int ,vector<double>> bnd = oldBounds.top(); oldBounds.pop();

	valueRowDual[row] = 0;
	double cost = postValue.top(); postValue.pop();
	double aij = getaij(row, col);
	double l = (get<1>(bnd))[0];
	double u = (get<1>(bnd))[1];
	double lrow = (get<1>(bnd))[2];
	double urow = (get<1>(bnd))[3];
	if ((aij*valuePrimal[col] - lrow) > tol && ( -aij*valuePrimal[col] + urow) > tol) {
		valueRowDual[row] = 0;
		//row is nonbasic

	}
	else {
		if  ((valuePrimal[col] > l && valuePrimal[col] < u && abs(valueColDual[col]) > tol ) ||
				(valuePrimal[col] == l && valuePrimal[col] < u && valueColDual[col] < -tol ) ||
				(valuePrimal[col] == u && valuePrimal[col] > l && valueColDual[col] > tol )) {

			valueColDual[col] = 0;
		}

		double sum = 0;
		for (int k=Astart[col]; k<Aend[col]; k++)
			if (flagRow[Aindex.at(k)])
				sum = sum + valueRowDual[Aindex.at(k)]*Avalue.at(k);

		flagRow[row] = 1;


		double y = (valueColDual[col] - cost - sum)/aij;
		if (y != 0) {
			if (iKKTcheck == 1)
				chk.addCost(col, cost);

			//aij * yi + sum + ci = zi
			if (urow != lrow)
				if  ((aij*valuePrimal[col] == lrow && y >  tol) ||
					(aij*valuePrimal[col] == urow && y < -tol)) {


					//bounds on y_row
					double loY = -HSOL_CONST_INF;
					double upY = HSOL_CONST_INF;

					if (y>tol)
						upY = 0;
					else //y < -tol
						loY = 0;

					//bounds on z_col
					double loZ = -HSOL_CONST_INF;
					double upZ = HSOL_CONST_INF;
					if  (valuePrimal[col] == l && l < u ) {
						loZ = 0;
					}
					else if (valuePrimal[col] == u && l < u) {
						upZ = 0;
					}
					else if (l==u) {
						loZ = 0;
						upZ = 0;
					}
					//aij * yi + sum + ci = zi

					double lo = -HSOL_CONST_INF;
					double up = HSOL_CONST_INF;
					//bounds on z by y
					if (aij > 0) {
						if (loY >  -HSOL_CONST_INF) {
							lo = sum + cost + aij * loY;
							if (lo>loZ)
								loZ = lo;
						}
						if (upY < HSOL_CONST_INF) {
							up = sum + cost + aij * upY;
							if (up<upZ)
								upZ = up;
						}
					}
					else if (aij < 0) {
						if (loY >  -HSOL_CONST_INF) {
							up = sum + cost + aij * loY;
							if (up<upZ)
								upZ = up;
						}
						if (upY < HSOL_CONST_INF) {
							lo = sum + cost + aij * upY;
							if (lo>loZ)
								loZ = lo;
						}
					}
					//bounds on y by z
					//aij * yi  = zi - sum - ci
					if (aij > 0) {
						if (loZ >  -HSOL_CONST_INF) {
							lo = (loZ - sum - cost )/aij;
							if (lo>loY)
								loY = lo;
						}
						if (upZ < HSOL_CONST_INF) {
							up = (upZ - sum - cost )/aij;
							if (up<upY)
								upY = up;
						}
					}
					else if (aij < 0) {
						if (loZ >  -HSOL_CONST_INF) {
							up = (loZ - sum - cost )/aij;
							if (up<upY)
								upY = up;
						}
						if (upZ < HSOL_CONST_INF) {
							lo = (loZ - sum - cost )/aij;
							if (lo>loY)
								loY = lo;
						}
					}

					cout<< "loY= "<< loY <<" upY= "<<upY<<"loZ= "<< loZ <<" upZ= "<<upZ<<endl;

					valueRowDual[row] = 0;

					sum = 0;
					for (int k=Astart[col]; k<Aend[col]; k++)
						if (flagRow[Aindex.at(k)]) {
							sum = sum + valueRowDual[Aindex.at(k)]*Avalue.at(k);
							//cout<<" row "<<Aindex.at(k)<<" dual "<<valueRowDual[Aindex.at(k)]<<" a_"<<Aindex.at(k)<<"_"<<j<<"\n";
						}
					double newz = cost + sum;
					if ((valueColDual[col] > 0 && newz < 0) ||
						(valueColDual[col] < 0 && newz > 0)) {
						//valueColDual[col] = 0;
						//update valueRowDual[row]
						//newz = 0 if cost + sum + aijyi = 0 so aijyi = - cost - sum

						valueRowDual[row] = (-cost - sum)/aij;
						valueColDual[col] = 0;
						if (iKKTcheck == 1) {
							chk.addChange(2, 0, col, valuePrimal[col], valueColDual[col], cost);
						}
						return;

					}

					valueColDual[col] = newz;
					return;
				}
		valueRowDual[row] = y;
		}
	}

	flagRow[row] = 1;
	//row is introduced so something needs to become basic :



	//check if x is at a bound forced by the singleton row: then x becomes basic and row nonbasic
	if (nonbasicFlag[col] == 1) {
		// x was not basic but is now
		if (valuePrimal[col]  != l && valuePrimal[col] != u ) {
			nonbasicFlag[col] = 0;
			nonbasicFlag[numColOriginal + row] = 1;
		}
		// x was not basic and is not now either, row is basic
		else
			nonbasicFlag[numColOriginal + row] = 0;
	}
	else if (nonbasicFlag[col] == 0)  // x is basic
		nonbasicFlag[numColOriginal + row] = 0; //row becomes basic too

}




void HPresolve::getDualsDoubletonEquation(int row, int col) {
	// colDual already set. need valuePrimal from stack. maybe change rowDual depending on bounds. old bounds kept in oldBounds.
	// variables j,k : we eliminated col(k)(c.col) and are left with changed bounds on j and no row.
	//                               y                                                 x

	pair< int ,vector<double>> p = oldBounds.top(); oldBounds.pop();
	vector<double> v = get<1>(p);
	int x = get<0>(p);
	double ubxNew = v[1];
	double lbxNew = v[0];
	double cxNew = v[2];
	p = oldBounds.top(); oldBounds.pop();
	v = get<1>(p);
	double ubxOld = v[1];
	double lbxOld = v[0];
	double cxOld = v[2];
	p = oldBounds.top(); oldBounds.pop();
	v = get<1>(p);
	double uby = v[1];
	double lby = v[0];
	double cy = v[2];

	int y = col;

	double b = postValue.top(); postValue.pop();
	double aky = postValue.top(); postValue.pop();
	double akx = postValue.top(); postValue.pop();
	double valueX = valuePrimal[x];

	// primal value and objective shift
	valuePrimal[y] = (b - akx*valueX)/aky;
	objShift += -cxNew*valueX + cxOld*valueX + cy*valuePrimal[y];

	//column cost of x
	colCostAtEl[x] = cxOld;


	//get nzCol[y] before unflag row as missing
	int nzy = Aend[y] - Astart[y];
	for (int kk=Astart[y]; kk<Aend[y]; kk++)
		if (!flagRow[Aindex[kk]])
			nzy--;


	double lo = -HSOL_CONST_INF;
	double up = HSOL_CONST_INF;

	getBoundOnLByZj(row, x, &lo, &up, lbxOld, ubxOld);
	getBoundOnLByZj(row, y, &lo, &up, lby,    uby);

		//calculate yi
	if (lo-up > tol)
		cout<<"PR: Error in postsolving doubleton equation "<<row <<" : inconsistent bounds for it's dual value.\n";

	if (lo<=0 && up>=0) {
		valueRowDual[row] = 0;
	}
	else if (lo>0) {
		valueRowDual[row] = lo;
	}
	else if (up<0) {
		valueRowDual[row] = up;
	}

	flagRow[row] = 1;
	valueColDual[y] = getColumnDualPost(y);
	valueColDual[x] = getColumnDualPost(x);

	if (iKKTcheck == 1)
		chk.colDual[x] = valueColDual[x];

	if ((nonbasicFlag[x] == 1 && valueX==ubxNew && ubxNew < ubxOld)  ||
		(nonbasicFlag[x] == 1 && valueX==lbxNew && lbxNew > lbxOld))   {
			nonbasicFlag[x] = 0;
	}
	else {
		//row becomes basic unless y is between bounds, in which case y is basic
		if (valuePrimal[y] - lby > tol  && uby - valuePrimal[y] > tol) {
			nonbasicFlag[y] = 0;
		}
		else if (abs(valueX-ubxNew ) < tol || abs(valueX-lbxNew ) < tol)
			nonbasicFlag[y] = 0;
		else
			nonbasicFlag[numColOriginal + row] = 0;
	}

	//flagRow[row] = true;
	flagCol[y] = 1;
}

/*
void HPresolve::updateRemovedColRow(int dim) {
	list<int>::iterator it;

	if (dim==1) {
		it = singRow.begin();
		while (it != singRow.end())
			if (!flagRow[*it]) {
					it = singRow.erase(it);
					}
			else
				it++;
				}
	else if (dim==2) {
		it = singCol.begin();
		while (it != singCol.end())
			if (!flagCol[*it]) {
					it = singCol.erase(it);
					}
			else
				it++;
	}
}

void HPresolve::updateRowsByNZ() {
	//for (int j=Astart)


	for (int i = 0; i < numRow; i++)
		if (flagRow.at(i)) {
			if (nzRow.at(i) == 1) {
				int singletonColIndex = testSingRows(i);
				if (singletonColIndex > 0)
					singRow.push_back(i);
			}

			if (nzRow.at(i) == 0)
				if (b.at(i) == 0) {
					if (iPrint > 0)
						cout<<"PR: Empty row "<<i<<" removed. "<<endl;
					flagRow.at(i) = false;
					addChange(0, i, 0);
				}
				else {
					if (iPrint > 0)
						cout<<"PR: Problem infeasible."<<endl;
					cout<<("NOT-OPT status = 1, detected on presolve.\n");
					exit(1);
				}
		}

}


bool HPresolve::checkDuplicateColumns(int i, int k) {
	int j, p;
	//check all entries in i are either in k too or singleton
	for (p=Astart.at(i); p<Aend.at(i); p++) {
		j = Aindex[p];
		if (!flagRow.at(j))
			continue;
		else if (isZeroA(j, k))
			return false;
	}
	//the other way about
	for (p=Astart.at(k); p<Aend.at(k); p++) {
		j = Aindex[p];
		if (!flagRow.at(j))
			continue;
		else if (isZeroA(j, i))
			return false;
	}

	return true;
}


bool HPresolve::checkDuplicateRows(int i, int k) {
	int j, p;
	//check all entries in i are either in k too or singleton
	for (p=ARstart.at(i); p<ARstart[i+1]; p++) {
		j = ARindex[p];
		if (nzCol.at(j)==1 || !flagCol.at(j))
			continue;
		else if (isZeroA(k, j))
			return false;
	}
	//the other way about
	for (p=ARstart.at(k); p<ARstart[k+1]; p++) {
		j = ARindex[p];
		if (nzCol.at(j)==1 || !flagCol.at(j))
			continue;
		else if (isZeroA(i, j))
			return false;
	}

	return true;
}


void HPresolve::findDuplicateRows() {
	timer.recordStart(HTICK_PRE_DUPLICATE_ROWS);
	int v = 0; //remaining potential duplicates
	int t = 1; //number of potential next set to be created

	//row pass
	vector<int>    s(numRow,-1);
	for (int i=0;i<numRow;i++)
		if (flagRow.at(i)){
			v++;
			s.at(i) = 0;
			}

	//matrix pass
	int n,t0, ind, r,r0;
	for (int j=0; j<numCol; j++)
		if (flagCol.at(j) && nzCol.at(j)>1) {
			n = 0;
			t0 = t;
			t++;
			for (int ind=Astart.at(j); ind<Aend.at(j); ind++)  {
				r = Aindex[ind];
				if (flagRow[r]) {
					if (s[r]==0) {
						r0 = r;
						s[r] = t0;
						n++;
					}
					else if (s[r]>0) {
						bool isSat1 = false;
						for (int ii=Astart.at(j); ii<Aend.at(j); ii++)  {
							int i = Aindex[ii];
							if (flagRow.at(i) && i!=r) {
								if (s.at(i) == s[r]) {
									s.at(i) = t;
									//s[r] = t;//
									isSat1 = true;
								}
							}
						}
						if (isSat1) {
							s[r] = t;
							t++;
						}
						else
							s[r] = -1;
					}
				}
			}
			if (n==1) {
				s[r0] = -1;
				v--;
				if (v==0) {
					break; break; }
			}
		}

	queue<double> pairWise;
	int col, indi, indj;
	//if s.at(i) == s.at(j) rows i and j have the same sparsity pattern in nonsingleton columns
	for (int i=0;i<numRow;i++) {
		if (s.at(i)>-1) {
			for (int j=i+1;j<numRow;j++)
				if (s.at(i) == s.at(j)) {
					bool same = true;
					double ratio = 0;
					//check coefficients
					indi = ARstart.at(i);
					while (indi<ARstart[i+1] && same) {
						col = ARindex[indi];
						if (!flagCol[col] || nzCol[col] <= 1) 	{
							indi++; continue; }
						// we've reached a nonsingleton column in row i, got to find the same in j
						indj = ARstart.at(j);
						while (ARindex[indj]!=col)
							indj++;
						if (indj >= ARstart[j+1]) {
							cout<<"Error in part 2 of finding duplicate rows: rows "<<i<<" and "<<j<<".\n";
							exit(18);
						}
						if (ratio == 0)
							ratio = ARvalue[indi]/ARvalue[indj];
						else if (ratio == ARvalue[indi]/ARvalue[indj]) {
							indi++; continue; }
						else {
							same = false;
							break;
						}
						indi++;
					}
					if (same && ratio!=0) {
						pairWise.push(i);
						pairWise.push(j);
						pairWise.push(ratio);
					}
				}
		}
	}

	while (!pairWise.empty()) {
			int r1 = (int) pairWise.front(); pairWise.pop();
			int r2 = (int) pairWise.front(); pairWise.pop();
			double rat = pairWise.front(); pairWise.pop();
			//if (iPrint > 0) cout<<"PR: Duplicate rows detected: "<<r1<<" and "<<r2<<". ratio = "<<rat<<"."<<endl;

			if (!checkDuplicateRows(r1,r2))
				cout<<"ROWS NOT DUPLICATE: "<<r1<<" and "<<r2<<endl;
			removeDuplicateRows(r1, r2, rat);
		}
	timer.recordFinish(HTICK_PRE_DUPLICATE_ROWS);
}


void HPresolve::removeDuplicateRows(int i, int k, double v) {
	if (!flagRow.at(i) || !flagRow.at(k))
		return;
	vector<double> coeffK, coeffI;
	vector<int> colK, colI;
	for (int l = ARstart.at(i); l<ARstart[i+1];l++)
		if (flagCol[ARindex[l]] && nzCol[ARindex[l]]==1) {
			coeffK.push_back(-(1/v)*ARvalue[l]);
			colK.push_back(ARindex[l]);
			coeffI.push_back(ARvalue[l]);
			colI.push_back(ARindex[l]);
		}
	for (int l = ARstart.at(k); l<ARstart[k+1];l++)
		if (flagCol[ARindex[l]] && nzCol[ARindex[l]]==1) {
			coeffK.push_back(ARvalue[l]);
			colK.push_back(ARindex[l]);
			coeffI.push_back(-v*ARvalue[l]);
			colI.push_back(ARindex[l]);
		}

	// coeffK = (row k)-v*(row i)
	int j = makeCheckForDuplicateRows(k, i, coeffK, colK, 1/v, 1);
	if (j>0) { //fill stack with linear transformation of A for postsolve
		postValue.push((double) k);
		postValue.push((double) i);
		postValue.push(1/v);
		cnts[9]++;
		}
	if (j==0) {
		// coeffI = (row i)-v*(row k)
		j = makeCheckForDuplicateRows(i, k, coeffI, colI, v, 2);
		if (j>0) { //fill stack with linear transformation of A for postsolve
			postValue.push((double) i);
			postValue.push((double) k);
			postValue.push(v);
			cnts[9]++;
		}
	}
}

/**
* return value: 0 no action, 1 empty row, 2 doubleton eq, 3 implied free
**/
/*
int HPresolve::makeCheckForDuplicateRows(int k, int i, vector<double>& coeff, vector<int>& colIndex, double v, int whichIsFirst) {
	double newRowKLower, newRowKUpper;
	int ii,kk;
	if (whichIsFirst==1) {
		ii = i;
		kk = k;
	}
	else {
		ii = k;
		kk = i;
	}

	if (v>0) {
		if (rowUpper[ii] == HSOL_CONST_INF || rowLower[kk] == -HSOL_CONST_INF)
			newRowKLower = -HSOL_CONST_INF;
		else
			newRowKLower = max(rowLower[kk] - v*rowUpper[ii], -HSOL_CONST_INF);
		if (rowUpper[kk] == HSOL_CONST_INF || rowLower[ii] == -HSOL_CONST_INF)
			newRowKUpper = HSOL_CONST_INF;
		else
			newRowKUpper = min(rowUpper[kk] - v*rowLower[ii], HSOL_CONST_INF);
	}
	else {
		if (rowLower[ii] == -HSOL_CONST_INF || rowLower[kk] == -HSOL_CONST_INF)
			newRowKLower = -HSOL_CONST_INF;
		else
			newRowKLower = max(rowLower[kk] - v*rowLower[ii], -HSOL_CONST_INF);
		if (rowUpper[kk] == HSOL_CONST_INF || rowUpper[ii] == HSOL_CONST_INF)
			newRowKUpper = HSOL_CONST_INF;
		else
			newRowKUpper = min(rowUpper[kk] - v*rowUpper[ii], HSOL_CONST_INF);
	}

	//empty
	if (coeff.size()==0) { return 0;
		if (newRowKLower <= tol && newRowKUpper >= -tol) {
			if (iPrint > 0)
				cout<<"PR: Row "<<ii<<" added to duplicate row "<<kk<<" resulted in empty row. Row "<<kk<<" removed."<<endl;
			//update bounds on row ii
			//add old bounds OF row ii to checker and for postsolve
			if (iKKTcheck == 1) {
				vector<pair<int, double>> bndsL, bndsU;
				bndsL.push_back( make_pair( ii, rowLower[ii]));
				bndsU.push_back( make_pair( ii, rowUpper[ii]));
				chk.rLowers.push(bndsL);
				chk.rUppers.push(bndsU);
			}
			vector<double> bnds;
			bnds.push_back(rowLower[ii]);
			bnds.push_back(rowUpper[ii]);


			if (rowUpper[kk] < HSOL_CONST_INF) {
				double ubOfVi;
				if ((rowUpper[ii]==HSOL_CONST_INF && v > 0) || (rowLower[ii]==-HSOL_CONST_INF && v < 0))
					ubOfVi=HSOL_CONST_INF;
				else if (v>0)
					ubOfVi = v*rowUpper[ii];
				else
					ubOfVi = v*rowLower[ii];
				if (v>0)
					rowUpper[ii] = min(ubOfVi, rowUpper[kk])/v;
				else
					rowLower[ii] = min(ubOfVi, rowUpper[kk])/v;
			}

			if (rowUpper[kk] < HSOL_CONST_INF) {
				double lbOfVi;
				if ((rowLower[ii]==HSOL_CONST_INF && v > 0) || (rowUpper[ii]==-HSOL_CONST_INF && v < 0))
					lbOfVi=-HSOL_CONST_INF;
				if (v>0)
					lbOfVi = v*rowLower[ii];
				else
					lbOfVi = v*rowUpper[ii];
				if (v>0)
					rowLower[ii] = max(lbOfVi, rowLower[kk])/v;
				else
					rowUpper[ii] = max(lbOfVi, rowLower[kk])/v;
			}

			bnds.push_back(rowLower[ii]);
			bnds.push_back(rowUpper[ii]);
			oldBounds.push(make_pair( ii, bnds));

			flagRow[kk] = false;
			addChange(11, kk, 0);
		}
		else {
			if (iPrint > 0)
				cout<<"PR: Problem infeasible."<<endl;
			cout<<("NOT-OPT status = 1, detected on presolve.\n");
			exit(1);
		}

		return 1;
	}
	/* /doubleton eq with a singleton equation
	else if (coeff.size()==2) { return 0;
		int col, j;
		double aij, aicol;

		if (isZeroA(i, colIndex[0])) {
			col = colIndex[0];
			j   = colIndex[1];
			aij = coeff[1];
			aicol = coeff[0];
		}
		else if (isZeroA(i, colIndex[1])) {
			j = colIndex[0];
			col   = colIndex[1];
			aij = coeff[0];
			aicol = coeff[1];
		}
		else return 0;
		cout<<"doubleton free"<<endl;
   		//modify bounds on variable j
		//set row 's bounds to the new ones before we eliminate it
		rowLower[kk] = newRowKLower;
		rowUpper[kk] = newRowKUpper;

		// additional check if it is indeed implied free: need
		// low and upp to be tighter than original bounds for variable col
		// so it is indeed implied free and we can remove it
		pair<double, double> p = getNewBoundsDoubletonConstraint(kk, j, col, aij, aicol);
		double low, upp;
		low = get<0>(p);
		upp = get<1>(p);
		if (!(colLower[col] <= low && colUpper[col] >= upp))
			return 0;


		postValue.push(aij);
		postValue.push(aicol);
		postValue.push(newRowKUpper);
		postValue.push(newRowKLower);

		//modify bounds on variable j, variable col (k) is substituted out
		//double aik = Avalue.at(k);
		//double aij = Avalue[kk];
		p = getNewBoundsDoubletonConstraint(kk, col, j, aicol, aij);
		low = get<0>(p);
		upp = get<1>(p);

		//add old bounds of xj to checker and for postsolve
		if (iKKTcheck == 1) {
			vector<pair<int, double>> bndsL, bndsU, costS;
			bndsL.push_back( make_pair( j, colLower.at(j)));
			bndsU.push_back( make_pair( j, colUpper.at(j)));
			costS.push_back( make_pair( j, colCost.at(j)));
			chk.cLowers.push(bndsL);
			chk.cUppers.push(bndsU);
			chk.costs.push(costS);
		}

		vector<double> bnds;
		bnds.push_back(colLower[col]);
		bnds.push_back(colUpper[col]);
		bnds.push_back(colCost[col]);
		oldBounds.push(make_pair( col, bnds));
		bnds.clear();
		bnds.push_back(colLower.at(j));
		bnds.push_back(colUpper.at(j));
		bnds.push_back(colCost.at(j));
		oldBounds.push(make_pair( j, bnds));

		if (low > colLower.at(j))
			colLower.at(j) = low;
		if (upp < colUpper.at(j))
			colUpper.at(j) = upp;

		//modify cost of xj
		colCost.at(j) = colCost.at(j) - colCost[col]*aij/aicol;

		//for postsolve: need the new bounds too
		//oldBounds.push_back(colLower.at(j)); oldBounds.push_back(colUpper.at(j));
		bnds.clear();
		bnds.push_back(colLower.at(j));
		bnds.push_back(colUpper.at(j));
		bnds.push_back(colCost.at(j));
		oldBounds.push(make_pair( j, bnds));

		if (iPrint > 0)
			cout<<"PR: Row "<<ii<<" added to duplicate row "<<kk<<" resulted in a doubleton equation. Variable "<<col<<" and row "<<kk<<" removed. variable left is "<<j<<endl;
		flagCol[col] = false;

		valueColDual[col] = 0;
		valueRowDual[kk] = -colCost[col]/aicol; //may be changed later, depending on bounds.



		removeRow(kk);
		addChange(12, kk, col);
		return 2;
	}
	//implied free
	else if (coeff.size()>2){
		for (int j=0;j<colIndex.size(); j++) {
			if (!flagRow.at(k))
				return 0;

			//check if singleton
			if (whichIsFirst == 1 && !isZeroA(i, colIndex.at(j)))
			 continue;

			if (whichIsFirst == 2 && !isZeroA(k, colIndex.at(j)))
			 continue;

			int elem = getaij(k, colIndex.at(j));
			bool res = removeIfImpliedFree(colIndex.at(j), k, elem);
			if (res) {
				cout<<"impl free"<<endl;
				addChange(11, k, colIndex.at(j));
				return 3;
			}
		}
	}* /
	return 0;
}
/*
void HPresolve::findDuplicateColumns() {
	int v = 0; //remaining potential duplicates
	int t = 1; //number of potential next set to be created

	//column pass
	vector<int> s(numCol,-1);
	for (int i=0;i<numCol;i++)
		if (flagCol.at(i)){
			v++;
			s.at(i) = 0;
			}

	//matrix pass
	int n,t0, ind, r,r0;
	for (int j=0; j<numRow; j++)
		if (flagRow.at(j)) { //row by row
			n = 0;
			t0 = t;
			t++;
			for (int ind=ARstart.at(j); ind<ARstart[j+1]; ind++)  {
				r = ARindex[ind];
				if (flagCol[r]) {
					if (s[r]==0) {
						r0 = r;
						s[r] = t0;
						n++;
					}
					else if (s[r]>0) {
						bool isSat1 = false;
						for (int ii=ARstart.at(j); ii<ARstart[j+1]; ii++)  {
							int i = ARindex[ii];
							if (flagCol.at(i) && i!=r) {
								if (s.at(i) == s[r]) {
									s.at(i) = t;
									s[r] = t;
									isSat1 = true;
								}
							}
						}
						if (isSat1)
							t++;
						else
							s[r] = -1;
					}
				}
			}
			if (n==1) {
				s[r0] = -1;
				v--;
				if (v==0) {
					break; break;	}
			}
		}

		queue<double> pairWise;
		int row, indi, indj;

		//if s.at(i) == s.at(j) columns i and j have the same sparsity pattern
		for (int i=0;i<numCol;i++) {
			if (s.at(i)>-1) {
				for (int j=i+1;j<numCol;j++)
					if (s.at(i) == s.at(j)) {
						bool same = true;
						double ratio = 0;
						//check coefficients
						indi = Astart.at(i);
						indj = Astart.at(j);
						while (indi<Aend.at(i) && same) {
							row = Aindex[indi];
							if (!flagRow[row]) 	{
								indi++; continue; }
							while (Aindex[indj]!=row)
								indj++;
							if (ratio == 0)
								ratio = Avalue[indi]/Avalue[indj];
							else if (ratio == Avalue[indi]/Avalue[indj]) {
								indi++; continue; }
							else {
								same = false;
								break;
							}
							indi++;
						}
						if (same) {
							pairWise.push(i);
							pairWise.push(j);
							pairWise.push(ratio);
						}
					}
			}
		}


		while (!pairWise.empty()) {
			int c1 = (int) pairWise.front(); pairWise.pop();
			int c2 = (int) pairWise.front(); pairWise.pop();
			double rat = pairWise.front(); pairWise.pop();
			//if (iPrint > 0) cout<<"PR: Duplicate columns detected: "<<c1<<" and "<<c2<<". ratio = "<<rat<<"."<<endl;
			if (checkDuplicateColumns(c1, c2))
				removeDuplicateColumns(c1,c2,rat);
		}


}

/*
void HPresolve::removeDuplicateColumns(int j, int k, double v) {
	//in case there's more than two duplicate columns
	if (!flagCol.at(j) || !flagCol.at(k))
		return;
	double val = colCost.at(j) - v*colCost.at(k);
	double xj;
	//fix a duplicate column
	if (val!=0) {
		//has only a lower bound
		if (colLower.at(k)>-HSOL_CONST_INF && colUpper.at(k) == HSOL_CONST_INF) {
			if (v>=0 && val>0) {
				//zj > 0, variable j is at lower bound.
				if (colLower.at(j) == -HSOL_CONST_INF) {
					if (iPrint > 0)
						cout<<"PR: Problem dual infeasible."<<endl;
					cout<<"NOT-OPT status = 2, detected on presolve.\n";
					exit(2);
				}
		        if (iKKTcheck == 1)   {
			        chk.rLowers.push(rowLower);
			        chk.rUppers.push(rowUpper);
		        }
				xj = colLower.at(j);
				setPrimalValue(j, xj);
				addChange(7, 0, j);//postsolve is like a fixed variable
				if (iPrint > 0)
					cout<<"PR: Duplicate column "<<j<<" fixed. Value := "<<xj<<endl;
			}
			else if (v<=0 && val<0) {
				//zj < 0, variable j is at upper bound.
				if (colUpper.at(j) == HSOL_CONST_INF) {
					if (iPrint > 0)
						cout<<"PR: Problem dual infeasible."<<endl;
					cout<<"NOT-OPT status = 2, detected on presolve.\n";
					exit(2);
				}
		        if (iKKTcheck == 1)   {
			        chk.rLowers.push(rowLower);
			        chk.rUppers.push(rowUpper);
		        }
				xj = colUpper.at(j);
				setPrimalValue(j, xj);
				addChange(7, 0, j);//postsolve is like a fixed variable
				if (iPrint > 0)
					cout<<"PR: Duplicate column "<<j<<" fixed. Value := "<<xj<<endl;
			}
		}
		//has only an upper bound
		if (colLower.at(k)==-HSOL_CONST_INF && colUpper.at(k) < HSOL_CONST_INF) {
			if (v<=0 && val>0) {
				//zj > 0, variable j is at lower bound.
				if (colLower.at(j) == -HSOL_CONST_INF) {
					if (iPrint > 0)
						cout<<"PR: Problem dual infeasible."<<endl;
					cout<<"NOT-OPT status = 2, detected on presolve.\n";
					exit(2);
				}
		        if (iKKTcheck == 1)   {
			        chk.rLowers.push(rowLower);
			        chk.rUppers.push(rowUpper);
		        }
				xj = colLower.at(j);
				setPrimalValue(j, xj);
				addChange(7, 0, j);//postsolve is like a fixed variable
				if (iPrint > 0)
					cout<<"PR: Duplicate column "<<j<<" fixed. Value := "<<xj<<endl;
			}
			else if (v>=0 && val<0) {
				//zj < 0, variable j is at upper bound.
				if (colUpper.at(j) == HSOL_CONST_INF) {
					if (iPrint > 0)
						cout<<"PR: Problem dual infeasible."<<endl;
					cout<<"NOT-OPT status = 2, detected on presolve.\n";
					exit(2);
				}
				xj = colUpper.at(j);
		        if (iKKTcheck == 1)   {
			        chk.rLowers.push(rowLower);
			        chk.rUppers.push(rowUpper);
		        }
				setPrimalValue(j, xj);
				addChange(7, 0, j);//postsolve is like a fixed variable

			}
		}
	}
	//replace two by one:
	else {
		oldBounds.push_back(v);
		oldBounds.push_back((double) k);
		oldBounds.push_back(colLower.at(k));
		oldBounds.push_back(colUpper.at(k));

		if (iKKTcheck == 1) {
			chk.cLowers.push(colLower);
			chk.cUppers.push(colUpper);
		}

		if (iPrint > 0)
			cout<<"PR: Duplicate columns "<<k<<" and "<<j<<" replaced by one at "<<k<<"."<<endl;
		if (v<0) {
		//deal with infinities.
			if (colLower.at(k)>-HSOL_CONST_INF && colUpper.at(j)<HSOL_CONST_INF) {
				colLower.at(k) = colLower.at(k) + v*colUpper.at(j);
				if (colLower.at(k) < -HSOL_CONST_INF)
					colLower.at(k) = -HSOL_CONST_INF;
			}
			else
				colLower.at(k)	= -HSOL_CONST_INF;

			if (colUpper.at(k)<HSOL_CONST_INF && colLower.at(j)> -HSOL_CONST_INF) {
				colUpper.at(k) = colUpper.at(k) + v*colLower.at(j);
				if (colUpper.at(k) > HSOL_CONST_INF)
					colUpper.at(k) = HSOL_CONST_INF;
			}
			else
				colUpper.at(k) = HSOL_CONST_INF;
		}
		else if (v>0) {
			if (colLower.at(k)>-HSOL_CONST_INF && colLower.at(j)>-HSOL_CONST_INF) {
				colLower.at(k) = colLower.at(k) + v*colLower.at(j);
				if (colLower.at(k) < -HSOL_CONST_INF)
					colLower.at(k) = -HSOL_CONST_INF;
			}
			else
				colLower.at(k)	= -HSOL_CONST_INF;

			if (colUpper.at(k)<HSOL_CONST_INF && colUpper.at(j)> -HSOL_CONST_INF) {
				colUpper.at(k) = colUpper.at(k) + v*colUpper.at(j);
				if (colUpper.at(k) > HSOL_CONST_INF)
					colUpper.at(k) = HSOL_CONST_INF;
			}
			else
				colUpper.at(k) = HSOL_CONST_INF;
		}
		addChange(14, k, j);
		flagCol.at(j) = false;

		//update nonzeros: two by one so for each row -1
		for (int kk=Astart.at(j); kk< Aend.at(j); kk++) {
			if (flagRow[Aindex[kk]])
				nzRow[Aindex[kk]]--;
				}
	}
}*/


/*/testing and dev
void HPresolve::setNBFfullproblem(vector<int>& nbfFull, vector<int>& bnFull) {
	nbffull = nbfFull;
	bindfull = bnFull;
}


void HPresolve::cmpNBF(int row, int col) {
	int i;
	bool print = false;


	int countBasic=0, nR=0;

     for (i=0; i< nonbasicFlag.size();i++) {
    	 if (nonbasicFlag.at(i) == 0)
    			 countBasic++;
     }

     for (i=0;i<numRowOriginal;i++)
		if (flagRow.at(i)) {
			nR++;
			}

     if (countBasic != nR)
    	 cout<<" Wrong count of basic variables:  numRow="<<nR<<" basic="<< countBasic << "\n";

     if (col>-1 && row>-1) {
		 i = col;
		 if (nbffull.at(i) != nonbasicFlag.at(i))
			 cout<<"		nbffull ["<<i<<"] ="<< nbffull.at(i) <<"		nonbasicFlag ["<<i<<"] = "<<  nonbasicFlag.at(i) <<endl;
		 i = numColOriginal + row;
		 if (nbffull.at(i) != nonbasicFlag.at(i)) {
			 cout<<"		nbffull ["<<i<<"] ="<< nbffull.at(i) <<"		nonbasicFlag ["<<i<<"] = "<<  nonbasicFlag.at(i) <<endl;
		 }
	 }



	// basis different.
	/*
	//nbffull
	if (nbffull.size() != nonbasicFlag.size())
		cout<<"nbf count diff"<<endl;

	else {
		for (i=0; i<nbffull.size(); i++) {
			if (nbffull.at(i) != nonbasicFlag.at(i)) {
				if (i<numColOriginal) {
					if (flagCol.at(i)) {
						cout<<"		nbffull ["<<i<<"] ="<< nbffull.at(i) <<"		nonbasicFlag ["<<i<<"] = "<<  nonbasicFlag.at(i) <<endl;
						print = true;
					}
				}
				else {
					if (flagRow[i - numColOriginal]) {
						cout<<"		nbffull ["<<i<<"] ="<< nbffull.at(i) <<"		nonbasicFlag ["<<i<<"] = "<<  nonbasicFlag.at(i) <<endl;
						print = true;
					}
				}
			}
		}

		//print = true;
		if (false) {
			for (i=0; i<nbffull.size(); i++) {
				cout<<nbffull.at(i) <<" ";
			}
			cout<<endl;
			for (i=0; i<nonbasicFlag.size(); i++) {
				cout<<nonbasicFlag.at(i) <<" ";
			}
		}
	}

	if (!print) {
		cout<<"nbf SAME"<<endl;
	}
	cout<<endl; * /

	//bindfull
//	if (bindfull.size() != basicIndex.size())
//		cout<<"bindfull count diff"<<endl;
//
//	else
//		for (i=0; i<bindfull.size(); i++) {
//			if (nbffull.at(i) != basicIndex.at(i))
//				cout<<"		bindfull ["<<i<<"] diff"<<endl;
//		}
}*/
