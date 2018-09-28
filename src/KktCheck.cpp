
#include "KktCheck.h"

#include <vector>
	

void KktCheck::printAR() {
	cout<<"N="<<numCol<<",  M="<<numRow<<",  NZ= "<<ARstart[numRow]<<'\n';	

	char buff [7];
	cout<<"\n-----cost-----\n";
	for (size_t i=0;i<colCost.size();i++) {
		sprintf(buff, "%2.1g ", colCost[i]);
		cout<<std::setw(7)<<buff; 
	}
	cout<<endl;
	cout<<"------AR | b----KktCheck-\n";
	for (i=0;i<numRow;i++) {
		for (j=0;j<numCol;j++) {

			int ind = ARstart[i];
			while (ARindex[ind]!=j && ind<ARstart[i+1])
				ind++;
			//if a_ij is nonzero print
			if (ARindex[ind]==j && ind<ARstart[i+1])
			{	
				sprintf(buff, "%2.1g ", ARvalue[ind]);
				cout<<std::setw(7)<<buff; 
				}
			else cout<<std::setw(7)<<"   ";
		
		}
		cout<<"  |   "<<std::setw(7)<<rowLower[i]<<" < < "<<rowUpper[i]<<endl;
	}
	cout<<endl;
	cout<<"------l------\n";
	for (int i=0;i<numCol;i++) { 
		if (colLower[i]>-HSOL_CONST_INF)
			sprintf(buff, "%2.1g ", colLower[i]);
		else 
			sprintf(buff, "-inf");
		cout<<setw(7)<<buff; 
	}
	cout<<endl;
	cout<<"------u------\n";
	for (int i=0;i<numCol;i++) { 
		if (colUpper[i]<HSOL_CONST_INF)
			sprintf(buff, "%2.1g ", colUpper[i]);
		else 
			sprintf(buff, "inf");
		cout<<setw(7)<<buff; 
	}
	cout<<endl;
}

void KktCheck::makeARCopy() {
	tol = 0.00001;
	// Make a AR copy
	vector<int> iwork(numRow, 0);
	ARstart.resize(numRow + 1, 0);
	int AcountX = Aindex.size();
	ARindex.resize(AcountX);
	ARvalue.resize(AcountX);
	for (k = 0; k < AcountX; k++)
	    iwork[Aindex[k]]++;
	for (i = 1; i <= numRow; i++)
	    ARstart[i] = ARstart[i - 1] + iwork[i - 1];
	for (i = 0; i < numRow; i++)
	    iwork[i] = ARstart[i];
	for (int iCol = 0; iCol < numCol; iCol++) {
	    for (k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
	        int iRow = Aindex[k];
	        int iPut = iwork[iRow]++;
	        ARindex[iPut] = iCol;
	        ARvalue[iPut] = Avalue[k];
	    }
	}
}

void KktCheck::chPrimalBounds() {
	for (i=0; i<numCol;i++) {
		if ((colLower[i] - colValue[i] > tol) ||  (colValue[i]- colUpper[i] > tol) ) {
			if (print == 1) 
				cout<<"Variable "<<cIndexRev[i]<<" infeasible: lb="<<colLower[i]<<", vaule="<<colValue[i]<<",  ub="<<colUpper[i]<<endl;
				//cout<<"Variable "<<i<<" infeasible: lb="<<colLower[i]<<", vaule="<<colValue[i]<<",  ub="<<colUpper[i]<<endl;
			istrueGlb = true;
		}
	}
}

void KktCheck::chPrimalFeas() {
	bool istrue = true;
	double rowV;
	//Ax = b
	for (i=0; i<numRow;i++) {
		rowV = 0;
		for (k=ARstart[i]; k<ARstart[i+1]; k++) 
			rowV = rowV + colValue[ARindex[k]]*ARvalue[k];

		if ( ((rowV-rowLower[i]) < 0) && (abs(rowV-rowLower[i]) > tol)) {
			if (print == 1) 
				cout<<"Row "<<rIndexRev[i]<<" infeasible: Row value="<<rowV<<"  L="<<rowLower[i]<<"  U="<<rowUpper[i]<<endl;
				//cout<<"Row "<<i<<" infeasible: Row value="<<rowV<<"  L="<<rowLower[i]<<"  U="<<rowUpper[i]<<endl;
			istrue = false;
		}

		if ( ((rowV-rowUpper[i]) > 0) && (abs(rowV-rowUpper[i]) > tol)) {
					if (print == 1)
						cout<<"Row "<<rIndexRev[i]<<" infeasible: Row value="<<rowV<<"  L="<<rowLower[i]<<"  U="<<rowUpper[i]<<endl;
						//cout<<"Row "<<i<<" infeasible: Row value="<<rowV<<"  L="<<rowLower[i]<<"  U="<<rowUpper[i]<<endl;
					istrue = false;
				}
	}
	
	if (istrue) {
		if (print == 1) 
			cout<<"Primal feasible.\n";
	}
	else {
		if (print == 1) 
			cout<<"KKT check error: Primal infeasible.\n";
		istrueGlb = true;
	}

}

void KktCheck::chDualFeas() {
	bool istrue = true;
	
	//check values of z_j are dual feasible
	for (i=0; i<numCol;i++) {
		// j not in L or U
		if (colLower[i] == -HSOL_CONST_INF && colUpper[i] == HSOL_CONST_INF) {
			if (abs(colDual[i]) > tol) {
				if (print == 1) 
					cout<<"Dual feasibility fail: l=-inf, x["<<cIndexRev[i]<<"]="<<colValue[i]<<", u=inf, z["<<i<<"]="<<colDual[i]<<endl;
					//cout<<"Dual feasibility fail: l=-inf, x["<<i<<"]="<<colValue[i]<<", u=inf, z["<<i<<"]="<<colDual[i]<<endl;
				istrue = false;
			}
		}
		// j in L: x=l and l<u
		else if (colValue[i] == colLower[i] && colLower[i] < colUpper[i]) {
			if (colDual[i]<0 && abs(colDual[i]) > tol) {
				if (print == 1)  
					cout<<"Dual feasibility fail: l["<<cIndexRev[i]<<"]="<<colLower[i]<<" = x["<<cIndexRev[i]<<"]="<<colValue[i]<<", z["<<cIndexRev[i]<<"]="<<colDual[i]<<endl;
					//cout<<"Dual feasibility fail: l["<<i<<"]="<<colLower[i]<<" = x["<<i<<"]="<<colValue[i]<<", z["<<i<<"]="<<colDual[i]<<endl;
				istrue = false;
			}
		}
		// j in U: x=u and l<u
		else if (colValue[i] == colUpper[i] && colLower[i] < colUpper[i]) {
			if (colDual[i]>tol) {
				if (print == 1)  
					cout<<"Dual feasibility fail: x["<<cIndexRev[i]<<"]="<<colValue[i]<<"=u["<<cIndexRev[i]<<"], z["<<cIndexRev[i]<<"]="<<colDual[i]<<endl;
					//cout<<"Dual feasibility fail: x["<<i<<"]="<<colValue[i]<<"=u["<<i<<"], z["<<i<<"]="<<colDual[i]<<endl;
				istrue = false;
			}
		}
	}
	
	//check values of y_i are dual feasible
	for (i=0; i<numRow;i++) {
		double rowV = 0;
		for (k=ARstart[i]; k<ARstart[i+1]; k++)
			rowV = rowV + colValue[ARindex[k]]*ARvalue[k];

		// L = Ax = U can be any sign
		if (abs(rowLower[i] - rowV) < tol && abs(rowUpper[i] - rowV) < tol) {
		}
		// L = Ax < U
		else if (abs(rowLower[i] - rowV) < tol && rowV < rowUpper[i]) {

			if (rowDual[i] > tol) {
				if (print == 1)
					cout<<"Dual feasibility fail for row "<<rIndexRev[i]<<": L= "<<rowLower[i] <<", Ax="<<rowV<<", U="<<rowUpper[i]<<", y="<<rowDual[i]<<endl;
					//cout<<"Dual feasibility fail for row "<<i<<": L= "<<rowLower[i] <<", Ax="<<rowV<<", U="<<rowUpper[i]<<", y="<<rowDual[i]<<endl;
				istrue = false;
			}
		}
		// L < Ax = U
		else if (rowLower[i] < rowV && abs(rowV - rowUpper[i]) < tol) {
			//cout<<"Dual feasibility fail for row "<<i<<": L= "<<rowLower[i] <<", Ax="<<rowV<<", U="<<rowUpper[i]<<", y="<<rowDual[i]<<endl;
			if (rowDual[i] < -tol) {
				if (print == 1)
					cout<<"Dual feasibility fail for row "<<rIndexRev[i]<<": L= "<<rowLower[i] <<", Ax="<<rowV<<", U="<<rowUpper[i]<<", y="<<rowDual[i]<<endl;
					//cout<<"Dual feasibility fail for row "<<i<<": L= "<<rowLower[i] <<", Ax="<<rowV<<", U="<<rowUpper[i]<<", y="<<rowDual[i]<<endl;
				istrue = false;
			}
		}
		// L < Ax < U
		else if ((rowLower[i] < (rowV+tol)) && (rowV < (rowUpper[i]+tol))) {
			if (abs(rowDual[i]) > tol) {
				if (print == 1)
					cout<<"Dual feasibility fail for row "<<rIndexRev[i]<<": L= "<<rowLower[i] <<", Ax="<<rowV<<", U="<<rowUpper[i]<<", y="<<rowDual[i]<<endl;
					//cout<<"Dual feasibility fail for row "<<i<<": L= "<<rowLower[i] <<", Ax="<<rowV<<", U="<<rowUpper[i]<<", y="<<rowDual[i]<<endl;istrue = false;
				istrue = false;
			}
		}
	}
	
	if (istrue) {
		if (print == 1) 
			cout<<"Dual feasible.\n";
	}
	else  {
		if (print == 1) 
			cout<<"KKT check error: Dual feasibility fail.\n";
		istrueGlb = true;
	}
	
}


void KktCheck::chComplementarySlackness() {
	bool istrue = true;
	
	for (i=0; i<numCol;i++) {
		if (colLower[i] > -HSOL_CONST_INF)  
			if (abs( (colValue[i]- colLower[i])*(colDual[i]) ) > tol && colValue[i]!=colUpper[i] && abs(colDual[i])> tol) {
				if (print == 1) 
					cout<<"Comp. slackness fail: "<<"l["<<cIndexRev[i]<<"]="<<colLower[i]<<", x["<<i<<"]="<<colValue[i]<<", z["<<i<<"]="<<colDual[i]<<endl;
					//cout<<"Comp. slackness fail: "<<"l["<<i<<"]="<<colLower[i]<<", x["<<i<<"]="<<colValue[i]<<", z["<<i<<"]="<<colDual[i]<<endl;
				istrue = false;
			}
		if (colUpper[i] < HSOL_CONST_INF) 
			if (abs( (colUpper[i] - colValue[i])*(colDual[i]) )> tol && colValue[i]!=colLower[i] && abs(colDual[i])> tol) {
				if (print == 1) 
					cout<<"Comp. slackness fail: x["<<cIndexRev[i]<<"]="<<colValue[i]<<", u["<<i<<"]="<<colUpper[i]<<", z["<<i<<"]="<<colDual[i]<<endl;
					//cout<<"Comp. slackness fail: x["<<i<<"]="<<colValue[i]<<", u["<<i<<"]="<<colUpper[i]<<", z["<<i<<"]="<<colDual[i]<<endl;
				istrue = false;
			}
	}
	
	if (istrue) {
		if (print == 1) 
			cout<<"Complementary Slackness.\n";}
	else {
		if (print == 1) 
			cout<<"KKT check error: Comp slackness fail.\n";
		istrueGlb = true;
	}
}




void KktCheck::printSol() {
	char buff [10];
	cout<<endl<<"Col value: ";
	for (size_t i=0;i<colValue.size();i++) {
		sprintf(buff, "%2.2f ", colValue[i]); 
		cout<<setw(5)<<buff; 
		}  
	cout<<endl<<"Col dual:  ";
	for (size_t i=0;i<colDual.size();i++) {
		sprintf(buff, "%2.2f ", colDual[i]);
		cout<<setw(5)<<buff; 
		}  
/*	cout<<endl<<"Row value: ";
	for (i=0;i<numRow;i++) {
	  	sprintf(buff, "%2.2f ", rowValue[i]);  
	  	cout<<setw(5)<<buff;
	  	}*/
	cout<<endl<<"Row dual:  ";
	for (size_t i=0;i<rowDual.size();i++) {
		sprintf(buff, "%2.2f ", rowDual[i]);  
		cout<<setw(5)<<buff; 
		}  
	cout<<endl<<endl;
}

void KktCheck::chStOfLagrangian() {
	bool istrue = true;
	double lagrV;
	//A'y + c - z = 0
	for (j=0; j<numCol;j++) {
		lagrV = colCost[j] -colDual[j];
		for (k=Astart[j]; k<Astart[j+1]; k++) 
				lagrV = lagrV + rowDual[Aindex[k]]*Avalue[k];
				
		if (abs(lagrV) > tol) {
			if (print == 1) 
				cout<<"Column "<<cIndexRev[j]<<" fails stationary of Lagrangian: dL/dx"<<j<<" = "<<lagrV<<", rather than zero."<<endl;
				//cout<<"Column "<<j<<" fails stationary of Lagrangian: dL/dx"<<j<<" = "<<lagrV<<", rather than zero."<<endl;
			istrue = false;
		}
	}
	
	if (istrue) {
		if (print == 1) 
			cout<<"Stationarity of Lagrangian.\n";
	}
	else {
		if (print == 1) 
			cout<<"KKT check error: Lagrangian is not stationary.\n";
		istrueGlb = true;
	}

}

void KktCheck::checkKKT() {
	if (numCol == 0)
		return;
	
	istrueGlb = false;
	

	makeARCopy();
	//printAR();printSol();
	chPrimalBounds() ;
	chPrimalFeas();
	chDualFeas();
	chComplementarySlackness();
	chStOfLagrangian();
	
	if (print == 2) {
		ofstream myfile;
  		myfile.open ("../experiments/out", ios::app ); 
		if (istrueGlb)	
			myfile << "           KKT fail      ";
		else
			myfile << "           KKT pass      ";
  		myfile.close();
	}
		
}



void KktCheck::passSolution(const vector<double>& colVal, const vector<double>& colDu,  const vector<double>& rDu) {
	colValue = colVal;
	colDual = colDu;
	rowDual = rDu;
}
//get DATA
void KktCheck::setMatrix(const vector<int>& Astart_, const  vector<int>& Aindex_, const  vector<double>& Avalue_) {
	Astart = Astart_;
	Aindex = Aindex_;
	Avalue = Avalue_;
}

void KktCheck::setBounds(const  vector<double>& colUpper_, const  vector<double>& colLower_) {
	colLower = colLower_;
	colUpper = colUpper_;
}

void KktCheck::setNumbersCostRHS(int nCol, int nRow, const vector<double>& rowLower_, const vector<double>& rowUpper_, const vector<double>& cost) {
	numCol = nCol;
	numRow = nRow;
	colCost = cost;
	rowLower = rowLower_;
	rowUpper = rowUpper_;
}

void KktCheck::setIndexVectors(vector<int>& rIndex, vector<int>& cIndex){
	rIndexRev.clear();
	cIndexRev.clear();

	for (size_t i=0;i<rIndex.size(); i++) {
		if (rIndex[i] != -1) {
			rIndexRev.push_back(i);
		}
	}
	for (size_t i=0;i<cIndex.size(); i++) {
		if (cIndex[i] != -1) {
			cIndexRev.push_back(i);
		}
	}

}

