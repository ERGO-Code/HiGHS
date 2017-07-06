/*
 * HinOut.h
 *
 *  Created on: 16 Jan 2017
 *      Author: s1131817
 */

#ifndef HINOUT_H_
#define HINOUT_H_

#include "HModel.h"
#include <cstring>

class HinOut : public HPreData {
public:
	int AcountX;
	string fileIn, fileOut;

	HinOut(string filenameIn, string filenameOut);
	void getData(HModel & ptr_model);
	void setData(HModel & ptr_model);
	void readDataPostsolve(HModel & ptr_model);
	void readDataColumnWise();
	void writeDataColumnWise();
	void clearData();
	void HinOutTestIO(HModel & ptr);
	void HinOutTestRead(HModel & ptr);
	void HinOutTestWrite(HModel & ptr);
	void compareData(int lvl);
	double getdiff(double v1,double v2);

	//data we are getting and printing
	int onumCol;
	int onumRow;
	int onumTot;

	vector<int> oAstart;
	vector<int> oAindex;
	vector<double> oAvalue;

	vector<double> ocolCost;
	vector<double> ocolLower;
	vector<double> ocolUpper;
	vector<double> orowLower;
	vector<double> orowUpper;

    int oAcountX;

};

#endif /* HINOUT_H_ */
