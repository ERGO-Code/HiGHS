#ifndef IO_HMPSIO_H_
#define IO_HMPSIO_H_

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
using namespace std;

int readMPS(const char* filename, int mxNumRow, int mxNumCol, int& numRow,
            int& numCol, int& objSense, double& objOffset, vector<int>& Astart,
            vector<int>& Aindex, vector<double>& Avalue,
            vector<double>& colCost, vector<double>& colLower,
            vector<double>& colUpper, vector<double>& rowLower,
            vector<double>& rowUpper, vector<int>& integerColumn);

int writeMPS(const char* filename, int& numRow, int& numCol, int& numInt,
             int& objSense, double& objOffset, vector<int>& Astart,
             vector<int>& Aindex, vector<double>& Avalue,
             vector<double>& colCost, vector<double>& colLower,
             vector<double>& colUpper, vector<double>& rowLower,
             vector<double>& rowUpper, vector<int>& integerColumn);

bool load_mpsLine(FILE* file, int& integerVar, int lmax, char* line, char* flag,
                  double* data);
bool hsol_isInfinity(double val);
int isspace(int c);
char* fgets(char* str, int num, FILE* stream);

inline const char* const BoolToString(bool b);

#endif /* IO_HMPSIO_H_ */
