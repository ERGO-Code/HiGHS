#ifndef HTESTER_H_
#define HTESTER_H_

#include "HConfig.h"
#ifdef HiGHSDEV
#include "HModel.h"
#include <string>
#include <vector>
using namespace std;

class HTester {
public:
    void setup(const char *pivotFile);
    void testUpdate(int item);
    void testCFT();
private:
    double solveTime;
    string modelName;
    int numPivot;

    vector<int> historyIn;
    vector<int> historyOut;
    vector<double> historyAlpha;

    HModel model;
};
#endif
#endif /* HTESTER_H_ */
