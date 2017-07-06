#ifndef HDUALROW_H_
#define HDUALROW_H_

#include "HVector.h"
#include "HModel.h"

#include <set>
#include <vector>
using namespace std;

/**
 * This class deal with ratio test
 * and some update dual/flip tasks.
 */

class HDualRow {
public:
    void setup(HModel *model);
    void setupSlice(HModel *model, int size);
    void clear();
    void choose_makepack(const HVector *row, const int offset);
    void choose_possible();
    void choose_joinpack(const HDualRow* otherRow);
    void choose_final();

    void update_flip(HVector *bfrtColumn);
    void update_dual(double theta);

    void create_Freelist();
    void create_Freemove(HVector *row_ep);
    void delete_Freemove();
    void delete_Freelist(int iColumn);
    void rp_hsol_pv_r();

    HModel *workModel;
    int workSize;
    const int *workRand;
    const int *workMove;
    const double *workDual;
    const double *workRange;

    set<int> freeList;

    int packCount;
    vector<int> packIndex;
    vector<double> packValue;

    double workDelta;
    double workAlpha;
    double workTheta;
    int workPivot;
    int workCount;

    vector<pair<int, double> > workData;
    vector<int> workGroup;
};

#endif /* HDUALROW_H_ */
