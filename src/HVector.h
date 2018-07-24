#ifndef HVECTOR_H_
#define HVECTOR_H_

#include <vector>
#include "HVectorUltra.h"
using namespace std;

class HVector {
public:
    // Main part
    void setup(int size_);
    void clear();
    int size;
    int count;   // count of non zeros
    vector<int> index;   // index of non zeros
    vector<double> array;   // array
    int pseudoTick;
    double fakeTick;

    // For update
    vector<char> cwork; // char working buffer
    vector<int> iwork;   // integer working buffer
    HVector *next;

    // Package
    void tight();
    void pack();
    bool packFlag;   // pack flag: do pack or not
    int packCount;   // pack count
    vector<int> packIndex;   // pack index
    vector<double> packValue;   // pack value

    // Advanced
    void copy(const HVector *from);
    void copyUltra(const HVectorUltra *from);
    double norm2();
    void saxpy(const double pivotX, const HVector *pivot);
};

typedef HVector* HVector_ptr;

#endif /* HVECTOR_H_ */
