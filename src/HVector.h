#ifndef HVECTOR_H_
#define HVECTOR_H_

#include <vector>
#include <set>
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

    //For Ultra-sparsity
    const int dfSparseDaStr = -1;
    const int p0SparseDaStr = 0;
    const int p1SparseDaStr = 1;
    const int p2SparseDaStr = 2;
    const int           mxSetP0 =  1024;//10;
    const unsigned char ilP1 =     255;//10;
    const unsigned short ilP2 =    65535;//20;
    int numEnSetP0; // Number of entries in setP0
    int pWd; // Bytes of pointer to values
    // -1 => Vanilla index-array sparse data structure
    //  0 => Ultra-sparse set of indices without pointers
    //  1 => Ultra-sparse set of indices with 1-byte pointers
    //  2 => Ultra-sparse set of indices with 2-byte pointers
    set<pair<int, double>> setP0; //set of index-array pairs 
    vector<unsigned char> valueP1;   // 1-byte pointers
    vector<unsigned short> valueP2;   // 2-byte pointers

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
    double norm2();
    void saxpy(const double pivotX, const HVector *pivot);
};

typedef HVector* HVector_ptr;

#endif /* HVECTOR_H_ */
