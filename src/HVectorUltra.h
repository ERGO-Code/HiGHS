#ifndef HVECTORULTRA_H_
#define HVECTORULTRA_H_

#include <vector>
using namespace std;

class HVectorUltra {
public:
    // Main part
    void setup(int size_);
    void clear();
    int size;
    int count;   // count of non zeros
    const unsigned char ilP1 = 255;
    const unsigned short ilP2 = 65535;
    const unsigned int ilP4 = 4294967295;
    int pWd; // Bytes of pointer to values [0 => No pointers]
    vector<int> index;   // index of non zeros
    vector<double> array;   // array
    vector<unsigned char> valueP1;   // 1-byte pointer to values
    vector<unsigned short> valueP2;   // 2-byte pointer to values
    vector<unsigned int> valueP4;   // 4-byte pointer to values
    int pseudoTick;
    double fakeTick;

    // For update
    vector<char> cwork; // char working buffer
    vector<int> iwork;   // integer working buffer
    HVectorUltra *next;

    // Package
    void tight();
    void pack();
    bool packFlag;   // pack flag: do pack or not
    int packCount;   // pack count
    vector<int> packIndex;   // pack index
    vector<double> packValue;   // pack value

    // Advanced
    void copy(const HVectorUltra *from);
    double norm2();
    void saxpy(const double pivotX, const HVectorUltra *pivot);
};

typedef HVectorUltra* HVectorUltra_ptr;

#endif /* HVECTORULTRA_H_ */
