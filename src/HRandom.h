#ifndef HRANDOM_H_
#define HRANDOM_H_

class HRandom {
public:
    HRandom() {
        random_mw = 1985;
        random_mz = 2012;
    }
    int intRandom() {
        random_mz = 36969 * (random_mz & 65535) + (random_mz >> 16);
        random_mw = 18000 * (random_mw & 65535) + (random_mw >> 16);
        unsigned result = (random_mz << 16) + random_mw;
        return result >> 1;
    }
    double dblRandom() {
        random_mz = 36969 * (random_mz & 65535) + (random_mz >> 16);
        random_mw = 18000 * (random_mw & 65535) + (random_mw >> 16);
        unsigned result = (random_mz << 16) + random_mw;
        return (result + 1.0) * 2.328306435454494e-10;
    }
private:
    unsigned random_mw;
    unsigned random_mz;
};

#endif /* HRANDOM_H_ */
