#ifndef FILEREADER_H
#define FILEREADER_H

#include "../lp_data/LpData.h"

enum class FilereaderRetcode {
    OKAY = 0,
    FILENOTFOUND = 1,
    PARSERERROR = 2
};

class Filereader {
    public:
    virtual FilereaderRetcode readModelFromFile(
        const char filename,
        LpData& model
    );
};


#endif