#ifndef FILEREADER_LP
#define FILEREADER_LP

#include "Filereader.h"

class FilereaderLp: public Filereader{
public:
    FilereaderRetcode readModelFromFile(
        const char filename,
        LpData& model
    );

};


#endif