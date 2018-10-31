#ifndef IO_FILEREADER_LP_H_
#define IO_FILEREADER_LP_H_

#include "Filereader.h"

class FilereaderLp: public Filereader{
public:
    FilereaderRetcode readModelFromFile(
        const char filename//,
       // HighsLp& model
    );

};


#endif