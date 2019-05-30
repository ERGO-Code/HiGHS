#include "Highs.h"
#include "io/Filereader.h"
#include "highs_c_api.h"

void* Highs_create() {
   return new Highs();
}

void Highs_destroy(void* highs) {
   delete (Highs*)highs;
}

int Highs_run(void* highs) {
   return (int)((Highs*)highs)->run();
}

int Highs_loadFromFile(void* highs, const char* filename) {
   Filereader* reader = Filereader::getFilereader(filename);
   HighsLp lp;
   HighsOptions options;
   options.filename = std::string(filename);
   reader->readModelFromFile(options, lp);
   return (int)((Highs*)highs)->initializeLp(lp);
}