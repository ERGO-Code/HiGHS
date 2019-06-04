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

int Highs_setHighsOptionValue(void* highs, const char* option, const char* value) {
   return (int)((Highs*)highs)->setHighsOptionValue(std::string(option), std::string(value));
}

int Highs_initializeFromFile(void* highs, const char* filename) {
   return (int)((Highs*)highs)->initializeFromFile(std::string(filename));
}

int Highs_writeToFile(void* highs, const char* filename) {
   return (int)((Highs*)highs)->writeToFile(std::string(filename));
}