#ifndef HIGHS_C_API
#define HIGHS_C_API

#ifdef __cplusplus
extern "C" {
#endif

void* Highs_create();

void Highs_destroy(void* highs);

int Highs_run(void* highs);

int Highs_loadFromFile(void* highs, const char* filename);

int Highs_setHighsOptionValue(void* highs, const char* option, const char* value);

int Highs_initializeFromFile(void* highs, const char* filename);

int Highs_writeToFile(void* highs, const char* filename);

#ifdef __cplusplus
}
#endif

#endif