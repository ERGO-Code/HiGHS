#ifndef HIGHS_IO_H
#define HIGHS_IO_H

enum HighsMessageType { INFO, WARNING, ERROR };
const char* const HighsMessageTypeTag[] = {"INFO", "WARNING", "ERROR"};

unsigned const int ML_VERBOSE = 1;
unsigned const int ML_DETAILLED = 2;
unsigned const int ML_MINIMAL = 4;

void HighsPrintMessage(unsigned int level, const char* format, ...);

void HighsLogMessage(HighsMessageType type, const char* format, ...);

#endif