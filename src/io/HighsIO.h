#ifndef HIGHS_IO_H
#define HIGHS_IO_H

enum class HighsMessageType { DEBUG, INFO, WARNING, ERROR };

void HighsPrintMessage(HighsMessageType type, const char* format, ...);

#endif