/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/HighsIO.h
 * @brief IO methods for HiGHS - currently just print/log messages
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef HIGHS_IO_H
#define HIGHS_IO_H

/**
 * @brief IO methods for HiGHS - currently just print/log messages
 */
enum HighsMessageType { INFO, WARNING, ERROR };
const char* const HighsMessageTypeTag[] = {"INFO", "WARNING", "ERROR"};

unsigned const int ML_VERBOSE = 1;
unsigned const int ML_DETAILED = 2;
unsigned const int ML_MINIMAL = 4;

/**
 * @brief Used to direct printed output to FILE* output, according
 * to whether the level bit is set in messageLevel
 */
void HighsPrintMessage(
		       unsigned int level, //!< The message level: "INFO", "WARNING", "ERROR"
		       const char* format, //!< Printing format: must contain exactly one "\n" at end of format 
		       ...
		       );

/**
 * @brief Used to direct _single-line_ logging output to FILE* logfile,
 * with timestamp and message type notification
 */
void HighsLogMessage(
		     HighsMessageType type, //!< The message type
		     const char* format,    //!< Printing format
		     ...
		     );  

#endif
