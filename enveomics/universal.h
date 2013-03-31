// enveomics/universal.h - Library for all enve-omics software
// @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
// @license artistic 2.0
// @version 1.0

#ifndef ENVEOMICS_H
#define ENVEOMICS_H

#include <iostream>
#include <string>
#include <stdarg.h>
#include <cstring>
#include <cstdlib>
#include <climits>

using namespace std;

/**
 * void error(const char *msg);
 * void error(const char *msg, const char *val);
 * void error(const char *msg, int val);
 * void error(const char *msg, unsigned int val);
 * void error(const char *msg, double val);
 * void error(const char *msg, string val);
 * Description:
 *   Produces an error message and terminates the program.
 * Input:
 *   const char *msg: The error message.
 *   mix val (optional): A value to be reported accompannying the error message.
 */
void error(const char *msg);
void error(const char *msg, const char *val);
void error(const char *msg, int val);
void error(const char *msg, unsigned int val);
void error(const char *msg, double val);
void error(const char *msg, string val);

/**
 * void set_verbosity(int v);
 * Description:
 *   Globally sets the verbosity level.  The higher the verbosity level, the noisier the
 *   program.  Some hallmarks in the verbosity:
 *    o 0: no output.
 *    o 1: only basic output.
 *    o 4: most important information is shown.
 *    o 7: all the information.
 *    o 9: debugging data is shown.
 * Input:
 *   int v: The new verbosity level.
 */
void set_verbosity(int v);

/**
 * int get_verbosity();
 * Description:
 *   Returns the globally set verbosity level.
 * Output:
 *   int: The current verbosity level.
 */
int get_verbosity();

/**
 * void open_log(char *logfile);
 * Description:
 *   Opens a log file to save all the messages passed via say().
 * Input:
 *   char *logfile: Path to the log file to be created.
 */
void open_log(char *logfile);

/**
 * void close_log();
 * Description:
 *   Closes the log file, if open.
 */
void close_log();

/**
 * bool log_is_open();
 * Description:
 *   Asserts whether the log is open or not.
 * Output:
 *   bool: Is it open?
 */
bool log_is_open();

/**
 * void say(const char *format, ...);
 * Description:
 *   Says something (prints to the stderr).
 * Input:
 *   const char *format:  A string describing the following parameters.  The format must be:
 *      /([0-9])(\!?)([siuf]*)([$^]?)/.  The elements captured by the regexp parenthesis are:
 *        o A digit determining the minimum level of verbosity at which this message must be printed.
 *        o A bang (optional).  If passed, indicates that the string must be printed "as is", without the
 *          time prefix.  This is useful if you want to print a message without carriage return, to be
 *          completed in a further call of say().
 *        o Zero or more characters indicating the type of arguments (in the same order as passed).  Supported
 *          types are:  int (i), size_t (u), double (f), char (c) and char* (s).
 *        o A finalization character (optional).  If the finalization character is '$', a carriage return is
 *          printed at the end.  If the finalization character is '^', the pointer in the terminal is moved to the
 *          left.  The '^' basically means that say() must print some spaces (to whipe out the line) and a '\r'
 *          character.  If the finalization character is '>', some spaces and a carriage return are printed.  This
 *          is useful to clean the line (with some closure message) after one or more messages finished with '^'.
 *          NOTE: If the verbosity is 9, all the termination characters are interpreted as '$'.
 */
void say(const char *format, ...);

#endif

