// enveomics/universal.h - Library for all enve-omics software
// @author Luis M. Rodriguez-R <lmrodriguezr at gmail dot com>
// @license artistic 2.0
// @version 1.0

#include <iostream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <stdarg.h>
#include <time.h>

#include "universal.h"

#define LARGEST_STRING 16384

using namespace std;

int		Verbosity = 1;
ofstream	LogH;
bool		OpenLog;

void error(const char *msg, const char *val){
   cerr << "Fatal error:" << endl << msg;
   if (val!=NULL && strlen(val) > 0) cerr << ": " << endl << val;
   cerr << endl;
   say("0ss", "Fatal error: ", msg);
   if (val!=NULL && strlen(val) > 0) say("0!ss", ": ", val);
   say("0!$");
   exit(1);
}

void error(const char *msg){
   return error(msg, (const char *)NULL);
}

void error(const char *msg, int val){
   char *valChr = new char[100];
   snprintf(valChr, 100, "%d", val);
   return error(msg, valChr);
}

void error(const char *msg, double val){
   char *valChr = new char[100];
   snprintf(valChr, 100, "%.5f", val);
   return error(msg, valChr);
}

void error(const char *msg, unsigned int val){
   char *valChr = new char[100];
   snprintf(valChr, 100, "%d", val);
   return error(msg, valChr);
}

void error(const char *msg, string val){
   char valChr[val.size()];
   for(size_t a=0; a<=val.size(); a++) valChr[a]=val[a];
   return error(msg, valChr);
}

void set_verbosity(int v){ Verbosity = v; }

int get_verbosity(){ return Verbosity; }

void open_log(char *logfile){
   if(OpenLog) close_log();
   LogH.open(logfile, ios::out);
   if(!LogH.is_open()) error("Cannot open the file", logfile);
   OpenLog = true;
}

void close_log(){
   if(!OpenLog) return;
   if(LogH && LogH.is_open()) LogH.close();
   OpenLog = false;
}

bool log_is_open() { return OpenLog; }

void say(const char *format, ...){
   va_list	arguments;
   va_start(arguments, format);
   int		level, startArg=1;
   char		*sArg, levelChr[2];
   string	out;
   
   levelChr[0] = format[0];
   levelChr[1] = '\0';
   level = atoi(levelChr);
   if(level > Verbosity) return;
   
   if(format[startArg]=='!') startArg++;
   else{
      sArg = new char[LARGEST_STRING];
      snprintf(sArg, LARGEST_STRING, " [% 9.1f] ",
               clock() / (60.0 * CLOCKS_PER_SEC));
      out.append(sArg);
      for(int i=0; i<level; i++) out.append(" ");
   }
   
   for(int i=startArg; format[i] != '\0'; i++){
      sArg = new char[LARGEST_STRING];
      // Known types
      if(format[i] == 'f') {
         snprintf(sArg, LARGEST_STRING, "%f", va_arg(arguments, double));
      } else if (format[i] == 'i') {
         snprintf(sArg, LARGEST_STRING, "%d", va_arg(arguments, int));
      } else if (format[i] == 'u') {
	 snprintf(sArg, LARGEST_STRING, "%u", va_arg(arguments, unsigned int));
      } else if (format[i] == 's') {
         sArg = va_arg(arguments, char *);
      } else if (format[i] == 'c') {
         snprintf(sArg, LARGEST_STRING, "%c", va_arg(arguments, int));
      
      } else if (format[i] == '$') {
         // Termination characters
         out.append("\n");
	 goto print;
      } else if (format[i] == '^' || format[i] == '>') {
         if (Verbosity == 9) {
	    out.append("\n");
	    goto print;
	 } else {
	    LogH << out << '\n';
	    cerr << out;
	    for(int a=0; a<30; a++) cerr << " ";
	    if(format[i]=='^') cerr << '\r';
	    else cerr << endl;
	    goto end;
	 }
      
      // Unknown type
      }else{
         sArg = (char *)"";
      }
      
      out.append(sArg);
   }

   print:
   if(OpenLog) LogH << out;
   cerr << out;

   end:
   va_end(arguments);
}

