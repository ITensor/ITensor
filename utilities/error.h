#ifndef _error_h
#define _error_h

#include <exception>
#include <iostream>
#include <cstdlib>
#include <iomanip>

using std::cout;
using std::cerr;
using std::endl;

void error(const char* s);
void error(const char* s, int line,const char* file);
#define Error(exp)  error(exp, __LINE__, __FILE__)

#endif
