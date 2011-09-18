#ifndef _error_h
#define _error_h

#include <stdexcept>
#include <iostream>
#include <sstream>

void error(const char* s);
void error(const char* s, int line,const char* file);
#define Error(exp)  error(exp, __LINE__, __FILE__)

class ITError : public std::runtime_error
{
    typedef std::runtime_error Parent;
public:
    explicit ITError(const std::string& message)
    : Parent(message)
    { }
}; //class ITError

#endif
