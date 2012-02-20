//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef _error_h
#define _error_h

#include <stdexcept>
#include <iostream>

void error(const std::string& s);
void error(const std::string& s, int line,const char* file);
#define Error(exp)  error(exp, __LINE__, __FILE__)

class ITError : public std::runtime_error
    {
    public:

    typedef std::runtime_error 
    Parent;

    explicit 
    ITError(const std::string& message)
        : Parent(message)
        { }

    }; //class ITError

#endif
