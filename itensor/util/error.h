//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef _error_h
#define _error_h

#include <cstdlib>
#include <stdexcept>
#include <iostream>

namespace itensor{
void error(const std::string& s);
void error(const std::string& s, int line,const char* file);
#define Error(exp)  error(exp, __LINE__, __FILE__)

struct ITError : std::runtime_error
    {
    explicit 
    ITError(const std::string& message = "")
      : std::runtime_error(message) { }
    virtual
    ~ITError() { }
    };

struct ResultIsZero : ITError
    {
    ResultIsZero(const std::string& message) 
        : ITError(message)
        { }
    };



inline std::ostream&
operator<<(std::ostream& s, const ITError& e)
    {
    s << e.what();
    return s;
    }

void inline
error(const std::string& s)
    {
    std::cerr << std::endl << s << std::endl;
    std::cout << std::endl << s << std::endl;
    std::cout.flush();
    abort();
    }

void inline
error(const std::string& s, int line, const char* file = 0)
    {
    std::cerr << "From line " << line;
    if(file != 0) std::cerr << ", file " << file;
    std::cerr << std::endl;

    std::cerr << std::endl << s << std::endl;
    std::cout << std::endl << s << std::endl;
    std::cout.flush();
    std::cerr.flush();
    abort();
    }

}

#endif
