//
// Distributed under the ITensor Library License, Version 1.2
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
    //std::string message_;
    public:
    
    using Parent = std::runtime_error;

    explicit 
    ITError(const std::string& message = "")
        : 
        Parent(message)
        { }

    //virtual 
    //const char* what() const throw()
    //    {
    //    return message_.c_str();
    //    }

    virtual
    ~ITError() { }
    }; //class ITError


inline 
std::ostream&
operator<<(std::ostream& s, const ITError& e)
    {
    s << e.what();
    return s;
    }


#endif
