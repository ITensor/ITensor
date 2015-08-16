//The following ifndef/define/endif pattern is called a 
//scope guard, and prevents the C++ compiler (actually, preprocessor)
//from including a header file more than once.
#ifndef __MY_CLASS_H_
#define __MY_CLASS_H_

#include <string>
#include <iostream>

class MyClass
    {
    // Data members: best to keep these private
    // in most designs to ensure object remains
    // in valid state and to separate interface
    // from implementation details
    int i_;
    std::string name_;
    public:
    
    //Default constructor
    MyClass();

    //Typical constructor
    MyClass(std::string name, int i);

    //Methods for read-only
    //access to data members:

    std::string const&
    name() const;  //<-  const at the
                   //    end means methods
    int            //    do not modify
    value() const; //<-  this object
                   //    and can be used
                   //    when object is const
    };

//
// Defining this method enables printing of MyClass objects
// using cout << m << endl; where m is a MyClass instance.
//
// It also allows printing using the print,println,printf, and
// printfln functions. Use the "%s" flag to print custom objects
// with printf and printfln.
//
// See myclass.cc for implementation.
//
std::ostream&
operator<<(std::ostream& s, MyClass const& m);
//                                    ^
// function takes MyClass as a        ^
// const& to avoid making a copy

// Implementation code is in myclass.cc ...

#endif
