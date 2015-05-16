//The following ifndef/define/endif pattern is called a 
//scope guard, and prevents the C++ compiler (actually, preprocessor)
//from including a header file more than once.
#ifndef __MY_CLASS_H_
#define __MY_CLASS_H_

#include <string>
#include <iostream>

class MyClass
    {
    //Data members: best to keep these private
    //in most designs to ensure object remains
    //in valid state and to separate interface
    //from implementation details
    int i_;
    std::string name_;
    public:
    
    //Default constructor
    MyClass();

    MyClass(std::string name, int i);

    const std::string&
    name() const;

    int
    value() const;
    };

//
//Defining this method enables printing of MyClass objects
//using cout << m << endl; where m is a MyClass instance.
//
//It also allows printing using the print,println,printf, and
//printfln functions. Use the "%s" flag to print custom objects
//with printf and printfln.
//
//See myclass.cc for implementation.
//
std::ostream&
operator<<(std::ostream& s, const MyClass& m);

//Implementation code is in myclass.cc ...

#endif
