//The following ifndef/define/endif pattern is called a 
//scope guard, and prevents the C++ compiler (actually, preprocessor)
//from including a header file more than once.
#ifndef __MY_CLASS_H_
#define __MY_CLASS_H_

#include <string>
#include <iostream>

class MyClass
    {
    public:
    
    //Default constructor
    MyClass();

    MyClass(std::string name, int i);

    const std::string&
    name() const;

    int
    i() const;

    private:

    ///////////

    int i_;
    std::string name_;

    //////////

    };

//
//Defining this method enables printing of MyClass objects
//using cout << m << endl; where m is a MyClass instance.
//See myclass.cc for implementation.
//
std::ostream&
operator<<(std::ostream& s, const MyClass& m);

//Implementation code is in myclass.cc ...

#endif
