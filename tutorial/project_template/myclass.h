//The following ifndef/define/endif pattern is called a 
//scope guard, and prevents the C++ compiler (actually, preprocessor)
//from including a header file more than once.
#ifndef __MY_CLASS_H_
#define __MY_CLASS_H_

#include <string>
#include <iostream>

class MyClass
    {
    int i_;
    std::string name_;
    public:
    
    //Default constructor
    MyClass() : i_(0) { }

    MyClass(std::string name, int i)
      : i_(i),
        name_(name)
        { }

    std::string const&
    name() const { return name_; }

    int
    value() const { return i_; }

    };

//
//Defining this method enables printing of MyClass objects
//using cout << m << endl; where m is a MyClass instance.
//
//It also allows printing using the print,println,printf, and
//printfln functions defined by ITensor. 
//Use the "%s" token to print custom objects such as a MyClass
//object with printf and printfln.
//
inline std::ostream&
operator<<(std::ostream& s, MyClass const& m)
    {
    s << "MyClass(" << m.name() << "," << m.value() << ")";
    return s;
    }


#endif //__MY_CLASS_H
