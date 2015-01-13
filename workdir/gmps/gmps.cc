#include "gmps.h"

MyClass::
MyClass()
    :
    i_(0)
    { } 

MyClass::
MyClass(std::string name, int i)
    :
    i_(i),
    name_(name)
    { }

const std::string& MyClass::
name() const
    {
    return name_;
    }

int MyClass::
value() const
    {
    return i_;
    }

std::ostream&
operator<<(std::ostream& s, const MyClass& m)
    {
    s << "MyClass(" << m.name() << "," << m.value() << ")";
    return s;
    }


