#include "myclass.h"

MyClass::
MyClass()
  : i_(0)
    { } 

MyClass::
MyClass(std::string name, int i)
  : i_(i),
    name_(name)
    { }

std::string const& MyClass::
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
operator<<(std::ostream& s, MyClass const& m)
    {
    s << "MyClass(" << m.name() << "," << m.value() << ")";
    return s;
    }


