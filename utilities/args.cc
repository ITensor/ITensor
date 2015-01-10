//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "args.h"
#include <cerrno>
#include <algorithm>
#include <iostream>
#include "error.h"

namespace itensor {

using std::string;
using std::cout;
using std::endl;
using std::ostream;
using std::istream;


Args::Val::
Val()
    :
    name_("Null"),
    type_(None),
    rval_(NAN)
    { }

Args::Val::
Val(const char* name)
    :
    name_(name),
    type_(Boolean),
    rval_(1.0)
    { }

Args::Val::
Val(const Name& name)
    :
    name_(name),
    type_(Boolean),
    rval_(1.0)
    { }

Args::Val::
Val(const Name& name, bool bval)
    :
    name_(name),
    type_(Boolean),
    rval_((bval ? 1.0 : 0.0))
    { }

Args::Val::
Val(const Name& name, const char* sval)
    :
    name_(name),
    type_(String),
    sval_(sval),
    rval_(NAN)
    { }

Args::Val::
Val(const Name& name, const string& sval)
    :
    name_(name),
    type_(String),
    sval_(sval),
    rval_(NAN)
    { }

Args::Val::
Val(const Name& name, int ival)
    :
    name_(name),
    type_(Numeric),
    rval_(ival)
    { }

Args::Val::
Val(const Name& name, Real rval)
    :
    name_(name),
    type_(Numeric),
    rval_(rval)
    { }

void Args::Val::
assertType(Type t) const
    {
    if(t != type_)
        throw ITError("Wrong value type for option " + name_);
    }

ostream& 
operator<<(ostream & s, const Args::Val& v)
    {
    s << v.name() << "=";
    if(v.type() == Args::Val::Boolean)
        {
        s << (v.boolVal() ? "true" : "false");
        }
    else
    if(v.type() == Args::Val::Numeric)
        {
        s << v.realVal();
        }
    else
    if(v.type() == Args::Val::String)
        {
        s << "\"" << v.stringVal() << "\"";
        }
    else
        {
        s << "(Null)";
        }
    return s;
    }



Args::
Args()
    { }

Args::
Args(const char* ostring)
    {
    processString(string(ostring));
    }

Args::
Args(const string& ostring)
    {
    processString(ostring);
    }

Args::
Args(const Args& other)
    { 
    if(!other.isGlobal())
        vals_ = other.vals_;
    }

Args::
Args(Args&& other)
    { 
    if(!other.isGlobal())
        vals_ = std::move(other.vals_);
    }

Args& Args::
operator=(const Args& other)
    { 
    if(!other.isGlobal())
        vals_ = other.vals_;
    return *this;
    }

Args& Args::
operator=(Args&& other)
    { 
    if(!other.isGlobal())
        vals_ = std::move(other.vals_);
    return *this;
    }

void Args::
add(const Name& name, bool bval) { add({name,bval}); }
void Args::
add(const Name& name, int ival) { add({name,ival}); }
void Args::
add(const Name& name, const std::string& sval) { add({name,sval}); }
void Args::
add(const Name& name, Real rval) { add({name,rval}); }

bool Args::
defined(const Name& name) const
    {
    for(const auto& x : vals_)
        {
        if(x.name() == name) return true;
        }

    if(isGlobal()) return false;

    //otherwise see if GlobalArgs contains it
    return Global().defined(name);
    }


void Args::
add(const Val& val)
    {
    if(!val) return;
    for(auto& x : vals_)
        //If already defined, replace
        if(x.name() == val.name()) 
            {
            x = val;
            return;
            }
    //Otherwise add to the end
    vals_.push_back(val);
    }

void Args::
add(const char* ostring)
    {
    processString(std::string(ostring));
    }

 
const Args::Val& Args::
get(const Name& name) const
    {
    for(const auto& x : vals_)
        {
        if(x.name() == name) return x;
        }
    //couldn't find the Val in this Args
    if(isGlobal())
        {
        throw ITError("Requested option " + name + " not found");
        }
    return Global().get(name);
    }

bool Args::
getBool(const Name& name) const
    {
    return get(name).boolVal();
    }

bool Args::
getBool(const Name& name, bool default_value) const
    {
    if(defined(name)) return get(name).boolVal();
    return default_value;
    }

 
const string& Args::
getString(const Name& name) const
    {
    return get(name).stringVal();
    }

const string& Args::
getString(const Name& name, const string& default_value) const
    {
    if(defined(name)) return get(name).stringVal();
    return default_value;
    }

int Args::
getInt(const Name& name) const
    {
    return get(name).intVal();
    }

int Args::
getInt(const Name& name, int default_value) const
    {
    if(defined(name)) return get(name).intVal();
    return default_value;
    }

Real Args::
getReal(const Name& name) const
    {
    return get(name).realVal();
    }

Real Args::
getReal(const Name& name, Real default_value) const
    {
    if(defined(name)) return get(name).realVal();
    return default_value;
    }

void Args::
processString(string ostring)
    {
    ostring.erase(std::remove(ostring.begin(), ostring.end(),' '), ostring.end());

    auto found = ostring.find_first_of(',');
    while(found != std::string::npos)
        {
        addByString(ostring.substr(0,found));
        ostring = ostring.substr(found+1);
        found = ostring.find_first_of(',');
        }

    addByString(ostring);
    }

void Args::
addByString(string ostring)
    {
    if(ostring.size() < 1) return;

    auto found = ostring.find_first_of('=');

    if(found == std::string::npos)
        {
        //if no '=' found, just create an Val by name only
        //which is the same as name=true
        add(Val(ostring));
        }
    else
        {
        string name = ostring.substr(0,found);
        string val = ostring.substr(found+1);
        if(name.size() < 1 || val.size() < 1) return;

        char f = val.at(0);
        if(f == '1' || f == '-' || isdigit(f) || f == '+')
            {
            //Try Real conversion
            errno = 0;
            char* end;
            Real d = strtod(val.c_str(),&end);

            if(errno == ERANGE)
                {
                throw ITError("Real out of range");
                }

            if(errno == 0) //success
                {
                add(Val(name,d));
                return;
                }
            }
        
        if(val == "false")     add(Val(name,false));
        else if(val == "true") add(Val(name,true));
        else                   add(Val(name,val));
        }
    }

Args& Args::
operator+=(const Args& args)
    {
    for(const auto& x : args.vals_)
        {
        add(x);
        }
    return *this;
    }


Args
operator+(Args args, const Args& other)
    {
    args += other;
    return args;
    }


Args 
operator+(Args args, const char* ostring)
    {
    args.add(ostring);
    return args;
    }

Args 
operator+(const char* ostring, Args args)
    {
    args.add(ostring);
    return args;
    }
 
ostream& 
operator<<(ostream & s, const Args& args)
    {
    if(args.isGlobal()) s << "Global Args:\n";
    else                s << "Args: (only showing overrides of global args)\n";

    for(const auto& opt : args.vals_)
        s << opt << "\n";

    return s;
    }

}; //namespace itensor
