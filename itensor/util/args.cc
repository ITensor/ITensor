//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include <cerrno>
#include <algorithm>
#include <iostream>
#include "itensor/util/args.h"
#include "itensor/util/error.h"
#include "itensor/util/readwrite.h"

namespace itensor {

using std::string;
using std::ostream;
using std::istream;

std::string
chopSpaceEq(std::string name)
    {
    auto s = name.size();
    while(s > 0)
        {
        if(name[s-1]=='=' || name[s-1]==' ')
            {
            --s;
            name.resize(s);
            }
        else
            {
            break;
            }
        }
    return name;
    }


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
    name_(chopSpaceEq(name)),
    type_(Boolean),
    rval_(1.0)
    { }

Args::Val::
Val(Name const& name)
    :
    name_(chopSpaceEq(name)),
    type_(Boolean),
    rval_(1.0)
    { }

Args::Val::
Val(Name const& name, bool bval)
    :
    name_(chopSpaceEq(name)),
    type_(Boolean),
    rval_((bval ? 1.0 : 0.0))
    { }

Args::Val::
Val(Name const& name, const char* sval)
    :
    name_(chopSpaceEq(name)),
    type_(String),
    sval_(sval),
    rval_(NAN)
    { }

Args::Val::
Val(Name const& name, const string& sval)
    :
    name_(chopSpaceEq(name)),
    type_(String),
    sval_(sval),
    rval_(NAN)
    { }

Args::Val::
Val(Name const& name, long ival)
    :
    name_(chopSpaceEq(name)),
    type_(Numeric),
    rval_(ival)
    { }


Args::Val::
Val(Name const& name, int ival)
    :
    name_(chopSpaceEq(name)),
    type_(Numeric),
    rval_(ival)
    { }

Args::Val::
Val(Name const& name, unsigned long ival)
    :
    name_(chopSpaceEq(name)),
    type_(Numeric),
    rval_(ival)
    { }

Args::Val::
Val(Name const& name, unsigned int ival)
    :
    name_(chopSpaceEq(name)),
    type_(Numeric),
    rval_(ival)
    { }

Args::Val::
Val(Name const& name, Real rval)
    :
    name_(chopSpaceEq(name)),
    type_(Numeric),
    rval_(rval)
    { }

void Args::Val::
read(std::istream& s)
    { 
    itensor::read(s, name_);
    itensor::read(s, type_);
    if(type_ == String)
        itensor::read(s, sval_);
    else
        itensor::read(s, rval_);
    }

void Args::Val::
write(std::ostream& s) const
    { 
    itensor::write(s, name_);
    itensor::write(s, type_);
    if(type_ == String)
        itensor::write(s, sval_);
    else
        itensor::write(s, rval_);
    }

void Args::Val::
assertType(Type t) const
    {
    if(t != type_)
        throw ITError("Wrong value type for option " + name_);
    }

ostream& 
operator<<(ostream & s, Args::Val const& v)
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
Args(Args const& other)
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
add(Name const& name, bool bval) { add({name,bval}); }
void Args::
add(Name const& name, long ival) { add({name,ival}); }
void Args::
add(Name const& name, int ival) { add({name,ival}); }
void Args::
add(Name const& name, const char* sval) { add({name,std::string(sval)}); }
void Args::
add(Name const& name, const std::string& sval) { add({name,sval}); }
void Args::
add(Name const& name, Real rval) { add({name,rval}); }

bool Args::
defined(Name const& name) const
    {
    for(auto& x : vals_)
        {
        if(x.name() == name) return true;
        }

    if(isGlobal()) return false;

    //otherwise see if global Args contains it
    return global().defined(name);
    }

// Remove an arg from the set - always succeeds
void Args::
remove(const Name& name)
    {
    for(auto it = vals_.begin(); it != vals_.end(); ++it)
        if(it->name() == name)
            {
            vals_.erase(it);
            break;
            }
    }


void Args::
add(Val const& val)
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
get(Name const& name) const
    {
    for(auto& x : vals_)
        {
        if(x.name() == name) return x;
        }
    //couldn't find the Val in this Args
    if(isGlobal())
        {
        throw ITError("Requested option " + name + " not found");
        }
    return global().get(name);
    }

bool Args::
getBool(Name const& name) const
    {
    return get(name).boolVal();
    }

bool Args::
getBool(Name const& name, bool default_value) const
    {
    if(defined(name)) return get(name).boolVal();
    return default_value;
    }

 
string const& Args::
getString(Name const& name) const
    {
    return get(name).stringVal();
    }

string const& Args::
getString(Name const& name, string const& default_value) const
    {
    if(defined(name)) return get(name).stringVal();
    return default_value;
    }

long Args::
getInt(Name const& name) const
    {
    return get(name).intVal();
    }

long Args::
getInt(Name const& name, long default_value) const
    {
    if(defined(name)) return get(name).intVal();
    return default_value;
    }

size_t Args::
getSizeT(Name const& name) const
    {
    return get(name).size_tVal();
    }

size_t Args::
getSizeT(Name const& name, long default_value) const
    {
    if(defined(name)) return get(name).size_tVal();
    return default_value;
    }

Real Args::
getReal(Name const& name) const
    {
    return get(name).realVal();
    }

Real Args::
getReal(Name const& name, Real default_value) const
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
operator+=(Args const& args)
    {
    for(auto& x : args.vals_)
        {
        add(x);
        }
    return *this;
    }

void Args::
read(std::istream& s)
    {
    itensor::read(s,vals_);
    }

void Args::
write(std::ostream& s) const
    {
    itensor::write(s,vals_);
    }

Args
operator+(Args args, Args const& other)
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
operator<<(ostream & s, Args const& args)
    {
    if(args.isGlobal()) s << "Global Args:\n";
    else                s << "Args: (only showing overrides of global args)\n";

    for(auto& opt : args.vals_)
        s << opt << "\n";

    return s;
    }

} //namespace itensor
