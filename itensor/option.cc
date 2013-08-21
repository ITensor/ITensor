//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#include "option.h"
#include <cerrno>

using namespace std;


Opt::
Opt()
    :
    name_("NullOpt"),
    bval_(false),
    rval_(NAN)
    { }

Opt::
Opt(const Name& name)
    :
    name_(name),
    bval_(true),
    rval_(NAN)
    { }

Opt::
Opt(const Name& name, bool bval)
    :
    name_(name),
    bval_(bval),
    rval_(NAN)
    { }

Opt::
Opt(const Name& name, const string& sval)
    :
    name_(name),
    bval_(true),
    sval_(sval),
    rval_(NAN)
    { }

Opt::
Opt(const Name& name, int ival)
    :
    name_(name),
    bval_(true),
    rval_(ival)
    { }

Opt::
Opt(const Name& name, Real rval)
    :
    name_(name),
    bval_(true),
    rval_(rval)
    { }

Opt::
Opt(const Name& name, 
       bool bval,
       const string& sval, 
       Real rval)
    :
    name_(name),
    bval_(bval),
    sval_(sval),
    rval_(rval)
    { }

ostream& 
operator<<(ostream & s, const Opt& opt)
    {
    s << "Opt \"" << opt.name() << "\"\n";
    s << "  boolVal   = " << (opt.boolVal() ? "true" : "false") << "\n";
    s << "  realVal   = " << opt.realVal() << "\n";
    s << "  stringVal = \"" << opt.stringVal() << "\"" 
      << endl;
    return s;
    }



OptSet::
OptSet()
    :
    is_global_(false)
    { }

OptSet::
OptSet(const Opt& opt)
    :
    is_global_(false)
    {
    add(opt);
    }

OptSet::
OptSet(const Opt& opt1, const Opt& opt2, 
       const Opt& opt3, const Opt& opt4)
    :
    is_global_(false)
    {
    add(opt1,opt2,opt3,opt4);
    }

OptSet::
OptSet(const char* ostring)
    :
    is_global_(false)
    {
    processString(string(ostring));
    }

OptSet::
OptSet(const string& ostring)
    :
    is_global_(false)
    {
    processString(ostring);
    }

OptSet::
OptSet(bool isGlobal)
    :
    is_global_(isGlobal)
    { }

OptSet::
OptSet(const OptSet& other)
    :
    is_global_(false)
    { 
    if(!other.is_global_)
        opts_ = other.opts_;
    }

bool OptSet::
defined(const Opt::Name& name) const
    {
    if(opts_.count(name) > 0)
        return true;

    if(is_global_) 
        return false;
    //else see if GlobalOpts contains it
    return GlobalOpts().defined(name);
    }

void OptSet::
add(const Opt& opt1, const Opt& opt2,
    const Opt& opt3, const Opt& opt4)
    {
    if(opt1.isNotNull())
        opts_[opt1.name()] = opt1;
    if(opt2.isNotNull())
        opts_[opt2.name()] = opt2;

    if(opt3.isNotNull())
        opts_[opt3.name()] = opt3;
    else
        return;

    if(opt4.isNotNull())
        opts_[opt4.name()] = opt4;
    }

 
const Opt& OptSet::
get(const Opt::Name& name) const
    {
    const_iterator it = opts_.find(name);
    if(it != opts_.end()) return it->second;
    //else, couldn't find the Opt
    if(is_global_)
        {
        Error("Requested option " + name + " not found");
        }
    return GlobalOpts().get(name);
    }

bool OptSet::
getBool(const Opt::Name& name) const
    {
    return get(name).boolVal();
    }

bool OptSet::
getBool(const Opt::Name& name, bool default_value) const
    {
    if(defined(name))
        return get(name).boolVal();
    else
        return default_value;
    }

 
const string& OptSet::
getString(const Opt::Name& name) const
    {
    return get(name).stringVal();
    }

const string& OptSet::
getString(const Opt::Name& name, const string& default_value) const
    {
    if(defined(name))
        return get(name).stringVal();
    else
        return default_value;
    }

int OptSet::
getInt(const Opt::Name& name) const
    {
    return get(name).intVal();
    }

int OptSet::
getInt(const Opt::Name& name, int default_value) const
    {
    if(defined(name))
        return get(name).intVal();
    else
        return default_value;
    }

Real OptSet::
getReal(const Opt::Name& name) const
    {
    return get(name).realVal();
    }

Real OptSet::
getReal(const Opt::Name& name, Real default_value) const
    {
    if(defined(name))
        return get(name).realVal();
    else
        return default_value;
    }

void OptSet::
processString(string ostring)
    {
    size_t found = ostring.find_first_of(',');
    while(found != std::string::npos)
        {
        addByString(ostring.substr(0,found));
        ostring = ostring.substr(found+1);
        found = ostring.find_first_of(',');
        }

    addByString(ostring);
    }

void OptSet::
addByString(string ostring)
    {
    if(ostring.size() < 1) return;

    const
    size_t found = ostring.find_first_of('=');

    if(found == std::string::npos)
        {
        add(Opt(ostring));
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
                add(Opt(name,d));
                return;
                }
            }
        
        if(val == "false")
            {
            add(Opt(name,false));
            }
        else
        if(val == "true")
            {
            add(Opt(name,true));
            }
        else
            {
            add(Opt(name,val));
            }

        }
    }




OptSet
operator&(const Opt& opt1, const Opt& opt2)
    {
    return OptSet(opt1,opt2);
    }

OptSet 
operator&(OptSet oset, const Opt& opt)
    {
    oset.add(opt);
    return oset;
    }

OptSet 
operator&(const Opt& opt, OptSet oset)
    {
    oset.add(opt);
    return oset;
    }
 
ostream& 
operator<<(ostream & s, const OptSet& oset)
    {
    typedef OptSet::const_iterator 
    const_it;

    if(oset.is_global_)
        s << "/- Global OptSet -------\n\n";
    else
        s << "/- OptSet (only showing overrides of global opts) -------\n\n";

    for(const_it it = oset.opts_.begin();
        it != oset.opts_.end(); ++it)
        {
        s << it->second << "\n";
        }

    s << "\\------------------" << endl;
    return s;
    }
