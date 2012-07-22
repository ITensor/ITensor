//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_OPTION_H
#define __ITENSOR_OPTION_H

#include <set>
#include "real.h"


class Option
    {
    public:

    typedef const std::string
    Name;

    Option();

    explicit
    Option(Name& name);

    Option(Name& name, bool bval);

    Option(Name& name, const std::string& sval);

    Option(Name& name, int ival);

    Option(Name& name, Real rval);

    Option(Name& name, 
           bool bval,
           const std::string& sval, 
           int ival, 
           Real rval);

    //
    // Operators for comparison and sorting
    //

    bool
    operator==(const Option& other) const
        { 
        return (name_ == other.name_ && 
                bval_ == other.bval_ &&
                sval_ == other.sval_ &&
                ival_ == other.ival_ &&
                rval_ == other.rval_);
        }

    bool
    operator<(const Option& other) const
        { 
        return (name_ < other.name_ || 
                bval_ < other.bval_ ||
                sval_ < other.sval_ ||
                ival_ < other.ival_ ||
                rval_ < other.rval_);
        }

    //
    // Accessor methods
    //

    bool
    boolVal() const { return bval_; }

    const std::string&
    stringVal() const { return sval_; }

    int
    intVal() const { return ival_; }

    Real
    realVal() const { return rval_.r; }

    bool
    isNull() const { return name_ == "NullOption"; }
    bool
    isNotNull() const { return name_ != "NullOption"; }

    Name
    name() const { return name_; }

    operator Name() const { return name_; }

    private:

    /////////////////////
    // Data Members

    Name name_;

    bool bval_;
    std::string sval_;
    int ival_;
    ApproxReal rval_;

    //
    /////////////////////

    };

inline Option::
Option()
    :
    name_("NullOption"),
    bval_(false),
    ival_(0),
    rval_(0)
    { }

inline Option::
Option(Name& name)
    :
    name_(name),
    bval_(false),
    ival_(0),
    rval_(0)
    { }

inline Option::
Option(Name& name, bool bval)
    :
    name_(name),
    bval_(bval),
    ival_(0),
    rval_(0)
    { }

inline Option::
Option(Name& name, const std::string& sval)
    :
    name_(name),
    bval_(false),
    sval_(sval),
    ival_(0),
    rval_(0)
    { }

inline Option::
Option(Name& name, int ival)
    :
    name_(name),
    bval_(false),
    ival_(ival),
    rval_(0)
    { }

inline Option::
Option(Name& name, Real rval)
    :
    name_(name),
    bval_(false),
    ival_(0),
    rval_(rval)
    { }

inline Option::
Option(Name& name, 
       bool bval,
       const std::string& sval, 
       int ival, 
       Real rval)
    :
    name_(name),
    bval_(bval),
    sval_(sval),
    ival_(ival),
    rval_(rval)
    { }


//
// OptionSet
//

class OptionSet
    {
    public:

    OptionSet();

    OptionSet(Option opt1);

    OptionSet(Option opt1, Option opt2);

    OptionSet(Option opt1, Option opt2, Option opt3);

    OptionSet(Option opt1, Option opt2, Option opt3,
              Option opt4);

    OptionSet(Option opt1, Option opt2, Option opt3,
              Option opt4, Option opt5);

    OptionSet(Option opt1, Option opt2, Option opt3,
              Option opt4, Option opt5, Option opt6);

    bool
    includes(const Option& val) const { return opts_.count(val) == 1; }
    bool
    includes(Option::Name& name) const;

    void
    insert(const Option& val) { if(val.isNotNull()) opts_.insert(val); }

    const Option&
    get(const Option& opt) const;
    const Option&
    get(Option::Name& name) const;

    bool
    boolVal(const Option& opt) const;
    bool
    boolVal(Option::Name& name) const;

    const std::string&
    stringVal(const Option& opt) const;
    const std::string&
    stringVal(Option::Name& name) const;

    int
    intVal(const Option& opt) const;
    int
    intVal(Option::Name& name) const;

    Real
    realVal(const Option& opt) const;
    Real
    realVal(Option::Name& name) const;

    private:

    std::set<Option> opts_;

    typedef std::set<Option>::iterator
    opts_it_;

    };

inline OptionSet::
OptionSet()
    { }

inline OptionSet::
OptionSet(Option opt1)
    {
    opts_.insert(opt1);
    }

inline OptionSet::
OptionSet(Option opt1, Option opt2)
    {
    insert(opt1);
    insert(opt2);
    }

inline OptionSet::
OptionSet(Option opt1, Option opt2, Option opt3)
    {
    insert(opt1);
    insert(opt2);
    insert(opt3);
    }

inline OptionSet::
OptionSet(Option opt1, Option opt2, Option opt3,
          Option opt4)
    {
    insert(opt1);
    insert(opt2);
    insert(opt3);
    insert(opt4);
    }

inline OptionSet::
OptionSet(Option opt1, Option opt2, Option opt3,
          Option opt4, Option opt5)
    {
    insert(opt1);
    insert(opt2);
    insert(opt3);
    insert(opt4);
    insert(opt5);
    }

inline OptionSet::
OptionSet(Option opt1, Option opt2, Option opt3,
          Option opt4, Option opt5, Option opt6)
    {
    insert(opt1);
    insert(opt2);
    insert(opt3);
    insert(opt4);
    insert(opt5);
    insert(opt6);
    }

bool inline OptionSet::
includes(Option::Name& name) const
    {
    Foreach(const Option& oo, opts_)
        {
        if(oo.name() == name) return true;
        }
    return false;
    }

inline const Option& OptionSet::
get(const Option& opt) const
    {
    opts_it_ it = opts_.find(opt);
    if(it == opts_.end())
        Error("OptionSet does not contain requested option");
    return *it;
    }

inline const Option& OptionSet::
get(Option::Name& name) const
    {
    Foreach(const Option& oo, opts_)
        {
        if(oo.name() == name) return oo;
        }
    Error("OptionSet does not contain requested option");
    return *(opts_.begin());
    }

bool inline OptionSet::
boolVal(const Option& opt) const
    {
    return get(opt).boolVal();
    }
bool inline OptionSet::
boolVal(Option::Name& name) const
    {
    return get(name).boolVal();
    }

inline const std::string& OptionSet::
stringVal(const Option& opt) const
    {
    return get(opt).stringVal();
    }
inline const std::string& OptionSet::
stringVal(Option::Name& name) const
    {
    return get(name).stringVal();
    }

int inline OptionSet::
intVal(const Option& opt) const
    {
    return get(opt).intVal();
    }
int inline OptionSet::
intVal(Option::Name& name) const
    {
    return get(name).intVal();
    }

Real inline OptionSet::
realVal(const Option& opt) const
    {
    return get(opt).realVal();
    }
Real inline OptionSet::
realVal(Option::Name& name) const
    {
    return get(name).realVal();
    }



//
// Convenience functions for
// creating Option instances

Option inline
Auto(bool val = true)
    {
    return Option("Auto",val);
    }

Option inline
Cutoff(int val)
    {
    return Option("Cutoff",val);
    }

Option inline
Cutoff(Real val = 0)
    {
    return Option("Cutoff",val);
    }

Option inline
DebugLevel(int level)
    {
    return Option("DebugLevel",level);
    }

Option inline
DoPinning(Real val = 1)
    {
    return Option("DoPinning",val);
    }

Option inline
NumCenter(int nc = 2)
    {
    return Option("NumCenter",nc);
    }

Option inline
Offset(int n = 0)
    {
    return Option("Offset",n);
    }

Option inline
PreserveShape()
    {
    return Option("PreserveShape");
    }

Option inline
Quiet(bool val = true)
    {
    return Option("Quiet",val);
    }

Option inline
UseWF()
    {
    return Option("UseWF");
    }

Option inline
Verbose(bool val = true)
    {
    return Option("Verbose",val);
    }

Option inline
Weight(Real w = 1)
    {
    return Option("Weight",w);
    }

#endif
