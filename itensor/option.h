//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_OPTION_H
#define __ITENSOR_OPTION_H

#include "matrix.h"
#include <map>


class Option
    {
    public:

    typedef std::string
    Name;

    Option();

    explicit
    Option(const Name& name);

    Option(const Name& name, bool bval);

    Option(const Name& name, const std::string& sval);

    Option(const Name& name, int ival);

    Option(const Name& name, Real rval);

    Option(const Name& name, 
           bool bval,
           const std::string& sval, 
           int ival, 
           Real rval);

    //
    // Operators for comparison and sorting
    //

    // Two Options are equal if they have the same name,
    // regardless of other fields that may be set.
    bool
    operator==(const Option& other) const
        { 
        return name_ == other.name_;
        }

    // Compares two Options based on their name, using
    // the < operator of std::string. Useful for sorting.
    bool
    operator<(const Option& other) const
        { 
        return name_ < other.name_;
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
    realVal() const { return rval_; }

    bool
    isNull() const { return name_ == "NullOption"; }
    bool
    isNotNull() const { return name_ != "NullOption"; }

    const Name&
    name() const { return name_; }

    operator const Name&() const { return name_; }

    friend std::ostream& 
    operator<<(std::ostream & s, const Option& opt);

    private:

    /////////////////////
    // Data Members

    Name name_;

    bool bval_;
    std::string sval_;
    int ival_;
    Real rval_;

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
Option(const Name& name)
    :
    name_(name),
    bval_(true),
    ival_(0),
    rval_(0)
    { }

inline Option::
Option(const Name& name, bool bval)
    :
    name_(name),
    bval_(bval),
    ival_(0),
    rval_(0)
    { }

inline Option::
Option(const Name& name, const std::string& sval)
    :
    name_(name),
    bval_(true),
    sval_(sval),
    ival_(0),
    rval_(0)
    { }

inline Option::
Option(const Name& name, int ival)
    :
    name_(name),
    bval_(true),
    ival_(ival),
    rval_(0)
    { }

inline Option::
Option(const Name& name, Real rval)
    :
    name_(name),
    bval_(true),
    ival_(0),
    rval_(rval)
    { }

inline Option::
Option(const Name& name, 
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

inline std::ostream& 
operator<<(std::ostream & s, const Option& opt)
    {
    s << "Option \"" << opt.name() << "\"\n";
    s << "  boolVal   = " << (opt.boolVal() ? "true" : "false") << "\n";
    s << "  intVal    = " << opt.intVal() << "\n";
    s << "  realVal   = " << opt.realVal() << "\n";
    s << "  stringVal = \"" << opt.stringVal() << "\"" 
      << std::endl;
    return s;
    }


//
// OptionSet
//

class OptionSet
    {
    public:

    OptionSet();

    OptionSet(const Option& opt1);

    OptionSet(const Option& opt1, const Option& opt2);

    OptionSet(const Option& opt1, const Option& opt2, const Option& opt3);

    OptionSet(const Option& opt1, const Option& opt2, 
              const Option& opt3, const Option& opt4);

    OptionSet(const Option& opt1, const Option& opt2, 
              const Option& opt3, const Option& opt4, const Option& opt5);

    OptionSet(const Option& opt1, const Option& opt2, const Option& opt3,
              const Option& opt4, const Option& opt5, const Option& opt6);

    //
    // Methods for accessing Options
    //

    bool
    defined(const Option::Name& name) const;
    bool
    defined(const Option& opt) const;

    void
    insert(const Option& opt) { if(opt.isNotNull()) opts_[opt.name()] = opt; }

    const Option&
    get(const Option::Name& name) const;
    const Option&
    get(const Option& opt) const;

    //
    // Methods for getting fields of a specific Option
    //

    bool
    boolVal(const Option::Name& name) const;
    bool
    boolVal(const Option& opt) const;

    bool
    boolOrDefault(const Option::Name& name, bool default_val) const;

    const std::string&
    stringVal(const Option::Name& name) const;
    const std::string&
    stringVal(const Option& opt) const;

    const std::string&
    stringOrDefault(const Option::Name& name, const std::string& default_val) const;

    int
    intVal(const Option::Name& name) const;
    int
    intVal(const Option& opt) const;

    int
    intOrDefault(const Option::Name& name, int default_val) const;

    Real
    realVal(const Option::Name& name) const;
    Real
    realVal(const Option& opt) const;

    Real
    realOrDefault(const Option::Name& name, Real default_val) const;

    friend std::ostream& 
    operator<<(std::ostream & s, const OptionSet& oset);

    private:

    typedef std::map<Option::Name,Option>
    storage_type;

    typedef storage_type::iterator
    iterator;

    typedef storage_type::const_iterator
    const_iterator;

    ///////////////
    //
    // Data Members

    storage_type opts_;

    //
    ///////////////

    };

inline OptionSet::
OptionSet()
    { }

inline OptionSet::
OptionSet(const Option& opt1)
    {
    insert(opt1);
    }

inline OptionSet::
OptionSet(const Option& opt1, const Option& opt2)
    {
    insert(opt1);
    insert(opt2);
    }

inline OptionSet::
OptionSet(const Option& opt1, const Option& opt2, const Option& opt3)
    {
    insert(opt1);
    insert(opt2);
    insert(opt3);
    }

inline OptionSet::
OptionSet(const Option& opt1, const Option& opt2, const Option& opt3,
          const Option& opt4)
    {
    insert(opt1);
    insert(opt2);
    insert(opt3);
    insert(opt4);
    }

inline OptionSet::
OptionSet(const Option& opt1, const Option& opt2, const Option& opt3,
          const Option& opt4, const Option& opt5)
    {
    insert(opt1);
    insert(opt2);
    insert(opt3);
    insert(opt4);
    insert(opt5);
    }

inline OptionSet::
OptionSet(const Option& opt1, const Option& opt2, const Option& opt3,
          const Option& opt4, const Option& opt5, const Option& opt6)
    {
    insert(opt1);
    insert(opt2);
    insert(opt3);
    insert(opt4);
    insert(opt5);
    insert(opt6);
    }

bool inline OptionSet::
defined(const Option::Name& name) const
    {
    return opts_.count(name);
    }

bool inline OptionSet::
defined(const Option& opt) const
    {
    return defined(opt.name());
    }

inline const Option& OptionSet::
get(const Option& opt) const
    {
    return get(opt.name());
    }

inline const Option& OptionSet::
get(const Option::Name& name) const
    {
    const_iterator it = opts_.find(name);
    if(it != opts_.end()) return it->second;
    //else, couldn't find the Option
    std::cout << "Option name = " << name << std::endl;
    Error("OptionSet does not contain requested option");
    return opts_.begin()->second;
    }

bool inline OptionSet::
boolVal(const Option& opt) const
    {
    return get(opt).boolVal();
    }
bool inline OptionSet::
boolVal(const Option::Name& name) const
    {
    return get(name).boolVal();
    }

bool inline OptionSet::
boolOrDefault(const Option::Name& name, bool default_value) const
    {
    if(defined(name))
        return get(name).boolVal();
    else
        return default_value;
    }

inline const std::string& OptionSet::
stringVal(const Option& opt) const
    {
    return get(opt).stringVal();
    }
inline const std::string& OptionSet::
stringVal(const Option::Name& name) const
    {
    return get(name).stringVal();
    }

inline const std::string& OptionSet::
stringOrDefault(const Option::Name& name, const std::string& default_value) const
    {
    if(defined(name))
        return get(name).stringVal();
    else
        return default_value;
    }

int inline OptionSet::
intVal(const Option& opt) const
    {
    return get(opt).intVal();
    }
int inline OptionSet::
intVal(const Option::Name& name) const
    {
    return get(name).intVal();
    }

int inline OptionSet::
intOrDefault(const Option::Name& name, int default_value) const
    {
    if(defined(name))
        return get(name).intVal();
    else
        return default_value;
    }

Real inline OptionSet::
realVal(const Option& opt) const
    {
    return get(opt).realVal();
    }
Real inline OptionSet::
realVal(const Option::Name& name) const
    {
    return get(name).realVal();
    }

Real inline OptionSet::
realOrDefault(const Option::Name& name, Real default_value) const
    {
    if(defined(name))
        return get(name).realVal();
    else
        return default_value;
    }

inline std::ostream& 
operator<<(std::ostream & s, const OptionSet& oset)
    {
    typedef OptionSet::const_iterator const_it;
    s << "/- OptionSet -------\n\n";
    for(const_it it = oset.opts_.begin();
        it != oset.opts_.end(); ++it)
        {
        s << it->second << "\n";
        }
    s << "\\------------------" << std::endl;
    return s;
    }

//
// Convenience functions for
// Options used within the library.
//

Option inline
Auto(bool val = true)
    {
    return Option("Auto",val);
    }

Option inline
Cutoff(int icut)
    {
    return Option("Cutoff",icut);
    }

Option inline
Cutoff(Real cut = 0)
    {
    return Option("Cutoff",cut);
    }

Option inline
DebugLevel(int level)
    {
    return Option("DebugLevel",level);
    }

Option inline
DoWrite(bool val = true)
    {
    return Option("DoWrite",val);
    }

Option inline
Pinning(Real val = 1)
    {
    return Option("Pinning",val);
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

Option inline
WriteM(int m)
    {
    return Option("WriteM",m);
    }

#endif
