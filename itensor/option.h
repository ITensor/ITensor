//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_OPTION_H
#define __ITENSOR_OPTION_H

#include <map>
#include "real.h"

class Opt
    {
    public:

    typedef std::string
    Name;

    Opt();

    explicit
    Opt(const Name& name);

    Opt(const Name& name, bool bval);

    Opt(const Name& name, const std::string& sval);

    Opt(const Name& name, int ival);

    Opt(const Name& name, Real rval);

    Opt(const Name& name, 
           bool bval,
           const std::string& sval, 
           int ival, 
           Real rval);

    //
    // Operators for comparison and sorting
    //

    // Two Opts are equal if they have the same name,
    // regardless of other fields that may be set.
    bool
    operator==(const Opt& other) const
        { 
        return name_ == other.name_;
        }

    // Compares two Opts based on their name, using
    // the < operator of std::string. Useful for sorting.
    bool
    operator<(const Opt& other) const
        { 
        return name_ < other.name_;
        }

    bool
    boolEquals(const Opt& other) const
        {
        return name_ == other.name_
            && bval_ == other.bval_;
        }

    bool
    intEquals(const Opt& other) const
        {
        return name_ == other.name_
            && ival_ == other.ival_;
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
    isNull() const { return name_ == "NullOpt"; }
    bool
    isNotNull() const { return name_ != "NullOpt"; }

    const Name&
    name() const { return name_; }

    //operator const Name&() const { return name_; }

    friend std::ostream& 
    operator<<(std::ostream & s, const Opt& opt);

    static Opt&
    Null()
        {
        static Opt null_;
        return null_;
        }

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

//
// OptSet
//

class OptSet
    {
    public:

    typedef Opt::Name
    Name;

    OptSet();

    OptSet(const Opt& opt1);

    OptSet(const Opt& opt1, 
           const Opt& opt2, 
           const Opt& opt3 = Opt::Null(),
           const Opt& opt4 = Opt::Null());


    OptSet(const OptSet& other);

    //
    // Methods for accessing Opts
    //

    bool
    defined(const Name& name) const;
    bool
    defined(const Opt& opt) const;

    void
    add(const Opt& opt) { if(opt.isNotNull()) opts_[opt.name()] = opt; }

    void
    add(const Opt& opt1, 
        const Opt& opt2,
        const Opt& opt3 = Opt::Null(), 
        const Opt& opt4 = Opt::Null());

    const Opt&
    get(const Name& name) const;
    const Opt&
    get(const Opt& opt) const;

    //
    // Methods for getting fields of a specific Opt
    //

    bool
    getBool(const Name& name) const;
    bool
    getBool(const Opt& opt) const;
    bool
    getBool(const Name& name, bool default_val) const;

    const std::string&
    getString(const Name& name) const;
    const std::string&
    getString(const Opt& opt) const;
    const std::string&
    getString(const Name& name, const std::string& default_val) const;

    int
    getInt(const Name& name) const;
    int
    getInt(const Opt& opt) const;
    int
    getInt(const Name& name, int default_val) const;

    Real
    getReal(const Name& name) const;
    Real
    getReal(const Opt& opt) const;
    Real
    getReal(const Name& name, Real default_val) const;

    friend std::ostream& 
    operator<<(std::ostream & s, const OptSet& oset);

    static OptSet&
    GlobalOpts()
        {
        const bool isGlobal = true;
        static OptSet gos_(isGlobal);
        return gos_;
        }

    private:

    typedef std::map<Name,Opt>
    storage_type;

    typedef storage_type::iterator
    iterator;

    typedef storage_type::const_iterator
    const_iterator;

    ///////////////
    //
    // Data Members

    storage_type opts_;

    const bool is_global_;

    //
    ///////////////

    OptSet(bool isGlobal);

    };

OptSet inline
operator&(const Opt& opt1, const Opt& opt2)
    {
    return OptSet(opt1,opt2);
    }

OptSet inline
operator&(OptSet oset, const Opt& opt)
    {
    oset.add(opt);
    return oset;
    }

OptSet inline
operator&(const Opt& opt, OptSet oset)
    {
    oset.add(opt);
    return oset;
    }

inline Opt::
Opt()
    :
    name_("NullOpt"),
    bval_(false),
    ival_(-1E5),
    rval_(NAN)
    { }

inline Opt::
Opt(const Name& name)
    :
    name_(name),
    bval_(true),
    ival_(-1E5),
    rval_(NAN)
    { }

inline Opt::
Opt(const Name& name, bool bval)
    :
    name_(name),
    bval_(bval),
    ival_(-1E5),
    rval_(NAN)
    { }

inline Opt::
Opt(const Name& name, const std::string& sval)
    :
    name_(name),
    bval_(true),
    sval_(sval),
    ival_(-1E5),
    rval_(NAN)
    { }

inline Opt::
Opt(const Name& name, int ival)
    :
    name_(name),
    bval_(true),
    ival_(ival),
    rval_(NAN)
    { }

inline Opt::
Opt(const Name& name, Real rval)
    :
    name_(name),
    bval_(true),
    ival_(-1E5),
    rval_(rval)
    { }

inline Opt::
Opt(const Name& name, 
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
operator<<(std::ostream & s, const Opt& opt)
    {
    s << "Opt \"" << opt.name() << "\"\n";
    s << "  boolVal   = " << (opt.boolVal() ? "true" : "false") << "\n";
    s << "  intVal    = " << opt.intVal() << "\n";
    s << "  realVal   = " << opt.realVal() << "\n";
    s << "  stringVal = \"" << opt.stringVal() << "\"" 
      << std::endl;
    return s;
    }



inline OptSet::
OptSet()
    :
    is_global_(false)
    { }

inline OptSet::
OptSet(const Opt& opt)
    :
    is_global_(false)
    {
    add(opt);
    }

inline OptSet::
OptSet(const Opt& opt1, const Opt& opt2, 
       const Opt& opt3, const Opt& opt4)
    :
    is_global_(false)
    {
    add(opt1,opt2,opt3,opt4);
    }

inline OptSet::
OptSet(bool isGlobal)
    :
    is_global_(isGlobal)
    { }

inline OptSet::
OptSet(const OptSet& other)
    :
    is_global_(false)
    { 
    if(!other.is_global_)
        opts_ = other.opts_;
    }

bool inline OptSet::
defined(const Opt::Name& name) const
    {
    if(opts_.count(name) > 0)
        return true;

    if(is_global_) 
        return false;
    //else see if GlobalOpts contains it
    return GlobalOpts().defined(name);
    }

bool inline OptSet::
defined(const Opt& opt) const
    {
    return defined(opt.name());
    }

void inline OptSet::
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

inline 
const Opt& OptSet::
get(const Opt& opt) const
    {
    return get(opt.name());
    }

inline 
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

bool inline OptSet::
getBool(const Opt& opt) const
    {
    return get(opt).boolVal();
    }
bool inline OptSet::
getBool(const Opt::Name& name) const
    {
    return get(name).boolVal();
    }

bool inline OptSet::
getBool(const Opt::Name& name, bool default_value) const
    {
    if(defined(name))
        return get(name).boolVal();
    else
        return default_value;
    }

inline 
const std::string& OptSet::
getString(const Opt& opt) const
    {
    return get(opt).stringVal();
    }
inline 
const std::string& OptSet::
getString(const Opt::Name& name) const
    {
    return get(name).stringVal();
    }

inline const std::string& OptSet::
getString(const Opt::Name& name, const std::string& default_value) const
    {
    if(defined(name))
        return get(name).stringVal();
    else
        return default_value;
    }

int inline OptSet::
getInt(const Opt& opt) const
    {
    return get(opt).intVal();
    }
int inline OptSet::
getInt(const Opt::Name& name) const
    {
    return get(name).intVal();
    }

int inline OptSet::
getInt(const Opt::Name& name, int default_value) const
    {
    if(defined(name))
        return get(name).intVal();
    else
        return default_value;
    }

Real inline OptSet::
getReal(const Opt& opt) const
    {
    return get(opt).realVal();
    }
Real inline OptSet::
getReal(const Opt::Name& name) const
    {
    return get(name).realVal();
    }

Real inline OptSet::
getReal(const Opt::Name& name, Real default_value) const
    {
    if(defined(name))
        return get(name).realVal();
    else
        return default_value;
    }

inline 
std::ostream& 
operator<<(std::ostream & s, const OptSet& oset)
    {
    typedef OptSet::const_iterator const_it;

    if(oset.is_global_)
        s << "/- Global OptSet -------\n\n";
    else
        s << "/- OptSet (only showing overrides of global opts) -------\n\n";

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
// Opts used within the library.
//

Opt inline
Auto(bool val = true)
    {
    return Opt("Auto",val);
    }

Opt inline
ConserveNf(bool val = true)
    {
    return Opt("ConserveNf",val);
    }

Opt inline
Cutoff(int icut)
    {
    return Opt("Cutoff",icut);
    }

Opt inline
Cutoff(Real cut = 0)
    {
    return Opt("Cutoff",cut);
    }

Opt inline
DebugLevel(int level)
    {
    return Opt("DebugLevel",level);
    }

Opt inline
DoNormalize(bool val = true)
    {
    return Opt("DoNormalize",val);
    }

Opt inline
Maxm(int val)
    {
    return Opt("Maxm",val);
    }

Opt inline
Minm(int val)
    {
    return Opt("Minm",val);
    }

Opt inline
NumCenter(int nc = 2)
    {
    return Opt("NumCenter",nc);
    }

Opt inline
Offset(int n = 0)
    {
    return Opt("Offset",n);
    }

Opt inline
Pinning(Real val = 1)
    {
    return Opt("Pinning",val);
    }


Opt inline
PreserveShape()
    {
    return Opt("PreserveShape");
    }

Opt inline
Quiet(bool val = true)
    {
    return Opt("Quiet",val);
    }

Opt inline
UseWF()
    {
    return Opt("UseWF");
    }

Opt inline
Verbose(bool val = true)
    {
    return Opt("Verbose",val);
    }

Opt inline
Weight(Real w = 1)
    {
    return Opt("Weight",w);
    }

Opt inline
WriteDir(const std::string& dirname)
    {
    if(dirname[dirname.length()-1] == '/')
        return Opt("WriteDir",dirname);
    else
        return Opt("WriteDir",dirname + "/");
    }

Opt inline
WriteM(int m)
    {
    return Opt("WriteM",m);
    }

#endif
