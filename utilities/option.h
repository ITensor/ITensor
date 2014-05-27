//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_OPTION_H
#define __ITENSOR_OPTION_H

#include <vector>
#include "math.h"
#include "flstring.h"
#include "safebool.h"
#include <limits>

#ifndef NAN
#define NAN (std::numeric_limits<Real>::quiet_NaN())
#endif

namespace itensor {

typedef double Real;

class Opt : public safe_bool<Opt>
    {
    public:

    typedef FLString<20>
    Name;

    enum Type { Boolean, Numeric, String, None };

    Opt();

    Opt(const char* name);

    Opt(const Name& name);

    Opt(const Name& name, bool bval);

    Opt(const Name& name, const char* sval);
    Opt(const Name& name, const std::string& sval);

    Opt(const Name& name, int ival);

    Opt(const Name& name, Real rval);

    //
    // Accessor methods
    //

    const Name&
    name() const { return name_; }

    bool
    boolVal() const { assertType(Boolean); return bool(rval_); }

    const std::string&
    stringVal() const { assertType(String); return sval_; }

    int
    intVal() const { assertType(Numeric); return int(rval_); }

    Real
    realVal() const { assertType(Numeric); return rval_; }

    bool
    valid() const { return type_ != None; }

    Type
    type() const { return type_; }

    //operator const Name&() const { return name_; }

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

    Type type_;

    std::string sval_;
    Real rval_;

    //
    /////////////////////

    void
    assertType(Type t) const;

    };

//
// OptSet
//

class OptSet
    {
    public:

    typedef Opt::Name
    Name;

    typedef std::vector<Opt>
    storage_type;

    typedef storage_type::value_type
    value_type;

    typedef storage_type::iterator
    iterator;

    typedef storage_type::const_iterator
    const_iterator;

    OptSet();

    OptSet(const Opt& opt1);

    OptSet(const Opt& opt1, 
           const Opt& opt2, 
           const Opt& opt3 = Opt::Null(),
           const Opt& opt4 = Opt::Null());

    OptSet(const char* ostring);
    OptSet(const std::string& ostring);

    OptSet(const OptSet& other);

#ifdef USE_CPP11
    template <typename T, typename... Args>
    OptSet(const char* name1, 
           const T& t1, 
           const Args&... rest)
        {
        initialize(name1,t1,rest...);
        }

    template <typename... Args>
    OptSet(const OptSet& other,
           const Args&... rest)
        {
        initialize(other,rest...);
        }
#endif

    //
    // Methods for accessing Opts
    //

    bool
    defined(const Opt& opt) const;

    void
    add(const Opt& opt);
    void
    add(const Name& name, bool bval) { add(Opt(name,bval)); }
    void
    add(const Name& name, int ival) { add(Opt(name,ival)); }
    void
    add(const Name& name, const std::string& sval) { add(Opt(name,sval)); }
    void
    add(const Name& name, Real rval) { add(Opt(name,rval)); }

    void
    add(const Opt& opt1, 
        const Opt& opt2,
        const Opt& opt3 = Opt::Null(), 
        const Opt& opt4 = Opt::Null());

    void
    add(const char* ostring);

    const Opt&
    get(const Name& name) const;

    //
    // Methods for getting fields of a specific Opt
    //

    bool
    getBool(const Name& name) const;
    bool
    getBool(const Name& name, bool default_val) const;

    const std::string&
    getString(const Name& name) const;
    const std::string&
    getString(const Name& name, const std::string& default_val) const;

    int
    getInt(const Name& name) const;
    int
    getInt(const Name& name, int default_val) const;

    Real
    getReal(const Name& name) const;
    Real
    getReal(const Name& name, Real default_val) const;

    //
    // Iteration
    //

    iterator
    begin() { return opts_.begin(); }
    iterator
    end() { return opts_.end(); }

    const_iterator
    begin() const { return opts_.begin(); }
    const_iterator
    end() const { return opts_.end(); }

    const_iterator
    cbegin() const { return opts_.begin(); }
    const_iterator
    cend() const { return opts_.end(); }

    OptSet&
    operator+=(const OptSet& other);
    OptSet&
    operator&=(const OptSet& other) { return operator+=(other); }

    bool
    isGlobal() const { return (this == &GlobalOpts()); }

    static OptSet&
    GlobalOpts()
        {
        static OptSet gos_;
        return gos_;
        }

    private:

    ///////////////

    storage_type opts_;

    ///////////////

    OptSet(bool isGlobal);

    void
    processString(std::string ostring);

    void
    addByString(std::string ostring);

#ifdef USE_CPP11
    template <typename T, typename... Args>
    void
    initialize(const char* name1, 
               const T& t1, 
               const Args&... rest)
        {
        //std::cout << "Adding " << name1 << "=" << t1 << std::endl;
        add(Opt(name1,t1));
        initialize(rest...);
        }

    template <typename... Args>
    void
    initialize(const OptSet& other,
               const Args&... rest)
        {
        operator+=(other);
        initialize(rest...);
        }

    void
    initialize() { }
#endif

    };

OptSet
operator+(const Opt& opt1, const Opt& opt2);

OptSet
operator+(OptSet oset, const Opt& opt);

OptSet&
operator+=(OptSet& oset, const Opt& opt);

OptSet
operator+(OptSet oset, const OptSet& other);

OptSet
operator+(const Opt& opt, OptSet oset);

OptSet
operator+(const Opt& opt, const char* ostring);

OptSet
operator+(const char* ostring, const Opt& opt);

OptSet
operator+(OptSet oset, const char* ostring);

OptSet
operator+(const char* ostring, OptSet oset);

///////////

OptSet inline
operator&(const Opt& opt1, const Opt& opt2) { return opt1+opt2; }

OptSet inline
operator&(OptSet oset, const Opt& opt) { return oset + opt; }

inline 
OptSet&
operator&=(OptSet& oset, const Opt& opt) { return oset += opt; }

OptSet inline
operator&(OptSet oset, const OptSet& other) { return oset + other; }

OptSet inline
operator&(const Opt& opt, OptSet oset) { return opt + oset; }

OptSet inline
operator&(const Opt& opt, const char* ostring) { return opt + ostring; }

OptSet inline
operator&(const char* ostring, const Opt& opt) { return ostring + opt; }

OptSet inline
operator&(OptSet oset, const char* ostring) { return oset + ostring; }

OptSet inline
operator&(const char* ostring, OptSet oset) { return ostring + oset; }

std::ostream& 
operator<<(std::ostream & s, const Opt& opt);

std::ostream& 
operator<<(std::ostream & s, const OptSet& oset);

}; //namespace itensor

#endif
