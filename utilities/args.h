//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_OPTION_H
#define __ITENSOR_OPTION_H

#include <vector>
#include "math.h"
#include <string>
#include "types.h"

namespace itensor {

class Args;
// Type names "OptSet" and "Opt" are
// aliases for Args for backwards compatibility.
using OptSet = Args;
using Opt = Args;

//
// Args
//

class Args
    {
    struct Val;
    public:

    using Name = std::string;
    using storage_type = std::vector<Val>;

    Args();

    //
    // Construct Args from a list of name-value pairs:
    // Args args("Name1",val1,"Name2",val2,"Name3",val3,...);
    //
    template <typename T, typename... Rest>
    Args(const char* name1, 
           const T& t1, 
           const Rest&... rest);

    //
    // Construct Args from another Args and a list of
    // name-value pairs:
    // Args args(other,"Name1",val1,"Name2",val2,"Name3",val3,...);
    //
    template <typename... Rest>
    Args(const Args& other,
         const Rest&... rest);


    //
    //Args("Name") is equivalent to Args("Name",true).
    //
    Args(const char* ostring);
    Args(const std::string& ostring);

    //
    // Copy and move constructors and assignment
    //
    Args(const Args& other);
    Args(Args&& other);
    Args&
    operator=(const Args& other);
    Args&
    operator=(Args&& other);


    //
    // Add a named value
    //
    void
    add(const Name& name, bool bval);
    void
    add(const Name& name, long ival);
    void
    add(const Name& name, int ival);
    void
    add(const Name& name, const char* sval);
    void
    add(const Name& name, const std::string& sval);
    void
    add(const Name& name, Real rval);
    void
    add(const char* ostring);

    //
    // Check if a specific name is defined in this Args instance
    //
    bool
    defined(const Name& name) const;

    //
    // Methods for getting values of named arguments
    //

    // Get value of bool-type argument, throws if not defined
    bool
    getBool(const Name& name) const;
    // Get value of bool-type argument, returns default_val if not defined
    bool
    getBool(const Name& name, bool default_val) const;

    // Get value of string-type argument, throws if not defined
    const std::string&
    getString(const Name& name) const;
    // Get value of string-type argument, returns default_val if not defined
    const std::string&
    getString(const Name& name, const std::string& default_val) const;

    // Get value of int-type argument, throws if not defined
    long
    getInt(const Name& name) const;
    // Get value of int-type argument, returns default_val if not defined
    long
    getInt(const Name& name, long default_val) const;

    // Get value of Real-type argument, throws if not defined
    Real
    getReal(const Name& name) const;
    // Get value of Real-type argument, returns default_val if not defined
    Real
    getReal(const Name& name, Real default_val) const;

    // Add contents of other to this
    Args&
    operator+=(const Args& other);

    // Check if this is the global Args object
    bool
    isGlobal() const { return (this == &Global()); }

    // Access the global Args object
    static Args&
    Global()
        {
        static Args gos_;
        return gos_;
        }

    private:

    ///////////////
    storage_type vals_;
    ///////////////

    void
    processString(std::string ostring);

    void
    addByString(std::string ostring);

    template <typename T, typename... Rest>
    void
    initialize(const char* name1, 
               const T& t1, 
               const Rest&... rest)
        {
        add(Val(name1,t1));
        initialize(rest...);
        }

    template <typename... Rest>
    void
    initialize(const Args& other,
               const Rest&... rest)
        {
        operator+=(other);
        initialize(rest...);
        }

    void
    initialize() { }

    void
    add(const Val& v);

    const Val&
    get(const Name& name) const;

    friend std::ostream& 
    operator<<(std::ostream & s, const Val& v);

    friend std::ostream& 
    operator<<(std::ostream & s, const Args& args);

    struct Val
        {
        enum Type { Boolean, Numeric, String, None };

        Val();

        Val(const char* name);

        Val(const Name& name);

        Val(const Name& name, bool bval);

        Val(const Name& name, const char* sval);
        Val(const Name& name, const std::string& sval);

        Val(const Name& name, long ival);
        Val(const Name& name, int ival);

        Val(const Name& name, Real rval);

        //
        // Accessor methods
        //

        const Name&
        name() const { return name_; }

        bool
        boolVal() const { assertType(Boolean); return bool(rval_); }

        const std::string&
        stringVal() const { assertType(String); return sval_; }

        long
        intVal() const { assertType(Numeric); return long(rval_); }

        Real
        realVal() const { assertType(Numeric); return rval_; }

        explicit operator bool() const { return valid(); }

        bool
        valid() const { return type_ != None; }

        Type
        type() const { return type_; }

        private:

        /////////////////
        Name name_;
        Type type_;
        std::string sval_;
        Real rval_;
        /////////////////

        void
        assertType(Type t) const;

        };

    public:
    
    //
    // Deprecated operator&= method. Use operator+= instead.
    //
    Args&
    operator&=(const Args& other) { return operator+=(other); }

    };


template <typename T, typename... Rest>
Args::
Args(const char* name1, 
       const T& t1, 
       const Rest&... rest)
    {
    initialize(name1,t1,rest...);
    }

template <typename... Rest>
Args::
Args(const Args& other,
     const Rest&... rest)
    {
    initialize(other,rest...);
    }

Args
operator+(Args args, const Args& other);

Args
operator+(Args args, const char* ostring);

Args
operator+(const char* ostring, Args args);




//
// Deprecated operator& methods. Use operator+ methods instead.
//

Args inline
operator&(Args args, const Args& other) { return args + other; }

Args inline
operator&(Args args, const char* ostring) { return args + ostring; }

Args inline
operator&(const char* ostring, Args args) { return ostring + args; }

}; //namespace itensor

#endif
