//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_OPTION_H
#define __ITENSOR_OPTION_H

#include <vector>
#include <string>
#include "math.h"
#include "itensor/types.h"
#include "itensor/util/infarray.h"

namespace itensor {

//
// Args - named argument system
//
// o An Args object holds a collection of name-value
//   pairs. To access an integer value for example,
//   use args.getInt("Name1"); or call
//   args.getInt("Name1",def); to provide a 
//   default value def.
// o To add a value do args.add("Name2",val);
// o To construct an Args with a given set of values
//   do Args args("Name1",val1,"Name2",val2,...);
// o There is a global Args, Args::global(). Values
//   not present in a given args instance will be
//   looked up in the global Args object before
//   either the default is selected or an error
//   thrown if no default is provided.
// o To have a function accept Args in read-only
//   mode, use the signature 
//   func(T1 t1, T2 t2, ..., const Args& args = Args::global());
//   which will incur essentially no overhead.
//   If you intend to add or modify the args set, take it by value.
//

class Args
    {
    class Val;
    public:
    using Name = std::string;
    using storage_type = InfArray<Val,7ul>;

    Args();

    //
    // Construct Args from a list of name-value pairs:
    // Args args("Name1",val1,"Name2",val2,"Name3",val3,...);
    //
    template <typename T, typename... Rest>
    Args(const char* name1, 
         T const& t1, 
         Rest const&... rest);

    //
    // Construct Args from another Args and a list of
    // name-value pairs:
    // Args args(other,"Name1",val1,"Name2",val2,"Name3",val3,...);
    //
    template <typename... Rest>
    Args(Args const& other,
         Rest const&... rest);


    //
    //Args("Name") is equivalent to Args("Name",true).
    //
    Args(const char* ostring);

    Args(std::string const& ostring);

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
    add(Name const& name, bool bval);
    void     
    add(Name const& name, long ival);
    void     
    add(Name const& name, int ival);
    void     
    add(Name const& name, const char* sval);
    void     
    add(Name const& name, std::string const& sval);
    void     
    add(Name const& name, Real rval);
    void
    add(const char* ostring);

    // Check if a specific name is defined in this Args instance
    bool
    defined(Name const& name) const;

    // Remove an arg from the set - always succeeds
    void
    remove(Name const& name);

    //
    // Methods for getting values of named arguments
    //

    // Get value of bool-type argument, throws if not defined
    bool
    getBool(Name const& name) const;
    // Get value of bool-type argument, returns default_val if not defined
    bool
    getBool(Name const& name, bool default_val) const;

    // Get value of string-type argument, throws if not defined
    const std::string&
    getString(Name const& name) const;
    // Get value of string-type argument, returns default_val if not defined
    const std::string&
    getString(Name const& name, std::string const& default_val) const;

    // Get value of int-type argument, throws if not defined
    long
    getInt(Name const& name) const;
    // Get value of int-type argument, returns default_val if not defined
    long
    getInt(Name const& name, long default_val) const;

    // Get value of int-type argument, throws if not defined
    size_t
    getSizeT(Name const& name) const;
    // Get value of int-type argument, returns default_val if not defined
    size_t
    getSizeT(Name const& name, long default_val) const;

    // Get value of Real-type argument, throws if not defined
    Real
    getReal(Name const& name) const;
    // Get value of Real-type argument, returns default_val if not defined
    Real
    getReal(Name const& name, Real default_val) const;

    // Add contents of other to this
    Args&
    operator+=(Args const& other);

    // Check if this is the global Args object
    bool
    isGlobal() const { return (this == &global()); }

    // Access the global Args object
    static Args&
    global()
        {
        static Args gos_;
        return gos_;
        }

    // Read Args object from binary
    void
    read(std::istream& s);

    // Write Args object to binary
    void
    write(std::ostream& s) const;

    private:

    void
    processString(std::string ostring);

    void
    addByString(std::string ostring);

    template <typename T, typename... Rest>
    void
    initialize(const char* name1, 
               T const& t1, 
               Rest const&... rest)
        {
        add(Val(name1,t1));
        initialize(rest...);
        }

    template <typename... Rest>
    void
    initialize(Args const& other,
               Rest const&... rest)
        {
        operator+=(other);
        initialize(rest...);
        }

    void
    initialize() { }

    void
    add(Val const& v);

    Val const&
    get(Name const& name) const;

    friend std::ostream& 
    operator<<(std::ostream & s, Val const& v);

    friend std::ostream& 
    operator<<(std::ostream & s, Args const& args);

    class Val
        {
        public:
        enum Type { Boolean, Numeric, String, None };
        private:
        Name name_;
        Type type_;
        std::string sval_;
        Real rval_;
        public:


        Val();

        Val(const char* name);

        Val(Name const& name);
                      
        Val(Name const& name, bool bval);
                      
        Val(Name const& name, const char* sval);
        Val(Name const& name, std::string const& sval);
                      
        Val(Name const& name, long ival);
        Val(Name const& name, int ival);
        Val(Name const& name, unsigned long ival);
        Val(Name const& name, unsigned int ival);
                      
        Val(Name const& name, Real rval);

        //
        // Accessor methods
        //

        Name const&
        name() const { return name_; }

        bool
        boolVal() const { assertType(Boolean); return bool(rval_); }

        std::string const&
        stringVal() const { assertType(String); return sval_; }

        long
        intVal() const { assertType(Numeric); return long(rval_); }

        size_t
        size_tVal() const { assertType(Numeric); return size_t(rval_); }

        Real
        realVal() const { assertType(Numeric); return rval_; }

        explicit operator bool() const { return valid(); }

        bool
        valid() const { return type_ != None; }

        Type
        type() const { return type_; }

        void
        read(std::istream& s);

        void
        write(std::ostream& s) const;
        
        private:

        void
        assertType(Type t) const;

        };

    storage_type vals_;

    };


template <typename T, typename... Rest>
Args::
Args(const char* name1, 
       T const& t1, 
       Rest const&... rest)
    {
    initialize(name1,t1,rest...);
    }

template <typename... Rest>
Args::
Args(Args const& other,
     Rest const&... rest)
    {
    initialize(other,rest...);
    }

Args
operator+(Args args, Args const& other);

Args
operator+(Args args, const char* ostring);

Args
operator+(const char* ostring, Args args);

} //namespace itensor

#endif
