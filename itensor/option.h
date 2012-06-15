//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_OPTION_H
#define __ITENSOR_OPTION_H

#include <set>


class Option
    {
    public:

    //
    // List of all Option types,
    // in alphabetical order
    // 
    enum OptionType
        {
        NullOption,
        PreserveShape,
        Quiet,
        UseWF,
        Verbose
        };

    Option();

    Option(OptionType type);

    Option(OptionType type, bool bval);

    Option(OptionType type, const std::string& sval);

    Option(OptionType type, int ival);

    Option(OptionType type, Real rval);

    Option(OptionType type, 
           bool bval,
           const std::string& sval, 
           int ival, 
           Real rval);

    //
    // Operators for comparison and sorting
    //

    bool
    operator==(const Option& other) const
        { return type_ == other.type_; }

    bool
    operator!=(const Option& other) const
        { return type_ != other.type_; }

    bool
    operator<(const Option& other) const
        { return type_ < other.type_; }

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
    isNull() const { return type_ == NullOption; }
    bool
    isNotNull() const { return type_ != NullOption; }

    private:

    /////////////////////
    // Data Members

    OptionType type_;

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
    type_(NullOption),
    bval_(false),
    ival_(0),
    rval_(0)
    { }

inline Option::
Option(OptionType type)
    :
    type_(type),
    bval_(false),
    ival_(0),
    rval_(0)
    { }

inline Option::
Option(OptionType type, bool bval)
    :
    type_(type),
    bval_(bval),
    ival_(0),
    rval_(0)
    { }

inline Option::
Option(OptionType type, const std::string& sval)
    :
    type_(type),
    bval_(false),
    sval_(sval),
    ival_(0),
    rval_(0)
    { }

inline Option::
Option(OptionType type, int ival)
    :
    type_(type),
    bval_(false),
    ival_(ival),
    rval_(0)
    { }

inline Option::
Option(OptionType type, Real rval)
    :
    type_(type),
    bval_(false),
    ival_(0),
    rval_(rval)
    { }

inline Option::
Option(OptionType type, 
       bool bval,
       const std::string& sval, 
       int ival, 
       Real rval)
    :
    type_(type),
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
    includes(const Option& val) const { return opts_.count(val) > 0; }

    void
    insert(const Option& val) { if(val.isNotNull()) opts_.insert(val); }

    private:

    std::multiset<Option> opts_;

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


//
// Convenience functions for
// creating Option instances

Option inline
PreserveShape()
    {
    return Option(Option::PreserveShape);
    }

Option inline
Quiet()
    {
    return Option(Option::Quiet);
    }

Option inline
UseWF()
    {
    return Option(Option::UseWF);
    }

Option inline
Verbose()
    {
    return Option(Option::Verbose);
    }

#endif
