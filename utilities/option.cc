//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#include "option.h"
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

/*              +---------------------------+               
>---------------|     Overview and Use      |---------------<
                +---------------------------+               */

//Users create an Args class either through its constructor:
//       auto myargs = Args("Name1",value1,"Name2",value2,"Name3",value3,...);
//where value1, value2, value3, etc. can be any mixture of strings, bools, ints, or Reals.
//Or after creating an Args, users can add additional named arguments by using the "add" function:
//       myargs.add("Name4",value4);
//       myargs.add("Name5",value5);
//Users retrieve the value of a named argument by calling either:
//       myargs.getBool("Name4");
//or
//       myargs.getBool("Name4",default_value);
//
//Note:  Here we assume "Name4" corresponds to a boolean value, but others can be included.
//The first version throws an error if "Name4" is not a valid name for one of the named arguments 
//(which can be added to the Arg class). The second version either returns the value corresponding 
//to "Name4", or if it is not defined, returns default_value.

/*              +--------------------+               
>---------------|      Contents      |---------------<
                +--------------------+               */

//1. Val Constructor
//2. Overloaded << function
//3. Args Constructor
//4. Is Args variable defined?
//5. Add Args variable
//6. get name of Args variable
//7. getBool
//8. getString
//9. getInt
//10. getReal
//11. processString
//12. addByString
//13. Overloading +,+=

/*              +--------------------+               
>---------------| Val Constructors  |---------------<
                +--------------------+               */
// Construct an Args by providing its name and value. The value can be a bool, string (const char* or std::string), int, or Real.
// If no value is provided, the Args will be of bool type with its value set to true.

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

//If type is undefined, then an error is tossed.

void Args::Val::
assertType(Type t) const
    {
    if(t != type_)
        throw ITError("Wrong value type for option " + name_);
    }

/*              +------------------------+               
>---------------| Overloaded << function |---------------<
                +------------------------+               */

//Allows for Args to be printed with the << operator.

ostream& 
operator<<(ostream & s, const Args::Val& v)
    {
    s << v.name() << "=";
    if(v.type() == Args::Val::Boolean)
        {
        s << (v.boolVal() ? "true" : "false");//formatting for Booleans
        }
    else
    if(v.type() == Args::Val::Numeric)
        {
        s << v.realVal();//formatting for Real values
        }
    else
    if(v.type() == Args::Val::String)
        {
        s << "\"" << v.stringVal() << "\"";//formatting for Strings
        }
    else
        {
        s << "(Null)";//Prints Null for undefined types
        }
    return s;
    }

/*              +------------------------+               
>---------------|      Args function     |---------------<
                +------------------------+               */
//Defines Args class that stores pairs of strings and values.  
//
//Notes:
//"add" initialized in option.h and defined below
//was OptSet but now Opt and OptSet (set of Opts) are combined (v0,v1)

Args::
Args()
    { }

Args::
Args(const char* ostring)
    {
    processString(string(ostring));//single string
    }

Args::
Args(const string& ostring)
    {
    processString(ostring);//referenced string
    }

Args::
Args(const Args& other)
    { 
    if(!other.isGlobal())
        vals_ = other.vals_;//variable defined in Args now is global
    }

//C++11 allows for a slick way to copy variables at a low memory cost.
//Instead of copying a variable to another variable and then deleting the temporary memory storage,
//one can use the && operator to swap the pointers and leave the original variable with a Null pointer.
//This means that the other variable is simply 'forgotten' (without a memory leak).  
//Good for quickly copying large data types.

Args::
Args(Args&& other)
    { 
    if(!other.isGlobal())
        vals_ = std::move(other.vals_);//rvalue reference
    }

//Same as the above but with equals signs.
//
//ex:  someargs.add("nsweeps",20)
//         [Args something] = nsweeps

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

/*        +------------------+               
>---------|    Add to Args   |---------------<
          +------------------+               */
//Add a named argument having either boolean, integer, string, or Real value to Args.

void Args::
add(const Name& name, bool bval) { add({name,bval}); }
void Args::
add(const Name& name, int ival) { add({name,ival}); }
void Args::
add(const Name& name, const char* sval) { add({name,std::string(sval)}); }
void Args::
add(const Name& name, const std::string& sval) { add({name,sval}); }
void Args::
add(const Name& name, Real rval) { add({name,rval}); }

/*        +---------------------------------+               
>---------|    Is Args variable defined?     |---------------<
          +---------------------------------+               */
//Return true or false depnding on if specifed variable is available.


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

/*           +-------------------------+               
>------------|Add Val variable to Args |---------------<
             +-------------------------+               */
//Can add variable to Args just by variable.
//
//ex:  Int nsweeps = 20;
//         someargs.add(nsweeps)

void Args::
add(const Val& val)
    {
    if(!val) return;
    for(auto& x : vals_)//function automatically scans through names in Args to match the input
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

/*           +-------------------------------+               
>------------| get name of variable in Args  |---------------<
             +-------------------------------+               */

//Gets value of Val in Args class with name "name".
//
//ex:  get("maxm") returns maxm variable value 
//
//Notes:
//Initialized in option.h
//Useful for printing purposes without writing out the name by hand
 
const Args::Val& Args::
get(const Name& name) const
    {
    for(const auto& x : vals_)
        {
        if(x.name() == name) return x;//scans through Args variable types
        }
    //couldn't find the Val in this Args
    if(isGlobal())
        {
        throw ITError("Requested option " + name + " not found");
        }
    return Global().get(name);
    }

/*              +--------------------+               
>---------------|       getBool      |---------------<
                +--------------------+               */
//Get the value of the named argument ("name") assuming it is of a Boolean type.
//Throws an error if no argument with name "name" is found

bool Args::
getBool(const Name& name) const
    {
    return get(name).boolVal();
    }

//The second version accepts a default argument.

bool Args::
getBool(const Name& name, bool default_value) const
    {
    if(defined(name)) return get(name).boolVal();
    return default_value;
    }

/*              +--------------------+               
>---------------|       getString    |--------------<
                +--------------------+               */
//Get the value of the named argument ("name") assuming it is a string.
//Throws an error if no argument with name "name" is found
 
const string& Args::
getString(const Name& name) const
    {
    return get(name).stringVal();
    }

//The second version accepts a default argument.

const string& Args::
getString(const Name& name, const string& default_value) const
    {
    if(defined(name)) return get(name).stringVal();
    return default_value;
    }

/*              +--------------------+               
>---------------|       getInt       |---------------<
                +--------------------+               */
//Get the value of the named argument ("name") assuming it is an Integer type.
//Throws an error if no argument with name "name" is found

int Args::
getInt(const Name& name) const
    {
    return get(name).intVal();
    }

//The second version accepts a default argument.

int Args::
getInt(const Name& name, int default_value) const
    {
    if(defined(name)) return get(name).intVal();
    return default_value;
    }

/*              +--------------------+               
>---------------|       getReal      |--------------<
                +--------------------+               */
//Get the value of the named argument ("name") assuming it is a Real number type.
//Throws an error if no argument with name "name" is found

Real Args::
getReal(const Name& name) const
    {
    return get(name).realVal();
    }

//The second version accepts a default argument.

Real Args::
getReal(const Name& name, Real default_value) const
    {
    if(defined(name)) return get(name).realVal();
    return default_value;
    }

/*              +--------------------+               
>---------------|   processString    |--------------<
                +--------------------+               */
//Add variables from a string directly into an Args.  
//This will get the name and value of several variables if separated by commas.
//
//ex:  "maxm=15,minm=1,..."

void Args::
processString(string ostring)
    {
    ostring.erase(std::remove(ostring.begin(), ostring.end(),' '), ostring.end());//remove any spaces in the string (functions internal to C++)

    auto found = ostring.find_first_of(',');
    while(found != std::string::npos)
        {
        addByString(ostring.substr(0,found));//defined below
        ostring = ostring.substr(found+1);
        found = ostring.find_first_of(',');
        }

    addByString(ostring);
    }

/*              +--------------------+               
>---------------|    addByString     |--------------<
                +--------------------+               */
//addByString adds variables to an Args by taking the values from a string.
//
//Notes:
//used in processString (above)


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

/*              +--------------------+               
>---------------| Overloading +,+=   |--------------<
                +--------------------+               */

//Allows for the addition of variables to an Args with +,+=
//
//ex:  foo += x //where "x" is an object of type Args
//
//Running example:
//     Int nsweeps = 20
//         someArgs += nsweeps

Args& Args::
operator+=(const Args& args)
    {
    for(const auto& x : args.vals_)
        {
        add(x);
        }
    return *this;
    }

//----------OR------------
//        someotherArgs = someArgs + nsweeps

Args
operator+(Args args, const Args& other)
    {
    args += other;
    return args;
    }

//----------OR-------------
//        someArgs + "nsweeps=20"

Args 
operator+(Args args, const char* ostring)
    {
    args.add(ostring);
    return args;
    }

//----------OR-------------
// ????

Args 
operator+(const char* ostring, Args args)
    {
    args.add(ostring);
    return args;
    }
 
//----------OR-------------
//       cout << nsweep

ostream& 
operator<<(ostream & s, const Args& args)
    {
    if(args.isGlobal()) s << "Global Args:\n";
    else                s << "Args: (only showing overrides of global args)\n";

    for(const auto& opt : args.vals_)
        s << opt << "\n";

    return s;
    }

} //namespace itensor
