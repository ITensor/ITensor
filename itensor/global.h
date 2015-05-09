//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_GLOBAL_H
#define __ITENSOR_GLOBAL_H
#include <cmath>
#include <cstdlib>
#include <sys/types.h>
#include <unistd.h>
#include <fstream>
#include <complex>
#include "assert.h"
#include "error.h"
#include "args.h"
#include "types.h"
#include <ctime>
#include <string.h>
#include <cstring>
#include "real.h"

namespace itensor {

enum Direction { Fromright, Fromleft, Both, None };

static const Real MIN_CUT = 1E-15;
static const int MAX_M = 5000;


static const Complex Complex_1 = Complex(1,0);
static const Complex Complex_i = Complex(0,1);

// The PAUSE macro is useful for debugging. 
// Prints the current line number and pauses
// execution until the enter key is pressed.
#define PAUSE { std::cout << "(Paused, Line " << __LINE__ << ")"; std::cin.get(); }


#ifndef DEBUG
#ifndef NDEBUG
#define NDEBUG //turn off asserts
#endif
#endif

#ifdef DEBUG
#define DO_IF_DEBUG(X) X
#else
#define DO_IF_DEBUG(X)
#endif


#ifdef DEBUG
#define GET(container,j) (container.at(j))
#else
#define GET(container,j) (container[j])
#endif	

enum Printdat { ShowData, HideData };

#define PrintEither(X,Y) \
    {\
    const bool savep = Global::printdat();\
    Global::printdat() = Y; \
    std::cout << "\n" << #X << " =\n" << X << std::endl; \
    Global::printdat() = savep;\
    }
#define Print(X)    PrintEither(X,false)
#define PrintDat(X) PrintEither(X,true)
#define PrintData(X) PrintEither(X,true)


bool inline
fileExists(const std::string& fname)
    {
    std::ifstream file(fname.c_str());
    return file.good();
    }

// Can overload read and write functions with
// signatures below for objects such as std::vector
// 
// For classes having member read/write functions, can 
// leave external read/write overloads undefined
// and the following template overloads will 
// be called
//

//Here we have to use a struct to implement the read(istream,T)
//function because function templates cannot be partially specialized
template<typename T, bool isPod = std::is_pod<T>::value>
struct DoRead
    {
    DoRead(std::istream& s, T& obj)
        {
        obj.read(s);
        }
    };
template<typename T>
struct DoRead<T, true>
    {
    DoRead(std::istream& s, T& val)
        {
        s.read((char*) &val, sizeof(val));
        }
    };
template<typename T>
void
read(std::istream& s, T& val)
    {
    DoRead<T>(s,val);
    }

template<typename T, typename... CtrArgs>
T
read(std::istream& s, CtrArgs&&... args)
    {
    T t(std::forward<CtrArgs>(args)...);
    DoRead<T>(s,t);
    return t;
    }

template<typename T, bool isPod = std::is_pod<T>::value>
struct DoWrite
    {
    DoWrite(std::ostream& s, const T& obj)
        {
        obj.write(s);
        }
    };
template<typename T>
struct DoWrite<T, true>
    {
    DoWrite(std::ostream& s, const T& val)
        {
        s.write((char*) &val, sizeof(val));
        }
    };
template<typename T>
void
write(std::ostream& s, const T& val)
    {
    DoWrite<T>(s,val);
    }

template<typename T>
void
write(std::ostream& s, const std::vector<T>& vec)
    {
    auto size = vec.size();
    s.write((char*)&size,sizeof(size));
    s.write((char*)vec.data(),sizeof(T)*size);
    }

template<typename T>
void
read(std::istream& s, std::vector<T>& vec)
    {
    auto size = vec.size(); //will overwrite
    s.read((char*)&size,sizeof(size));
    vec.resize(size);
    s.read((char*)vec.data(),sizeof(T)*size);
    }

void inline
read(std::istream& s, Cplx& z)
    {
    auto &r = reinterpret_cast<Real(&)[2]>(z)[0];
    auto &i = reinterpret_cast<Real(&)[2]>(z)[1];
    s.read((char*)&r,sizeof(r));
    s.read((char*)&i,sizeof(i));
    }

void inline
write(std::ostream& s, const Cplx& z)
    {
    auto &r = reinterpret_cast<const Real(&)[2]>(z)[0];
    auto &i = reinterpret_cast<const Real(&)[2]>(z)[1];
    s.write((char*)&r,sizeof(r));
    s.write((char*)&i,sizeof(i));
    }

template<class T> 
void
readFromFile(const std::string& fname, T& t) 
    { 
    std::ifstream s(fname.c_str(),std::ios::binary);
    if(!s.good()) 
        Error("Couldn't open file \"" + fname + "\" for reading");
    read(s,t); 
    s.close(); 
    }


template<class T, typename... InitArgs>
T
readFromFile(const std::string& fname, InitArgs&&... iargs)
    { 
    std::ifstream s(fname.c_str(),std::ios::binary); 
    if(!s.good()) 
        Error("Couldn't open file \"" + fname + "\" for reading");
    T t(std::forward<InitArgs>(iargs)...);
    read(s,t); 
    s.close(); 
    return t;
    }


template<class T> 
void
writeToFile(const std::string& fname, const T& t) 
    { 
    std::ofstream s(fname.c_str(),std::ios::binary); 
    if(!s.good()) 
        Error("Couldn't open file \"" + fname + "\" for writing");
    write(s,t); 
    s.close(); 
    }

//Given a prefix (e.g. pfix == "mydir")
//and an optional location (e.g. locn == "/var/tmp/")
//creates a temporary directory and returns its name
//without a trailing slash
//(e.g. /var/tmp/mydir_SfqPyR)
std::string inline
mkTempDir(const std::string& pfix,
          const std::string& locn = "./")
    {
    //Construct dirname
    std::string dirname = locn;
    if(dirname[dirname.length()-1] != '/')
        dirname += '/';
    //Add prefix and template string of X's for mkdtemp
    dirname += pfix + "_XXXXXX";

    //Create C string version of dirname
    auto cstr = std::unique_ptr<char[]>(new char[dirname.size()+1]);
    strcpy(cstr.get(),dirname.c_str());

    //Call mkdtemp
    char* retval = mkdtemp(cstr.get());
    //Check error condition
    if(retval == NULL) throw ITError("mkTempDir failed");

    //Prepare return value
    std::string final_dirname(retval);

    return final_dirname;
    }


/*
*
* The Arrow enum is used to label how indices
* transform under a particular symmetry group. 
* Indices with an Out Arrow transform as vectors
* (kets) and with an In Arrow as dual vectors (bras).
*
* Conventions regarding arrows:
*
* * Arrows point In or Out, never right/left/up/down.
*
* * The Site indices of an MPS representing a ket point Out.
*
* * Conjugation switches arrow directions.
*
* * All arrows flow Out from the ortho center of an MPS 
*   (assuming it's a ket - In if it's a bra).
*
* * IQMPOs are created with the same arrow structure as if they are 
*   orthogonalized to site 1, but this is just a default since they 
*   aren't actually ortho. If position is called on an IQMPO it follows 
*   the same convention as for an MPS except Site indices point In and 
*   Site' indices point Out.
*
* * Local site operators have two IQIndices, one unprimed and pointing In, 
*   the other primed and pointing Out.
*
*/

enum Arrow { In = -1, Out = 1, Neither = 0 };

Arrow inline
operator-(Arrow dir)
    {
#ifdef DEBUG
    if(dir == Neither)
        Error("Cannot reverse Arrow direction 'Neither'");
#endif
    return (dir == In ? Out : In);
    }

inline std::ostream& 
operator<<(std::ostream& s, Arrow D)
    { 
    switch(D)
        {
        case In:
            s << "In";
            return s;
        case Out:
            s << "Out";
            return s;
        case Neither:
            s << "Neither";
            return s;
        default:
            Error("Missing Arrow case");
        }
    return s; 
    }

////////
///////


class Global
    {
    public:
    static bool& 
    checkArrows()
        {
        static bool checkArrows_ = true;
        return checkArrows_;
        }
    static bool& 
    debug1()
        {
        static bool debug1_ = false;
        return debug1_;
        }
    static bool& 
    debug2()
        {
        static bool debug2_ = false;
        return debug2_;
        }
    static bool& 
    debug3()
        {
        static bool debug3_ = false;
        return debug3_;
        }
    static bool& 
    debug4()
        {
        static bool debug4_ = false;
        return debug4_;
        }
    //Global named arguments
    static Args&
    args()
        {
        return Args::Global();
        }
    void static
    args(const Args::Name& name, bool bval)
        {
        Args::Global().add(name,bval);
        }
    void static
    args(const Args::Name& name, int ival)
        {
        Args::Global().add(name,ival);
        }
    void static
    args(const Args::Name& name, Real rval)
        {
        Args::Global().add(name,rval);
        }
    void static
    args(const Args::Name& name, const std::string& sval)
        {
        Args::Global().add(name,sval);
        }
    static bool& 
    printdat()
        {
        static bool printdat_ = false;
        return printdat_;
        }
    static Real& 
    printScale()
        {
        static Real printScale_ = 1E-10;
        return printScale_;
        }
    static Real
    random()//int seed = 0)
        {
        int seed = 0;
        using Generator = mt19937;
        using Distribution = uniform_real_distribution<Real>;

        static Generator rng(std::time(NULL)+getpid());
        static Distribution dist(0,1);

        if(seed != 0)  //reseed rng
            {
            rng = Generator(seed);
            }

        return dist(rng);
        }
    void static
    warnDeprecated(const std::string& message)
        {
        static int depcount = 1;
        if(depcount <= 10)
            {
            println("\n\n",message,"\n");
            ++depcount;
            }
        }
    };


class ResultIsZero : public ITError
    {
    public:

    using Parent = ITError;

    ResultIsZero(const std::string& message) 
        : Parent(message)
        { }
    };

class ArrowError : public ITError
    {
    public:

    using Parent = ITError;

    ArrowError(const std::string& message) 
        : Parent(message)
        { }
    };


} //namespace itensor

#endif
