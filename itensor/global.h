//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_GLOBAL_H
#define __ITENSOR_GLOBAL_H
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <error.h> //utilities
#include "option.h"
#include "assert.h"
#include "boost/array.hpp"
#include "boost/format.hpp"

#include "boost/foreach.hpp"
#define Foreach BOOST_FOREACH

using namespace std::rel_ops;

static const int NMAX = 8;
static const Real MIN_CUT = 1E-20;
static const int MAX_M = 5000;


#ifndef DEBUG

#ifndef NDEBUG
#define NDEBUG //turn off asserts
#endif

#ifndef BOOST_DISABLE_ASSERTS
#define BOOST_DISABLE_ASSERTS
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
    bool savep = Global::printdat();\
    Global::printdat() = Y; \
    std::cout << "\n" << #X << " =\n" << X << std::endl; \
    Global::printdat() = savep;\
    }
#define Print(X)    PrintEither(X,false)
#define PrintDat(X) PrintEither(X,true)
#define PrintIndices(T) { T.printIndices(#T); }

/*
template<class T> std::vector<T>& 
operator*=(std::vector<T>& v1, const std::vector<T>& v2) 
    {
    const unsigned int sz = v1.size();
    assert(v2.size() == sz);
    for(unsigned int n = 0; n < sz; ++n) 
	v1[n] *= v2[n];
    return v1;
    }

template<class T> std::vector<T> 
operator*(const std::vector<T>& v1, const std::vector<T>& v2) 
    { 
    std::vector<T> res(v1); 
    res *= v2; 
    return res; 
    }

template<class T> std::vector<T>& 
operator*=(std::vector<T>& v1, const std::vector<T*>& v2) 
    {
    const unsigned int sz = v1.size();
    assert(v2.size() == sz);
    for(unsigned int n = 0; n < sz; ++n) 
	v1[n] *= *(v2[n]);
    return v1;
    }

template<class T> std::vector<T> 
operator*(const std::vector<T>& v1, const std::vector<T*>& v2) 
    { 
    std::vector<T> res(v1); 
    res *= v2; 
    return res; 
    }

template<class T> std::vector<T> 
operator*(const std::vector<const T*>& v1, const std::vector<const T*>& v2) 
    { 
    const size_t sz = v1.size();
    assert(v2.size() == sz);
    std::vector<T> res(sz); 
    for(size_t n = 0; n < sz; ++n) 
	res[n] = *(v1[n]) * *(v2[n]);
    return res; 
    }

template<class T> std::ostream& 
operator<<(std::ostream& s, const std::vector<T>& v)
    { 
    if(v.size() == 0) 
	s << "(Empty vector)\n";
    for(size_t n = 0; n < v.size(); ++n) 
	s << n << ": " << GET(v,n) << "\n";
    return s; 
    }

template<class T> T& 
operator*=(T& t1, const T* pt2) 
    { 
    t1 *= *(pt2); 
    return t1; 
    }

template<class T> T 
operator*(const T& t1, const T* pt2) 
    { 
    T res(t1); 
    res *= *(pt2); 
    return res; 
    }
*/

bool inline
fileExists(const std::string& fname)
    {
    std::ifstream file(fname.c_str());
    return file.good();
    }
bool inline
fileExists(const boost::format& fname)
    {
    return fileExists(fname.str());
    }

template<class T> 
void inline
readFromFile(const std::string& fname, T& t) 
    { 
    std::ifstream s(fname.c_str()); 
    if(!s.good()) 
        Error("Couldn't open file \"" + fname + "\" for reading");
    t.read(s); 
    s.close(); 
    }

template<class T> 
void inline
readFromFile(const boost::format& fname, T& t) 
    { 
    readFromFile(fname.str(),t);
    }

template<class T> 
void inline
writeToFile(const std::string& fname, const T& t) 
    { 
    std::ofstream s(fname.c_str()); 
    if(!s.good()) 
        Error("Couldn't open file \"" + fname + "\" for writing");
    t.write(s); 
    s.close(); 
    }

template<class T> 
void inline
writeToFile(const boost::format& fname, const T& t) 
    { 
    writeToFile(fname.str(),t); 
    }

void inline
writeVec(std::ostream& s, const Vector& V)
    {
    int m = V.Length();
    s.write((char*)&m,sizeof(m));
    Real val;
    for(int k = 1; k <= m; ++k)
        {
        val = V(k);
        s.write((char*)&val,sizeof(val));
        }
    }

void inline
readVec(std::istream& s, Vector& V)
    {
    int m = 1;
    s.read((char*)&m,sizeof(m));
    V.ReDimension(m);
    Real val;
    for(int k = 1; k <= m; ++k)
        {
        s.read((char*)&val,sizeof(val));
        V(k) = val;
        }
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
    char* cstr;
    cstr = new char[dirname.size()+1];
    strcpy(cstr,dirname.c_str());

    //Call mkdtemp
    char* retval = mkdtemp(cstr);
    //Check error condition
    if(retval == NULL)
        {
        delete[] cstr;
        throw ITError("mkTempDir failed");
        }

    //Prepare return value
    std::string final_dirname(retval);

    //Clean up
    delete[] cstr;

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

enum Arrow { In = -1, Out = 1 };

Arrow inline
operator*(const Arrow& a, const Arrow& b)
    { 
    return (int(a)*int(b) == int(In)) ? In : Out; 
    }
const Arrow Switch = In*Out;

inline std::ostream& 
operator<<(std::ostream& s, Arrow D)
    { 
    s << (D == In ? "In" : "Out");
    return s; 
    }

////////
///////


class Global
    {
    public:

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
    static Vector& 
    lastd()
        {
        static Vector lastd_(1);
        return lastd_;
        }
    static bool& 
    checkArrows()
        {
        static bool checkArrows_ = true;
        return checkArrows_;
        }
    static OptionSet&
    options()
        {
        static OptionSet oset_;
        return oset_;
        }
    };


class ResultIsZero : public ITError
    {
    public:

    typedef ITError
    Parent;

    ResultIsZero(const std::string& message) 
        : Parent(message)
        { }
    };

class ArrowError : public ITError
    {
    public:

    typedef ITError
    Parent;

    ArrowError(const std::string& message) 
        : Parent(message)
        { }
    };

Real ran1();

//void reportnew() { }

#endif
