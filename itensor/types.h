#ifndef __TYPES_H
#define __TYPES_H
#include <cmath>
#include <cstdlib>
#include "matrix.h"
#include <error.h> //utilities
#include <vector>
#include <iostream>
#include <fstream>
#include "assert.h"
#include "boost/foreach.hpp"
#define foreach BOOST_FOREACH
using namespace std::rel_ops;

//#define NMAX 8
static const int NMAX = 8;
static const Real MAX_CUT = 1E-15;
static const int MAX_M = 5000;

//----------------------------------
//For bounds checking - can remove once implementation is tested
//#define ITENSOR_USE_AT

#ifdef  ITENSOR_USE_AT
#define GET(container,j) (container.at(j))
#else
#define GET(container,j) (container[j])
#endif
//----------------------------------

#ifndef NDEBUG
#define DO_IF_DEBUG(X) X
#else
#define DO_IF_DEBUG(X)
#endif

//---------------------------------------
#define ENABLE_INTRUSIVE_PTR(ClassName) \
friend inline void intrusive_ptr_add_ref(ClassName* p) { ++(p->numref); } \
friend inline void intrusive_ptr_release(ClassName* p) { if(--(p->numref) == 0){ delete p; } } \
int count() const { return numref; }
//---------------------------------------

enum Printdat { ShowData, HideData };

#define Print(X) { printdat = false; cerr << "\n" << #X << " =\n" << X << "\n"; }
#define PrintDat(X) { printdat = true; cerr << "\n" << #X << " =\n" << X << "\n"; printdat = false; }


template<class T, class Op> void for_all(T& a, Op f) { for_each(a.begin(),a.end(),f); }

template<class T> std::vector<T>& operator*=(std::vector<T>& v1, const std::vector<T>& v2) 
{
    const unsigned int sz = v1.size();
    assert(v2.size() == sz);
    for(unsigned int n = 0; n < sz; ++n) v1[n] *= v2[n];
    return v1;
}
template<class T> std::vector<T> operator*(const std::vector<T>& v1, const std::vector<T>& v2) 
{ std::vector<T> res(v1); res *= v2; return res; }

template<class T> std::vector<T>& operator*=(std::vector<T>& v1, const std::vector<T*>& v2) 
{
    const unsigned int sz = v1.size();
    assert(v2.size() == sz);
    for(unsigned int n = 0; n < sz; ++n) v1[n] *= *(v2[n]);
    return v1;
}
template<class T> std::vector<T> operator*(const std::vector<T>& v1, const std::vector<T*>& v2) 
{ std::vector<T> res(v1); res *= v2; return res; }

template<class T> std::vector<T> operator*(const std::vector<const T*>& v1, const std::vector<const T*>& v2) 
{ 
    const size_t sz = v1.size();
    assert(v2.size() == sz);
    std::vector<T> res(sz); 
    for(size_t n = 0; n < sz; ++n) res[n] = *(v1[n]) * *(v2[n]);
    return res; 
}

template<class T>
ostream& operator<<(ostream& s, const std::vector<T>& v)
{ 
    if(v.size() == 0) s << "(Empty vector)\n";
    for(size_t n = 0; n < v.size(); ++n) { s << n << ": " << GET(v,n) << "\n"; } 
    return s; 
}

template<class T> T& operator*=(T& t1, const T* pt2) 
{ t1 *= *(pt2); return t1; }
template<class T> T operator*(const T& t1, const T* pt2) 
{ T res(t1); res *= *(pt2); return res; }

template<class T>
inline void readFromFile(const char* fname, T& t) 
    { std::ifstream s(fname); t.read(s); s.close(); }

template<class T>
inline void writeToFile(const char* fname, const T& t) 
    { std::ofstream s(fname); t.write(s); s.close(); }


extern bool printdat;
extern bool debug1, debug2, debug3, debug4;
extern Vector lastd;
Real ran1();
#ifdef THIS_IS_MAIN
void reportnew() {}
Real ran1(int);
bool printdat = false;
bool debug1 = false, debug2 = false, debug3 = false, debug4 = false;
Vector lastd(1);
#endif

#endif
