// tarray1.h -- template class for array with one index, starting at offset. 

#ifndef _tarray1_h
#define _tarray1_h
#define SP << " " <<

#include "error.h"

#include <iostream>
#include "indent.h"

//---------------------------
#define HEADER_DEFS
//---------------------------

//#define BOUNDS	        /* Define this if you want bounds checking */

typedef void copyarrayfun( void *, const void *, int);
typedef void copyvaluefun( void *, const void *, int);
typedef void* donewfun(int); 
typedef void dodeletefun(void *);
typedef void outputfun(std::ostream&,void *,int,int);

class FunPoint
    {
public:
    copyarrayfun * copyarray;
    copyvaluefun * copyvalue;
    donewfun * donew;
    dodeletefun * dodelete;
    outputfun * outputarray;
    FunPoint(copyarrayfun * a, copyvaluefun * b, donewfun * c, 
	dodeletefun * d, outputfun * e) : copyarray(a), copyvalue(b), 
	donew(c), dodelete(d), outputarray(e) { }
    };

enum {SPECIAL = -398777};

class Array1_base
    {
    Array1_base(const Array1_base&);
    void operator=(const Array1_base&);
public:
    int size;
    int storesize;
    short offset;
    int numref;			// number of references to this storage
    void* store;

    Array1_base();
    void init();
    void deletestore(dodeletefun * dodel);
    void checkindex(int i);
#ifdef HEADER_DEFS
    static Array1_base& Nullbase()
        {
        static Array1_base Nullbase_;
        return Nullbase_;
        }
#endif
    };

class ArrayDo
    {
    void decideMake(int length, int store_is_zero,
                    int& dodel, int& do_new, int& zerostore);
public:
    FunPoint * fp;
    Array1_base ** prep;
    int inuse;

    ArrayDo();
    int Size();
    int StoreSize();
    int Lower();
    int Upper();
    int Numref();

    void ReSize(int limit);	// Index starts at 0
    void ReDimension(int off,int limit = SPECIAL);
    void EnLarge(int off,int limit = SPECIAL);
    void SetSize(int s);	// Reduce size without ReDimensioning
    void ReduceDimension(int off,int limit = SPECIAL);

    std::ostream& outputarray(std::ostream& s);

    void make(int off,int length);
    void deleterep();
    void OwnCopy();
    void OwnNullCopy();
    };

ArrayDo * makedo(Array1_base ** prep, FunPoint * funp);

#ifndef HEADER_DEFS
extern Array1_base Nullbase;
#endif

// *************************************************************

template<class T> class Array1
    {
    Array1_base * rep;

public:
    inline ArrayDo * operator->() const
#ifdef HEADER_DEFS
	{ return makedo((Array1_base **)&rep,T::pfunpoint()); }
#else
	{ return makedo((Array1_base **)&rep,T::pfunpoint); }
#endif

    inline Array1(int off,int limit = SPECIAL)	// Usual usage: Array1(0,lim-1)
	{
	rep = new Array1_base;
	(*this)->ReDimension(off,limit);
	}

    inline Array1(const Array1<T>& other)
	{ rep = other.rep; rep->numref++; }

    inline Array1<T>& operator=(const Array1<T>& other)
	{
	if(rep != other.rep)
	    {
	    (*this)->deleterep();
	    rep = other.rep; 
	    rep->numref++;
	    }
	return *this;
	}

    inline const T& operator()(int i) const
	{
#ifdef BOUNDS
	rep->checkindex(i);
#endif
	const T * store = (T*) rep->store;
	return store[i-rep->offset];
	}

    inline T& operator[](int i)
	{
#ifdef BOUNDS
	rep->checkindex(i);
#endif
	if(rep->numref > 1)
	    (*this)->OwnCopy();
	T * store = (T*) rep->store;
	return store[i-rep->offset];
	}

    inline int Size() const
	{ return rep->size; }

    inline T* Store() const
	{ return (T*) rep->store; }

    inline void operator=(const T& val)
	{
	(*this)->OwnCopy();
#ifdef HEADER_DEFS
	T::pfunpoint()->copyvalue(rep->store,(void*)&val,rep->size);
#else
	T::pfunpoint->copyvalue(rep->store,(void*)&val,rep->size);
#endif
	}

    inline Array1()
#ifdef HEADER_DEFS
	{ rep = &Array1_base::Nullbase(); rep->numref++; }
#else
	{ rep = &Nullbase; rep->numref++; }
#endif

    inline ~Array1()
	{ 
	(*this)->deleterep();
	}

    inline friend std::ostream& operator<<(std::ostream& s, const Array1<T>& V)
	{ return V->outputarray(s); }
    };


#ifdef HEADER_DEFS

#define ARRAY1H_DEFS(T) \
static FunPoint* pfunpoint();

#define ARRAY1CC_DEFS(T) \
inline void T##copyarray( void * vA, const void * vB, int n)\
    {\
    T *A = (T *) vA;\
    const T *B = (const T *) vB;\
    for (int i = 0; i < n; i++)\
	A[i] = B[i];\
    } \
inline void T##copyvalue( void * vA, const void * pvalue, int n)\
    {\
    T *A = (T *) vA;\
    const T& value = *((const T*)pvalue);\
    for (int i = 0; i < n; i++)\
	A[i] = value;\
    } \
inline void * T##donew(int len)\
    { return new T[len]; }\
inline void T##dodelete(void * vA)\
    { delete [] ((T *)vA); }\
inline void T##outputarray( std::ostream& s, void * vA, int n, int offset)\
    {\
    T *A = (T *) vA;\
    for (int i = 0; i < n; i++)\
	s << i + offset << " " << A[i] << iendl;\
    }\
inline FunPoint* T::pfunpoint()\
    {\
    static FunPoint fp_(T##copyarray,T##copyvalue,T##donew,\
                        T##dodelete,T##outputarray);\
    return &fp_;\
    }

#else //HEADER_DEFS undefined

#define ARRAY1H_DEFS(T) \
static FunPoint * pfunpoint;

#define ARRAY1CC_DEFS(T) \
void T##copyarray( void * vA, const void * vB, int n)\
    {\
    T *A = (T *) vA;\
    const T *B = (const T *) vB;\
    for (int i = 0; i < n; i++)\
	A[i] = B[i];\
    } \
void T##copyvalue( void * vA, const void * pvalue, int n)\
    {\
    T *A = (T *) vA;\
    const T& value = *((const T*)pvalue);\
    for (int i = 0; i < n; i++)\
	A[i] = value;\
    } \
void * T##donew(int len)\
    { return new T[len]; }\
void T##dodelete(void * vA)\
    { delete [] ((T *)vA); }\
void T##outputarray( std::ostream& s, void * vA, int n, int offset)\
    {\
    T *A = (T *) vA;\
    for (int i = 0; i < n; i++)\
	s << i + offset << " " << A[i] << iendl;\
    }\
FunPoint T##funpoint(T##copyarray,T##copyvalue,T##donew,\
		    T##dodelete,T##outputarray);\
FunPoint * T::pfunpoint = &T##funpoint;

#endif //HEADER_DEFS

#define ARRAYTYPE(t,T) \
class T\
    {\
public:\
    t a;\
    operator t() const\
	{ return a; }\
    T(t aa) { a = aa; }\
    T() {a = 0; }\
    ARRAY1H_DEFS(T)\
    };\
std::istream & operator >> (std::istream &s, T & x);

#ifdef HEADER_DEFS

#define ARRAYTYPECC(T) \
inline std::istream & operator >> (std::istream &s, T & x)\
    { return s >> x.a; }

#else //HEADER_DEFS undefined

#define ARRAYTYPECC(T) \
std::istream & operator >> (std::istream &s, T & x)\
    { return s >> x.a; }

#endif //HEADER_DEFS

ARRAYTYPE(int,Int)

class IntArray1 : public Array1<Int>
    {
public:
    IntArray1(int n) : Array1<Int>(n) {}
    IntArray1(int n1,int n2) : Array1<Int>(n1,n2) {}
    IntArray1() : Array1<Int>() {}
    IntArray1(const IntArray1& other) : Array1<Int>(other) {}
    void operator=(int n)
	{ Array1<Int>::operator=(Int(n)); }
    IntArray1& operator=(const IntArray1& other)
	{ return (IntArray1&)Array1<Int>::operator=(other); }

    int& operator[](int i)
	{ return Array1<Int>::operator[](i).a; }
    int operator()(int i) const
	{ return Array1<Int>::operator()(i).a; }
    ARRAY1H_DEFS(IntArray1)
    };

#ifdef HEADER_DEFS

ARRAYTYPECC(Int)

ARRAY1CC_DEFS(Int)

ARRAY1CC_DEFS(IntArray1)

#else //ifndef HEADER_DEFS

#ifdef THIS_IS_MAIN

ARRAYTYPECC(Int)

ARRAY1CC_DEFS(Int)
Array1_base Nullbase;

ARRAY1CC_DEFS(IntArray1)
#endif

#endif //HEADER_DEFS

#endif
