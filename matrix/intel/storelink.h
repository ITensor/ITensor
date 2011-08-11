// storelink.h -- Header file for the StoreLink class-- S.R. White 9/94 

#ifndef _storelink_h
#define _storelink_h

#include <iostream>
using namespace std;

typedef double Real;

// The Matrix/Vector Ref classes have a StoreLink, which prevents the 
// storage on which they are based from being deleted prematurely. 
// StoreLink utilizes reference counting. The ref classes never 
// allocate storage. The actual storage classes utilize makestorage, 
// etc. for allocation.

class StoreReport;

// Actual StoreLink structure
struct storerep
    {
    int storage;			// Size of storage
    int numref;				// Number of references 
    storerep() : storage(0), numref(1) {}
    };

class StoreLink
    {
public:
// The following are used both by ref and storage classes.
// The only ways to make a new StoreLink from an old one increment numref.
// Reserve "=" for copies of storage pointed to, defined in later classes
// using StoreLink. Do ref copy (rebinding) using <<
    inline Real * Store() const;	// Don't ever delete this on your own.
    inline StoreLink(const StoreLink &);		// New link.
    inline StoreLink & operator<<(const StoreLink &);	// Remove old link,
							// copy new one.
// Commands for new storage, used by storage classes only.
    inline StoreLink(int);		// Negative int treated as 0.
    inline void makestorage(int);	// Resize storage to int.
    inline void increasestorage(int);	// Increase size to int, no reduce.
    	
    int defragment(int newsize);	// Tries to move storage to lower 
    					// place in heap. Returns 1 if
					// successful, 0 otherwise.
    inline int memory() const;		// Memory allocated by this object.
    	
// Miscellaneous functions
    inline int NumRef() const;
    inline int Storage() const;
    inline StoreLink();
    inline ~StoreLink();
    inline static int NumObjects();
    inline static int TotalStorage();
    friend class StoreReport;
private:
    storerep *p;			// Only data member

    static int storageinuse;
    static int numberofobjects;		// Number of new's - no. of deletes
    static storerep nullrep; 
    static storerep * pnullrep; 

    enum { offset = (sizeof(storerep)-1) / sizeof(Real) + 1 };
    inline void donew(int s);
    inline void dodelete();
// " =" is private, not allowed.  Put in to replace default shallow copy.
    inline StoreLink & operator = (const StoreLink &); 
    };
class StoreReport
    {
    int i;
public:
    friend class StoreLink;
    StoreReport() 
    	{ i = StoreLink::TotalStorage(); }
    ~StoreReport()
    	{ 
	/*
	cout << "StoreReport: When created, total storage = " << i 
	     << ", now total storage = " << 
	     StoreLink::TotalStorage() << endl; 
	*/
	}
    };

// Inline functions for StoreLink

inline void StoreLink::donew(int s)
    {
    if (s > 0)
	{
	p = (storerep *) new Real[s + offset];
	p->numref = 1; p->storage = s; storageinuse += s;
	numberofobjects++;
	// cout << "Making storage address " << (long)(p) << endl;
	}
    else  
	{ p = pnullrep; p->numref++; }
    }

inline void StoreLink::dodelete()
    { 
    if(--p->numref == 0) 
	{
	// cout << "Deleting storage address " << (long)(p) << endl;
	storageinuse -= p->storage; numberofobjects--;
	delete [] ((Real *) p);
//	if(storageinuse <= 0)
//	    cout << "Storage in use is now " << storageinuse << endl;
	}
    }

inline StoreLink::StoreLink() : p(pnullrep)
    { p->numref++; }

inline Real * StoreLink::Store() const
    { return ((Real *)p)+offset; }

inline int StoreLink::NumRef() const { return p->numref; }

inline int StoreLink::Storage() const { return p->storage; }

inline StoreLink::~StoreLink() { dodelete(); }

inline StoreLink::StoreLink(const StoreLink & S) : p(S.p)
    { p->numref++; }

inline StoreLink & StoreLink::operator<<(const StoreLink & S)		
    { 			
    if(this != &S) { dodelete(); p = S.p; p->numref++; }
    return *this; 
    }

inline StoreLink::StoreLink(int s) 
    { donew(s); }

inline void StoreLink::makestorage(int s)	// Negative s treated as 0
    {
    if(p->storage != s) { dodelete(); donew(s); }
    }

inline void StoreLink::increasestorage(int s)
    {
    if(p->storage < s) { dodelete(); donew(s); }
    }

inline int StoreLink::memory() const
    { return sizeof(Real)*(Storage()+offset); }

inline int StoreLink::TotalStorage() { return storageinuse; }

inline int StoreLink::NumObjects() { return numberofobjects; }

inline StoreLink & StoreLink::operator = (const StoreLink & other)
    { return *this << other; } 		// private member function!

// Static member definitions
#ifdef THIS_IS_MAIN

int StoreLink::storageinuse = 0;
int StoreLink::numberofobjects = 0;
// The following has one initial reference, so it never goes away.
storerep StoreLink::nullrep;
storerep * StoreLink::pnullrep = &StoreLink::nullrep;
StoreReport doreport;

#endif

#endif

