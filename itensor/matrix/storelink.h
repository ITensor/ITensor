// storelink.h -- Header file for the StoreLink class-- S.R. White 9/94 

#ifndef _storelink_h
#define _storelink_h

#include <vector>
#include "itensor/util/types.h"


namespace itensor {

// The Matrix/Vector Ref classes have a StoreLink, which prevents the 
// storage on which they are based from being deleted prematurely. 
// StoreLink utilizes reference counting through a shared_ptr. The ref classes never 
// allocate storage. The actual storage classes utilize makestorage, 
// etc. for allocation, implemented with a vector<Real>

typedef double Real;

#ifdef USE_CPP11
#define STATIC_VAR(type,name) static type& name() { static type t{}; return t; }
#else
#define STATIC_VAR(type,name) static type& name() { static type t(0); return t; }
#endif

class StoreLink
    {
    shared_ptr<std::vector<Real> > p; // Only data member
    public:

    typedef long int lint;

    // The following are used both by ref and storage classes.
    // Reserve "=" for copies of storage pointed to, defined in later classes
    // using StoreLink. Do ref copy (rebinding) using <<

#ifdef USE_CPP11
    StoreLink(const StoreLink & S) = default;
#endif

    // Commands for new storage, used by storage classes only.
    StoreLink(lint s = 0) // Negative lint treated as 0.
        { 
        donew(s); 
        }


    ~StoreLink() 
        { 
        dodelete(); 
        }

    Real* 
    Store() const	// Don't ever delete this on your own.
        {
#ifdef USE_CPP11
        if(!p || p->empty()) return nullptr;
#else
        if(!p || p->empty()) return NULL;
#endif
        return &(p->front());
        }

    StoreLink& 
    operator<<(const StoreLink &S) // Remove old link, copy new one.
        { 
        if(this != &S) { dodelete(); p = S.p; }
        return *this; 
        }

    void 
    makestorage(lint s)	// Resize storage to lint. Negative s treated as 0.
        { 
        if(Storage() != s) 
            { 
            dodelete(); 
            donew(s); 
            } 
        }

    void 
    increasestorage(lint s) // Increase size to lint, no reduce.
        { 
        if(Storage() <  s) 
            { 
            dodelete(); 
            donew(s); 
            } 
        }

    lint 
    defragment(lint newsize); // Tries to move storage to lower 
                              // place in heap. Returns 1 if
                              // successful, 0 otherwise.

    // Information functions
    lint 
    NumRef() const 
        { return p ? p.use_count() : 0; }

    lint 
    Storage() const 
        { return p ? p->size() : 0; }

    lint memory() const			// Memory allocated by this object.
        { 
        return sizeof(Real)*(Storage()) + sizeof(std::vector<Real>) + sizeof(p); 
        }

    lint static
    TotalStorage() { return storageinuse(); }

    lint static
    NumObjects() { return numberofobjects(); }

    private:

    StoreLink& operator=(const StoreLink &);

    STATIC_VAR(lint,storageinuse)
    STATIC_VAR(lint,numberofobjects)

    void donew(lint s)
        {
        if(s <= 0) return;
        p = make_shared<std::vector<Real> >(s);
        storageinuse() += s;
        ++numberofobjects();
        }

    void 
    dodelete()
        { 
        if(p && p.unique())
            {
            storageinuse() -= p->size(); 
            --numberofobjects();
            }
        p.reset();
        }
    };

} //namespace itensor

#undef STATIC_VAR

#endif

