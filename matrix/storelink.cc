// storelink.cc -- Code for StoreLink class

#include <stdlib.h>
#include <memory.h>
#include "storelink.h"
#include "minmax.h"

int 
StoreLink::defragment(int newsize)	// move object to lower place in heap 
    {
    if(NumRef() != 1) return 0;
    newsize = min(newsize,Storage());
    if(newsize < 1) return 0;
    int rsize = sizeof(Real)*newsize;
    const int stacksize = 10000;
    if(newsize <= stacksize)
    	{
	Real tempstore[stacksize];
	memcpy((void *) tempstore, (void*) Store(), rsize);
	dodelete();
	donew(newsize);
	memcpy((void *) Store(), (void*) tempstore, rsize);
	return 1;
	}
    else
    	{
	StoreLink S(newsize);
	if(S.Store() != 0 && S.Store() < Store())
	    {
	    memcpy((void*) S.Store(), (void*) Store(), rsize);
	    *this << S;
	    return 1;
	    }
	return 0;
	}
    }
