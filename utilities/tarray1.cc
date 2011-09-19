#include "tarray1.h"

// Private Null defs
Array1_base:: Array1_base(const Array1_base&) {}
void Array1_base::operator=(const Array1_base&) {} 

Array1_base::Array1_base() 
    { init(); }

void Array1_base::init() 
    { size = storesize = 0; offset = 1; store = 0; numref = 1; }

void Array1_base::deletestore(dodeletefun * dodel)
    {
    if(store != 0)
	dodel(store);
    size = storesize = 0;
    }

void Array1_base::checkindex(int i) 
    {
    if (i < offset || i > offset + size - 1)
	{
	std::cout << "i, offset, size are " << i SP offset SP size << std::endl;
	error("Array1rep<T>: Index out of bounds");
	}
    }

ArrayDo::ArrayDo() : fp(0), prep(0), inuse(0) {}

void ArrayDo::deleterep()
    {
    if((*prep)->numref-- == 1)
	{
	(*prep)->deletestore(fp->dodelete);
	delete *prep;
	}
    inuse--;
    }

void ArrayDo::make(int off,int length)
    {
    int dodel, zerostore, do_new;
    decideMake(length,(*prep)->store == 0,dodel,do_new,zerostore);
    if(dodel) fp->dodelete((*prep)->store);
    if(do_new) (*prep)->store = fp->donew(length);
    if(zerostore) (*prep)->store = 0;
    (*prep)->offset = off;
    inuse--;
    }

void ArrayDo::decideMake(int length, int store_is_zero,
                int& dodel, int& do_new, int& zerostore)
    {
    if (length < 0)
        error("ArrayDo::make: length < 0");
    dodel = zerostore = do_new = 0;
    if(store_is_zero)
        do_new = (length > 0);
    else if(length == 0)
        {
        dodel = 1;
        zerostore = 1;
	(*prep)->storesize = 0;
        }
    else if(length != (*prep)->size && length != (*prep)->storesize)
        {
        dodel = 1;
        do_new = 1;
        }
    (*prep)->size = length;
    if(do_new) (*prep)->storesize = length;
    }

int ArrayDo::Size() 
    { inuse--; return (*prep)->size; }

int ArrayDo::StoreSize() 
    { inuse--; return (*prep)->storesize; }

int ArrayDo::Lower() 
    { inuse--; return (*prep)->offset; }

int ArrayDo::Upper() 
    { inuse--; return (*prep)->offset+(*prep)->size-1; }

int ArrayDo::Numref() 
    { inuse--; return (*prep)->numref; }

void ArrayDo::ReSize(int limit)		// Index starts at 0
    { ReDimension(0,limit-1); }

void ArrayDo::OwnCopy()
    {			      // make sure it has it's own copy
    Array1_base * rep = *prep;
    if(rep->numref > 1)
	{
	rep->numref--;
	*prep = new Array1_base;
	inuse++;
	make(rep->offset,rep->size);
	fp->copyarray((*prep)->store,rep->store,rep->size);
	}
    inuse--;
    }

void ArrayDo::OwnNullCopy()
    {
    if((*prep)->numref > 1)
	{
	(*prep)->numref--;
	*prep = new Array1_base;
	}
    inuse--;
    }

void ArrayDo::ReDimension(int off,int limit)
    {
    if(limit == SPECIAL)
	{ limit = off; off = 1; }
    inuse++;
    OwnNullCopy();
    inuse++;
    make(off,limit-off+1);
    inuse--;
    }

// Reduce size without ReDimensioning
void ArrayDo::SetSize(int s)
    {
    if (s <= (*prep)->size && s >= 0) (*prep)->size = s;
    else error("ArrayDo::SetSize: Bad new size");
    inuse--;
    }

void ArrayDo::EnLarge(int off,int limit) // Assume offset stays same
    {
    int length;
    if(limit == SPECIAL)
	{
	length = off;
	off = (*prep)->offset;
	limit = length + off - 1;
	}
    else
	{
	length = limit-off+1;
	}

    if(length < 0)
	error("length < 0 in Array1::EnLarge");
    if(length == 0)
	{
	inuse++;
	deleterep();
	(*prep) = &Nullbase;
	(*prep)->numref++;
	inuse--;
	return;
	}
    if(length <= (*prep)->storesize && (*prep)->numref == 1)	// no new storage needed
	{
	(*prep)->size = length;
	(*prep)->offset = off;
	}
    else
	{
	Array1_base * rep = new Array1_base;
	rep->size = length;
	rep->offset = off;
	rep->storesize = length;
	rep->store = fp->donew(length);
	int ns = (*prep)->size; 
	int lim = (length < ns ? length : ns);
	fp->copyarray(rep->store,(*prep)->store,lim);
	inuse++;
	deleterep();
	*prep = rep;
	}
    inuse--;
    }
	
void ArrayDo::ReduceDimension(int off,int limit)
    {
    if(limit == SPECIAL)
	{ limit = off; off = 1; }
    inuse++;
    OwnNullCopy();
    int length = limit-off+1;
    if(length <= (*prep)->storesize)
	{
	(*prep)->size = length;
	(*prep)->offset = off;
	}
    else
	{
	inuse++;
	make(off,limit-off+1);
	}
    inuse--;
    }

std::ostream& ArrayDo::outputarray(std::ostream& s)
    {
    //long f = s.flags();
    fp->outputarray(s,(*prep)->store,(*prep)->size,(*prep)->offset);
    s << iendl;
    //s.flags(f);
    inuse--;
    return s;
    }

class ArrayDoNext
    { public: ArrayDo ado; ArrayDoNext *next; ArrayDoNext(){} };
class ADNArray
    { public: ArrayDoNext adnarray[300]; ADNArray(); };

ADNArray::ADNArray()
    {
    for(int i=0; i < 299; i++)
	adnarray[i].next = adnarray + i + 1;
    adnarray[299].next = adnarray;
    }

ArrayDo * makedo(Array1_base ** prep, FunPoint * funp)
    {
    static ADNArray ADN;
    static ArrayDoNext * curadn = ADN.adnarray;

#ifdef BOUNDS
    int cnt = 0;
#endif
    while(curadn->ado.inuse > 0)
	{
	curadn = curadn->next;
#ifdef BOUNDS
	if(cnt++ > 300)
	    error("bad index cnt in makedo");
#endif
	}
#ifdef BOUNDS
    if(curadn->ado.inuse < 0)
	error("bad inuse in makedo");
#endif
    curadn->ado.fp = funp;
    curadn->ado.prep = prep;
    curadn->ado.inuse = 1;
    return (ArrayDo *)&curadn->ado;
    }

