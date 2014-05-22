//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_PERMUTATION_H
#define __ITENSOR_PERMUTATION_H
#include "global.h"

namespace itensor {

//
// Tell where each index will go, 
// if(p.dest(2) == 1) then 2 -> 1, etc.
//
class Permutation
    {
    public:
    typedef array<int,NMAX+1> 
    int9;
    
    Permutation();

    Permutation(int i1, int i2 = 2, int i3 = 3, int i4 = 4, 
                int i5 = 5, int i6 = 6,int i7 = 7, int i8 = 8);

    bool 
    isTrivial() const { return trivial; }

    void 
    fromTo(int j, int k);

    int 
    dest(int j) const { return GET(ind_,j); }

    bool 
    check(int d);

    const int9& 
    ind() const { return ind_; }

    private:

    ///////////
    int9 ind_;

    bool trivial;
    //////////

    void 
    set8(int9 *n, int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8);

    };

Permutation inline
inverse(const Permutation& P)
    {
    Permutation inv;
    for(int n = 1; n <= NMAX; ++n) 
        inv.fromTo(P.dest(n),n);
    return inv;
    }

inline Permutation::
Permutation() 
    : trivial(true) 
    { set8(&ind_,1,2,3,4,5,6,7,8); }

inline Permutation::
Permutation(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8)
    : trivial(i1==1 && i2==2 && i3==3 && i4==4 && i5==5 && i6==6 && i7==7 && i8==8)
	{ set8(&ind_,i1,i2,i3,i4,i5,i6,i7,i8); }

void inline Permutation::
fromTo(int j, int k) 
    { 
    if(j!=k) { trivial = false; } 
    GET(ind_,j) = k; 
    }

bool inline Permutation::
check(int d)
	{
    for(int i = 1; i <= d; i++)
        {
        if(ind_[i] > d || ind_[i] < 1) 
            {
            std::cerr << "\nbad Permutation level 1\n\n";
            return false;
            }
        }

    for(int i = 1; i <= d; i++)
    for(int j = 1; j <= d; j++)
        {
        if(i == j) continue;
        if(ind_[i] == ind_[j]) 
            {
            std::cerr << "\nbad Permutation level 2\n\n";
            return false;
            }
        }

    return true;
	}

void inline Permutation::
set8(int9 *n, int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8)
    {
    (*n)[1] = i1; (*n)[2] = i2; (*n)[3] = i3; (*n)[4] = i4;
    (*n)[5] = i5; (*n)[6] = i6; (*n)[7] = i7; (*n)[8] = i8;
    }

inline std::ostream& 
operator<<(std::ostream& s, const Permutation& p)
    {
    for(int i = 1; i <= NMAX; ++i) 
        s << "(" << i << "," << p.dest(i) << ")";
    return s;
    }

}; //namespace itensor

#endif
