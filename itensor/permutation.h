//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __PERMUTATION_H
#define __PERMUTATION_H
#include "global.h"

//
// Tell where each index will go, 
// if(p.dest(2) == 1) then 2 -> 1, etc.
//
class Permutation
{
public:
    typedef boost::array<int,NMAX+1> int9;
    
    inline const int9& 
    ind() const { return ind_; }

    inline bool 
    is_trivial() const { return trivial; }

    Permutation();

    Permutation(int i1, int i2 = 2, int i3 = 3, int i4 = 4, 
                int i5 = 5, int i6 = 6,int i7 = 7, int i8 = 8);

    void 
    from_to(int j, int k);

    inline int 
    dest(int j) const { return GET(ind_,j); }

    bool 
    check(int d);

    friend inline std::ostream& 
    operator<<(std::ostream& s, const Permutation& p);

private:
    void set8(int9 *n, int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8);

    int9 ind_;

    bool trivial;
};

inline Permutation inverse(const Permutation& P)
    {
    Permutation inv;
    for(int n = 1; n <= NMAX; ++n) 
        inv.from_to(P.dest(n),n);
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

inline void Permutation::
from_to(int j, int k) 
    { 
    if(j!=k) { trivial = false; } 
    GET(ind_,j) = k; 
    }

inline bool Permutation::
check(int d)
	{
        for(int i = 1; i <= d; i++)
        if(ind_[i] > d || ind_[i] < 1) 
        {
        std::cerr << "\nbad Permutation level 1\n\n";
        return false;
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

inline void Permutation::
set8(int9 *n, int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8)
    {
    (*n)[1] = i1; (*n)[2] = i2; (*n)[3] = i3; (*n)[4] = i4;
    (*n)[5] = i5; (*n)[6] = i6; (*n)[7] = i7; (*n)[8] = i8;
    }

inline std::ostream& 
operator<<(std::ostream& s, const Permutation& p)
    {
    for(int i = 1; i <= NMAX; i++) s << "(" << i << "," << p.ind_[i] << ")";
    return s;
    }

#endif
