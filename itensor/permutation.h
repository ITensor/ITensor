#ifndef __PERMUTATION_H
#define __PERMUTATION_H
#include "types.h"
#include <error.h> //utilities
#include "boost/array.hpp"

class Permutation // Tell where each index will go, p(2,1,3) says 1 -> 2, 2 -> 1, 3 -> 3
{
public:
    typedef boost::array<int,NMAX+1> int9;
private:
    void set8(int9 *n, int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8)
    {
        (*n)[1] = i1; (*n)[2] = i2; (*n)[3] = i3; (*n)[4] = i4;
        (*n)[5] = i5; (*n)[6] = i6; (*n)[7] = i7; (*n)[8] = i8;
    }
    int9 ind_;
    bool trivial;
public:
    const int9& ind() const { return ind_; }
    bool is_trivial() const { return trivial; }

    Permutation() : trivial(true) { set8(&ind_,1,2,3,4,5,6,7,8); }
    Permutation(int i1, int i2 = 2, int i3 = 3, int i4 = 4, int i5 = 5, int i6 = 6,
	    int i7 = 7,  int i8 = 8)
    : trivial(i1==1 && i2==2 && i3==3 && i4==4 && i5==5 && i6==6 && i7==7 && i8==8)
	{ set8(&ind_,i1,i2,i3,i4,i5,i6,i7,i8); }

    void from_to(int j, int k) { if(j!=k) { trivial = false; } GET(ind_,j) = k; }
    inline int dest(int j) const { return GET(ind_,j); }

    void check(int d)
	{
        for(int i = 1; i <= d; i++)
        if(ind_[i] > d || ind_[i] < 1) Error("bad Permutation level 1");

        for(int i = 1; i <= d; i++)
        for(int j = 1; j <= d; j++)
        if(ind_[i] == ind_[j] && i != j) Error("bad Permutation level 2");
	}

    friend inline ostream& operator<<(ostream& s, const Permutation& p)
    {
        for(int i = 1; i <= NMAX; i++) s << "(" << i << "," << p.ind_[i] << ")";
        return s;
    }
};

inline Permutation inverse(const Permutation& P)
{
    Permutation inv;
    for(int n = 1; n <= NMAX; ++n) inv.from_to(P.dest(n),n);
    return inv;
}

#endif
