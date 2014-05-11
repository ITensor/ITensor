//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __COMBINER_H
#define __COMBINER_H

#include "itensor.h"

namespace itensor {

/*
Combine several indices into one, use * to convert tensors efficiently
   \
    \
  ---C====
    /
   /

*/

class Combiner;

Combiner
primed(Combiner C, int inc = 1);

class Combiner
    {
    public:

    typedef boost::array<Index,NMAX+1>::const_iterator 
    left_it;

    //Accessor Methods ----------------------------------------------

    inline const Index& 
    right() const 
        { init(); return right_; }

    int 
    numLeft() const { return rl_; }

    const Index& 
    left(int j) const { return GET(left_,j); }

    const std::pair<left_it,left_it> 
    left() const 
        { 
        return std::make_pair(left_.begin()+1,left_.begin()+rl_+1); 
        }

    //Constructors --------------------------------------------------

    Combiner() : rl_(0), initted(false) {}

    Combiner(const Index& l1, const Index& l2 = Index::Null(), 
             const Index& l3 = Index::Null(), const Index& l4 = Index::Null(), 
             const Index& l5 = Index::Null(), const Index& l6 = Index::Null(), 
             const Index& l7 = Index::Null(), const Index& l8 = Index::Null() );

    Combiner(const Combiner& other) { operator=(other); }

    //Operators -----------------------------------------------------

    ITensor 
    operator*(const ITensor& t) const 
        { ITensor res; product(t,res); return res; }

    Combiner&
    operator=(const Combiner& other);

    //Index Methods -------------------------------------------------

    inline bool 
    isInit() const { return initted; }

    void 
    reset();

    void 
    addleft(const Index& l);// Include another left index

    //Call addleft on a vector of Indices ls[0]..ls[ls.size()-1]
    void
    addleft(const std::vector<Index>& ls);

    //Initialize after all lefts are added and before being used
    void 
    init(std::string rname = "cmb", 
         IndexType type = Link, 
         Arrow dir = Out,
         int primelevel = 0) const;

    void
    prime(int inc = 1);

    void
    prime(IndexType type, int inc = 1);

    friend Combiner
    primed(Combiner C, int inc);


    //Other Methods -------------------------------------------------

    Real
    uniqueReal() const;

    operator ITensor() const;

    friend inline std::ostream& 
    operator<<(std::ostream & s, const Combiner & c);

    void 
    conj() { init(); }

    void 
    product(const ITensor& t, ITensor& res) const;

    //For interface compatibility with IQCombiner
    void 
    doCondense(bool) { } 

    private:

    /////////

    boost::array<Index,NMAX+1> left_; // max dim is 8
    mutable Index right_;
    int rl_; //Number of m>1 'left' indices (indices to be combined into one)
    mutable bool initted;

    ///////

    }; //class Combiner

Combiner inline
conj(Combiner res) { res.conj(); return res; }

ITensor inline
operator*(const ITensor& t, const Combiner& c) { return c*t; }

inline
Combiner& Combiner::
operator=(const Combiner& other)
    {
    other.init();
    left_ = other.left_;
    right_ = other.right_;
    rl_ = other.rl_;
    initted = other.initted;
    return *this;
    }

inline
Combiner::
Combiner(const Index& l1, const Index& l2,
         const Index& l3, const Index& l4, 
         const Index& l5, const Index& l6, 
         const Index& l7, const Index& l8)
    : 
    rl_(0), 
    initted(false)
	{
    boost::array<const Index*,NMAX+1> ll 
    = {{ &Index::Null(), &l1, &l2, &l3, &l4, &l5, &l6, &l7, &l8 }};

    do { ++rl_; left_[rl_] = *ll[rl_]; } 
    while(rl_ < NMAX && *ll[rl_+1] != Index::Null());

    assert(rl_ == NMAX || left_[rl_+1] == Index::Null());
    assert(left_[rl_] != Index::Null());
	}

inline
void Combiner::
reset()
    {
    rl_ = 0;
    initted = false;
    }

void inline Combiner::
addleft(const Index& l)// Include another left index
    { 
    initted = false;
    if(rl_ == NMAX) 
        Error("Combiner: already reached max number of left indices.");
    ++rl_;
    left_[rl_] = l; 
    }

void inline Combiner::
addleft(const std::vector<Index>& ls)
    { 
    initted = false;
    if(rl_+int(ls.size()) > NMAX) 
        Error("Combiner: too many left indices.");
    for(size_t j = 0; j < ls.size(); ++j)
        left_[++rl_] = ls[j]; 
    }

inline
void Combiner::
init(std::string rname, IndexType type, Arrow dir, int primelevel) const
    {
    if(initted) return;
    int m = 1; 
    for(int i = 1; i <= rl_; ++i) 
        { m *= left_[i].m(); }
    right_ = Index(rname,m,type,primelevel); 
    initted = true;
    }

void inline Combiner::
prime(int inc)
    {
    for(int j = 1; j <= rl_; ++j) 
        {
        left_[j].prime(inc);
        }
    if(initted)
        {
        right_.prime(inc);
        }
    }

void inline Combiner::
prime(IndexType type, int inc)
    {
    for(int j = 1; j <= rl_; ++j) 
        {
        left_[j].prime(type,inc);
        }
    if(initted)
        {
        right_.prime(type,inc);
        }
    }

Combiner inline
primed(Combiner C, int inc)
    {
    C.prime(All,inc);
    return C;
    }

inline
Combiner::
operator ITensor() const
    {
    /*
    if(right_.m() > 16) 
    { 
        std::cerr << "\n\n" 
        << "WARNING: too large of an m in Combiner to ITensor!\n\n"; 
    }
    */

    //Use a kronecker delta tensor to convert this Combiner into an Tensor
    ITensor res = operator*(ITensor(right_,primed(right_,5),1));
    res.prime(primed(right_,5),-5);
    return res;
    }

inline
void Combiner::
product(const ITensor& t, ITensor& res) const
    {
    init();

    if(hasindex(t,right_))
        {
        IndexSet<Index> nind;
        Foreach(const Index& I, t.indices())
            {
            if(I == right_)
                {
                for(int i = 1; i <= rl_; ++i)
                    nind.addindex(left_[i]);
                }
            else
                {
                nind.addindex(I);
                }
            }
        res = ITensor(nind,t);
        return;
        }

    t.groupIndices(left_,rl_,right_,res);
    }

Real inline Combiner::
uniqueReal() const
    {
    Real ur = 0;
    for(int j = 1; j <= rl_; ++j)
        ur += left_[j].uniqueReal();
    return ur;
    }

//
// Combiner helper method
//

bool inline
hasindex(const Combiner& C, const Index& I)
    {
    for(int j = 1; j <= C.numLeft(); ++j) 
        if(C.left(j) == I) return true;
    return false;
    }



inline 
std::ostream& 
operator<<(std::ostream & s, const Combiner & c)
    {
    if(c.isInit())
        s << "\nRight index: " << c.right() << "\n";
    else
        s << "\nRight index not initialized" << "\n";
    s << "Left indices:\n";
    Foreach(const Index& l, c.left()) s << " " << l << "\n";
    return s;
    }

}; //namespace itensor


#endif
