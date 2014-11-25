//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __COMBINER_H
#define __COMBINER_H

#include "itensor.h"
#include "iterpair.h"

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
prime(Combiner C, int inc = 1);

class Combiner
    {
    public:

    //Constructors --------------------------------------------------

    Combiner() : initted(false) {}

    template<typename... Inds>
    Combiner(const Index& l1, 
             const Inds&... inds);

    Combiner(const Combiner& other) { operator=(other); }

    //Accessor Methods ----------------------------------------------

    const Index& 
    right() const { init(); return right_; }

    int 
    numLeft() const { return left_.size(); }

    //1-indexed access to left indices
    const Index& 
    left(int j) const { return GET(left_,j-1); }

    const std::vector<Index>&
    left() const { return left_; }

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


    //Other Methods -------------------------------------------------

    ITensor
    toITensor() const;

    void 
    dag() { init(); }

    void 
    product(const ITensor& t, ITensor& res) const;

    //For interface compatibility with IQCombiner
    void 
    doCondense(bool) { } 

    private:

    /////////

    mutable bool initted;
    std::vector<Index> left_;
    mutable Index right_;

    ///////

    }; //class Combiner

Combiner inline
dag(Combiner res) { res.dag(); return res; }

ITensor inline
operator*(const ITensor& t, const Combiner& c) { return c*t; }

inline
Combiner& Combiner::
operator=(const Combiner& other)
    {
    other.init();
    left_ = other.left_;
    right_ = other.right_;
    initted = other.initted;
    return *this;
    }

template<typename... Inds>
Combiner::
Combiner(const Index& l1, 
         const Inds&... inds);
    : 
    initted(false),
    left_{l1,inds...}
	{
	}

inline
void Combiner::
reset()
    {
    left_.clear();
    initted = false;
    }

void inline Combiner::
addleft(const Index& l)// Include another left index
    { 
    initted = false;
    left_.push_back(l); 
    }

void inline Combiner::
addleft(const std::vector<Index>& ls)
    { 
    initted = false;
    left_.insert(left_.end(),ls.begin(),ls.end());
    }

inline
void Combiner::
init(std::string rname, IndexType type, Arrow dir, int primelevel) const
    {
    if(initted) return;
    int m = 1; 
    for(const auto& ll : left_) m *= ll.m();
    right_ = Index(rname,m,type,primelevel); 
    initted = true;
    }

void inline Combiner::
prime(int inc)
    {
    for(auto& ll : left_) ll.prime(inc);
    if(initted) right_.prime(inc);
    }

void inline Combiner::
prime(IndexType type, int inc)
    {
    for(auto& ll : left_) ll.prime(type,inc);
    if(initted) right_.prime(type,inc);
    }

Combiner inline
prime(Combiner C, int inc)
    {
    C.prime(All,inc);
    return C;
    }

ITensor inline
Combiner::
toITensor() const
    {
    /*
    if(right_.m() > 16) 
    { 
        std::cerr << "\n\n" 
        << "WARNING: too large of an m in Combiner to ITensor!\n\n"; 
    }
    */

    //Use a kronecker delta tensor to convert this Combiner into an Tensor
    Index rP = right_;
    rP.prime(5);
    ITensor res = operator*(ITensor(right_,rP,1));
    res.prime(rP,-5);
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
        for(const Index& I : t.indices())
            {
            if(I == right_)
                {
                for(const auto& ll : left_)
                    nind.addindex(ll);
                }
            else
                {
                nind.addindex(I);
                }
            }
        res = ITensor(nind,t);
        return;
        }

    Error("Need to change groupIndices to assume 0-indexed container");
    t.groupIndices(left_,left_.size(),right_,res);
    }

//
// Combiner helper method
//

int inline
hasindex(const Combiner& C, const Index& I)
    {
    for(int j = 1; j <= C.numLeft(); ++j) 
        if(C.left(j) == I) return j;
    return 0;
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
    for(const Index& l : c.left()) s << " " << l << "\n";
    return s;
    }

//Deprecated, for backwards compatibility only:
Combiner inline
primed(Combiner C, int inc = 1)
    {
    return prime(C,inc);
    }

}; //namespace itensor


#endif
