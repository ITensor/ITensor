//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_LOCAL_OP
#define __ITENSOR_LOCAL_OP
#include "itensor/iqtensor.h"

namespace itensor {

//
// The LocalOp class represents
// an MPO or other operator that
// has been projected into the
// reduced Hilbert space of 
// two sites of an MPS.
//
//   .-              -.
//   |    |      |    |
//   L - Op1 -- Op2 - R
//   |    |      |    |
//   '-              -'
//
// (Note that L, Op1, Op2 and R
//  are not required to have this
//  precise structure. L and R
//  can even be null in which case
//  they will not be used.)
//


template <class Tensor>
class LocalOp
    {
    public:

    using IndexT = typename Tensor::index_type;

    //
    // Constructors
    //

    LocalOp(const Args& args = Global::args());

    LocalOp(const Tensor& Op1, const Tensor& Op2,
            const Args& args = Global::args());

    LocalOp(const Tensor& Op1, const Tensor& Op2, 
            const Tensor& L, const Tensor& R,
            const Args& args = Global::args());

    //
    // Sparse Matrix Methods
    //

    void
    product(const Tensor& phi, Tensor& phip) const;

    Real
    expect(const Tensor& phi) const;

    Tensor
    deltaRho(const Tensor& rho, 
             const Tensor& combine, Direction dir) const;

    Tensor
    diag() const;

    long
    size() const;

    //
    // Accessor Methods
    //

    void
    update(const Tensor& Op1, const Tensor& Op2);

    void
    update(const Tensor& Op1, const Tensor& Op2, 
           const Tensor& L, const Tensor& R);

    const Tensor&
    Op1() const 
        { 
        if(!(*this)) Error("LocalOp is default constructed");
        return *Op1_;
        }

    const Tensor&
    Op2() const 
        { 
        if(!(*this)) Error("LocalOp is default constructed");
        return *Op2_;
        }

    const Tensor&
    L() const 
        { 
        if(!(*this)) Error("LocalOp is default constructed");
        return *L_;
        }

    const Tensor&
    R() const 
        { 
        if(!(*this)) Error("LocalOp is default constructed");
        return *R_;
        }

    const Tensor&
    bondTensor() const 
        { 
        return (*Op1_) * (*Op2_);
        }

    explicit operator bool() const { return bool(Op1_); }

    bool
    LIsNull() const;

    bool
    RIsNull() const;

    void
    operator=(const LocalOp& other)
        {
        Op1_ = other.Op1_;
        Op2_ = other.Op2_;
        L_ = other.L_;
        R_ = other.R_;
        }

    private:

    /////////////////
    const Tensor *Op1_, *Op2_; 
    const Tensor *L_, *R_; 
    mutable long size_;
    /////////////////

    void
    makeBond() const;

    };

template <class Tensor>
inline LocalOp<Tensor>::
LocalOp(const Args& args)
    :
    Op1_(nullptr),
    Op2_(nullptr),
    L_(nullptr),
    R_(nullptr),
    size_(-1)
    { 
    }

template <class Tensor>
inline LocalOp<Tensor>::
LocalOp(const Tensor& Op1, const Tensor& Op2,
        const Args& args)
    : 
    Op1_(nullptr),
    Op2_(nullptr),
    L_(nullptr),
    R_(nullptr),
    size_(-1)
    {
    update(Op1,Op2);
    }

template <class Tensor>
inline LocalOp<Tensor>::
LocalOp(const Tensor& Op1, const Tensor& Op2, 
        const Tensor& L, const Tensor& R,
        const Args& args)
    : 
    Op1_(nullptr),
    Op2_(nullptr),
    L_(nullptr),
    R_(nullptr),
    size_(-1)
    {
    update(Op1,Op2,L,R);
    }

template <class Tensor>
void inline LocalOp<Tensor>::
update(const Tensor& Op1, const Tensor& Op2)
    {
    Op1_ = &Op1;
    Op2_ = &Op2;
    L_ = nullptr;
    R_ = nullptr;
    size_ = -1;
    }

template <class Tensor>
void inline LocalOp<Tensor>::
update(const Tensor& Op1, const Tensor& Op2, 
       const Tensor& L, const Tensor& R)
    {
    update(Op1,Op2);
    L_ = &L;
    R_ = &R;
    }

template <class Tensor>
bool inline LocalOp<Tensor>::
LIsNull() const
    {
    if(L_ == nullptr) return true;
    return !bool(*L_);
    }

template <class Tensor>
bool inline LocalOp<Tensor>::
RIsNull() const
    {
    if(R_ == nullptr) return true;
    return !bool(R_);
    }

template <class Tensor>
void inline LocalOp<Tensor>::
product(const Tensor& phi, Tensor& phip) const
    {
    if(!(*this)) Error("LocalOp is null");

    const Tensor& Op1 = *Op1_;
    const Tensor& Op2 = *Op2_;

    if(LIsNull())
        {
        phip = phi;

        if(!RIsNull()) 
            phip *= R(); //m^3 k d

        phip *= Op2; //m^2 k^2
        phip *= Op1; //m^2 k^2
        }
    else
        {
        phip = phi * L(); //m^3 k d

        phip *= Op1; //m^2 k^2
        phip *= Op2; //m^2 k^2

        if(!RIsNull()) 
            phip *= R();
        }

    phip.mapprime(1,0);
    }

template <class Tensor>
Real inline LocalOp<Tensor>::
expect(const Tensor& phi) const
    {
    Tensor phip;
    product(phi,phip);
    return Dot(phip,phi);
    }

template <class Tensor>
Tensor inline LocalOp<Tensor>::
deltaRho(const Tensor& AA, const Tensor& combine, Direction dir) const
    {
    Tensor delta(AA);
    if(dir == Fromleft)
        {
        if(!LIsNull()) delta *= L();
        delta *= (*Op1_);
        }
    else //dir == Fromright
        {
        if(!RIsNull()) delta *= R();
        delta *= (*Op2_);
        }

    delta.noprime();
    delta = combine * delta;
    auto ci = commonIndex(combine,delta);
    
    delta *= dag(prime(delta,ci));

    return delta;
    }


template <class Tensor>
Tensor inline LocalOp<Tensor>::
diag() const
    {
    if(!(*this)) Error("LocalOp is null");

    auto& Op1 = *Op1_;
    auto& Op2 = *Op2_;

    IndexT toTie;
    bool found = false;
    for(auto& s : Op1.inds())
        {
        if(s.primeLevel() == 0 && s.type() == Site) 
            {
            toTie = s;
            found = true;
            break;
            }
        }
    if(!found) 
        {
        Print(Op1);
        Error("Couldn't find Index");
        }
    auto Diag = noprime(Op1 * diagTensor(1,toTie,prime(toTie)),toTie);

    found = false;
    for(auto& s : Op2.inds())
        {
        if(s.primeLevel() == 0 && s.type() == Site) 
            {
            toTie = s;
            found = true;
            break;
            }
        }
    if(!found) Error("Couldn't find Index");
    Diag *= noprime(Op2 * diagTensor(1,toTie,prime(toTie)),toTie);

    if(!LIsNull())
        {
        found = false;
        for(auto& ll : L().inds())
            {
            if(ll.primeLevel() == 0 && hasindex(L(),prime(ll)))
                {
                toTie = ll;
                found = true;
                break;
                }
            }
        if(found)
            Diag *= noprime(L()*diagTensor(1,toTie,prime(toTie)),toTie);
        else
            Diag *= L();
        }

    if(!RIsNull())
        {
        found = false;
        for(auto& rr : R().inds())
            {
            if(rr.primeLevel() == 0 && hasindex(R(),prime(rr)))
                {
                toTie = rr;
                found = true;
                break;
                }
            }
        if(found)
            Diag *= noprime(R()*diagTensor(1,toTie,prime(toTie)),toTie);
        else
            Diag *= R();
        }

    Diag.dag();
    Diag.takeReal(); //Diag must be real since operator assumed Hermitian
    return Diag;
    }

template <class Tensor>
long inline LocalOp<Tensor>::
size() const
    {
    if(!(*this)) Error("LocalOp is default constructed");
    if(size_ == -1)
        {
        //Calculate linear size of this 
        //op as a square matrix
        size_ = 1;
        if(!LIsNull()) 
            {
            for(auto& I : L().inds())
                {
                if(I.primeLevel() > 0)
                    {
                    size_ *= I.m();
                    break;
                    }
                }
            }
        if(!RIsNull()) 
            {
            for(auto& I : R().inds())
                {
                if(I.primeLevel() > 0)
                    {
                    size_ *= I.m();
                    break;
                    }
                }
            }

        size_ *= findtype(*Op1_,Site).m();
        size_ *= findtype(*Op2_,Site).m();
        }
    return size_;
    }

} //namespace itensor

#endif
