//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_LOCAL_OP
#define __ITENSOR_LOCAL_OP
#include "itensor/iqtensor.h"
//#include "itensor/util/print_macro.h"

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
    Tensor const* Op1_;
    Tensor const* Op2_;
    Tensor const* L_;
    Tensor const* R_;
    mutable long size_;
    public:

    using IndexT = typename Tensor::index_type;

    //
    // Constructors
    //

    LocalOp();

    LocalOp(Tensor const& Op1, 
            Tensor const& Op2,
            Args const& args = Global::args());

    LocalOp(Tensor const& Op1, 
            Tensor const& Op2, 
            Tensor const& L, 
            Tensor const& R,
            Args const& args = Global::args());

    //
    // Sparse Matrix Methods
    //

    void
    product(Tensor const& phi, Tensor & phip) const;

    Real
    expect(Tensor const& phi) const;

    Tensor
    deltaRho(Tensor const& rho, 
             Tensor const& combine, 
             Direction dir) const;

    Tensor
    diag() const;

    long
    size() const;

    //
    // Accessor Methods
    //

    void
    update(Tensor const& Op1, Tensor const& Op2);

    void
    update(Tensor const& Op1, 
           Tensor const& Op2, 
           Tensor const& L, 
           Tensor const& R);

    Tensor const&
    Op1() const 
        { 
        if(!(*this)) Error("LocalOp is default constructed");
        return *Op1_;
        }

    Tensor const&
    Op2() const 
        { 
        if(!(*this)) Error("LocalOp is default constructed");
        return *Op2_;
        }

    Tensor const&
    L() const 
        { 
        if(!(*this)) Error("LocalOp is default constructed");
        return *L_;
        }

    Tensor const&
    R() const 
        { 
        if(!(*this)) Error("LocalOp is default constructed");
        return *R_;
        }

    explicit operator bool() const { return bool(Op1_); }

    bool
    LIsNull() const;

    bool
    RIsNull() const;

    };

template <class Tensor>
inline LocalOp<Tensor>::
LocalOp()
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
    return !bool(*R_);
    }

template <class Tensor>
void inline LocalOp<Tensor>::
product(Tensor const& phi, 
        Tensor      & phip) const
    {
    if(!(*this)) Error("LocalOp is null");

    auto& Op1 = *Op1_;
    auto& Op2 = *Op2_;

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
    return (dag(phip) * phi).real();
    }

template <class Tensor>
Tensor inline LocalOp<Tensor>::
deltaRho(Tensor const& AA, 
         Tensor const& combine, 
         Direction dir) const
    {
    auto drho = AA;
    if(dir == Fromleft)
        {
        if(!LIsNull()) drho *= L();
        drho *= (*Op1_);
        }
    else //dir == Fromright
        {
        if(!RIsNull()) drho *= R();
        drho *= (*Op2_);
        }
    drho.noprime();
    drho = combine * drho;
    auto ci = commonIndex(combine,drho);
    drho *= dag(prime(drho,ci));

    //Expedient to ensure drho is Hermitian
    drho = drho + dag(swapPrime(drho,0,1));
    drho /= 2.;

    return drho;
    }


template <class Tensor>
Tensor inline LocalOp<Tensor>::
diag() const
    {
    if(!(*this)) Error("LocalOp is null");

    auto& Op1 = *Op1_;
    auto& Op2 = *Op2_;

    //lambda helper function:
    auto findIndPair = [](Tensor const& T) {
        for(auto& s : T.inds())
            {
            if(s.primeLevel() == 0 && hasindex(T,prime(s))) 
                {
                return s;
                }
            }
        return IndexT();
        };

    auto toTie = noprime(findtype(Op1,Site));
    auto Diag = Op1 * delta(toTie,prime(toTie),prime(toTie,2));
    Diag.noprime();

    toTie = noprime(findtype(Op2,Site));
    auto Diag2 = Op2 * delta(toTie,prime(toTie),prime(toTie,2));
    Diag *= noprime(Diag2);

    if(!LIsNull())
        {
        toTie = findIndPair(L());
        if(toTie)
            {
            auto DiagL = L() * delta(toTie,prime(toTie),prime(toTie,2));
            Diag *= noprime(DiagL);
            }
        else
            {
            Diag *= L();
            }
        }

    if(!RIsNull())
        {
        toTie = findIndPair(R());
        if(toTie)
            {
            auto DiagR = R() * delta(toTie,prime(toTie),prime(toTie,2));
            Diag *= noprime(DiagR);
            }
        else
            {
            Diag *= R();
            }
        }

    Diag.dag();
    //Diag must be real since operator assumed Hermitian
    Diag.takeReal();

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
