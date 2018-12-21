//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_LOCAL_OP
#define __ITENSOR_LOCAL_OP
#include "itensor/itensor.h"
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


class LocalOp
    {
    ITensor const* Op1_;
    ITensor const* Op2_;
    ITensor const* L_;
    ITensor const* R_;
    mutable size_t size_;
    public:


    //
    // Constructors
    //

    LocalOp();

    LocalOp(ITensor const& Op1, 
            ITensor const& Op2,
            Args const& args = Args::global());

    LocalOp(ITensor const& Op1, 
            ITensor const& Op2, 
            ITensor const& L, 
            ITensor const& R,
            Args const& args = Args::global());

    //
    // Sparse Matrix Methods
    //

    void
    product(ITensor const& phi, ITensor & phip) const;

    Real
    expect(ITensor const& phi) const;

    ITensor
    deltaRho(ITensor const& rho, 
             ITensor const& combine, 
             Direction dir) const;

    ITensor
    diag() const;

    size_t
    size() const;

    //
    // Accessor Methods
    //

    void
    update(ITensor const& Op1, ITensor const& Op2);

    void
    update(ITensor const& Op1, 
           ITensor const& Op2, 
           ITensor const& L, 
           ITensor const& R);

    ITensor const&
    Op1() const 
        { 
        if(!(*this)) Error("LocalOp is default constructed");
        return *Op1_;
        }

    ITensor const&
    Op2() const 
        { 
        if(!(*this)) Error("LocalOp is default constructed");
        return *Op2_;
        }

    ITensor const&
    L() const 
        { 
        if(!(*this)) Error("LocalOp is default constructed");
        return *L_;
        }

    ITensor const&
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

inline LocalOp::
LocalOp()
    :
    Op1_(nullptr),
    Op2_(nullptr),
    L_(nullptr),
    R_(nullptr),
    size_(-1)
    { 
    }

inline LocalOp::
LocalOp(const ITensor& Op1, const ITensor& Op2,
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

inline LocalOp::
LocalOp(const ITensor& Op1, const ITensor& Op2, 
        const ITensor& L, const ITensor& R,
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

void inline LocalOp::
update(const ITensor& Op1, const ITensor& Op2)
    {
    Op1_ = &Op1;
    Op2_ = &Op2;
    L_ = nullptr;
    R_ = nullptr;
    size_ = -1;
    }

void inline LocalOp::
update(const ITensor& Op1, const ITensor& Op2, 
       const ITensor& L, const ITensor& R)
    {
    update(Op1,Op2);
    L_ = &L;
    R_ = &R;
    }

bool inline LocalOp::
LIsNull() const
    {
    if(L_ == nullptr) return true;
    return !bool(*L_);
    }

bool inline LocalOp::
RIsNull() const
    {
    if(R_ == nullptr) return true;
    return !bool(*R_);
    }

void inline LocalOp::
product(ITensor const& phi, 
        ITensor      & phip) const
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

    phip.mapPrime(1,0);
    }

Real inline LocalOp::
expect(const ITensor& phi) const
    {
    ITensor phip;
    product(phi,phip);
    return (dag(phip) * phi).real();
    }

ITensor inline LocalOp::
deltaRho(ITensor const& AA, 
         ITensor const& combine, 
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
    drho.noPrime();
    drho = combine * drho;
    auto ci = commonIndex(combine,drho);
    drho *= dag(prime(drho,ci));

    //Expedient to ensure drho is Hermitian
    drho = drho + dag(swapPrime(drho,0,1));
    drho /= 2.;

    return drho;
    }


ITensor inline LocalOp::
diag() const
    {
    if(!(*this)) Error("LocalOp is null");

    auto& Op1 = *Op1_;
    auto& Op2 = *Op2_;

    //lambda helper function:
    auto findIndPair = [](ITensor const& T) {
        for(auto& s : T.inds())
            {
            if(s.primeLevel() == 0 && hasIndex(T,prime(s))) 
                {
                return s;
                }
            }
        return Index();
        };

    auto toTie = noPrime(findIndex(Op1,"Site"));
    auto Diag = Op1 * delta(toTie,prime(toTie),prime(toTie,2));
    Diag.noPrime();

    toTie = noPrime(findIndex(Op2,"Site"));
    auto Diag2 = Op2 * delta(toTie,prime(toTie),prime(toTie,2));
    Diag *= noPrime(Diag2);

    if(!LIsNull())
        {
        toTie = findIndPair(L());
        if(toTie)
            {
            auto DiagL = L() * delta(toTie,prime(toTie),prime(toTie,2));
            Diag *= noPrime(DiagL);
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
            Diag *= noPrime(DiagR);
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

size_t inline LocalOp::
size() const
    {
    if(!(*this)) Error("LocalOp is default constructed");
    if(size_ == size_t(-1))
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

        size_ *= findIndex(*Op1_,"Site").m();
        size_ *= findIndex(*Op2_,"Site").m();
        }
    return size_;
    }

} //namespace itensor

#endif
