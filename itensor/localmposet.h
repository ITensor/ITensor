//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_LOCALMPOSET
#define __ITENSOR_LOCALMPOSET
#include "mpo.h"
#include "localmpo.h"

template <class Tensor>
class LocalMPOSet
    {
    public:

    LocalMPOSet();

    LocalMPOSet(const std::vector<MPOt<Tensor> >& Op);

    typedef typename Tensor::IndexT
    IndexT;

    typedef typename Tensor::CombinerT
    CombinerT;

    typedef LocalMPO<Tensor>
    LocalMPOT;

    void
    product(const Tensor& phi, Tensor& phip) const;

    Real
    expect(const Tensor& phi) const;

    Tensor
    deltaRho(const Tensor& AA, 
             const CombinerT& comb, Direction dir) const;

    Tensor
    diag() const;

    template <class MPSType>
    void
    position(int b, const MPSType& psi);

    bool
    combineMPO() const { return lmpo_.at(1).combineMPO(); }
    void
    combineMPO(bool val);

    int
    numCenter() const { return lmpo_.at(1).numCenter(); }
    void
    numCenter(int val);

    int
    size() const { return lmpo_.at(1).size(); }

    bool
    isNull() const { return Op_ == 0; }
    bool
    isNotNull() const { return Op_ != 0; }

    static LocalMPOSet& Null()
        {
        static LocalMPOSet Null_;
        return Null_;
        }

    private:

    /////////////////
    //
    // Data Members
    //

    const std::vector<MPOt<Tensor> >* Op_;
    std::vector<LocalMPO<Tensor> > lmpo_;

    //
    /////////////////

    };

template <class Tensor>
inline LocalMPOSet<Tensor>::
LocalMPOSet()
    : 
    Op_(0)
    { }

template <class Tensor>
inline LocalMPOSet<Tensor>::
LocalMPOSet(const std::vector<MPOt<Tensor> >& Op)
    : 
    Op_(&Op),
    lmpo_(Op.size())
    { 
    for(size_t n = 1; n < lmpo_.size(); ++n)
        {
        lmpo_[n] = LocalMPOT(Op.at(n));
        }
    }

template <class Tensor>
void inline LocalMPOSet<Tensor>::
product(const Tensor& phi, Tensor& phip) const
    {
    lmpo_.at(1).product(phi,phip);

    Tensor phi_n;
    for(size_t n = 2; n < lmpo_.size(); ++n)
        {
        lmpo_.at(n).product(phi,phi_n);
        phip += phi_n;
        }
    }

template <class Tensor>
Real inline LocalMPOSet<Tensor>::
expect(const Tensor& phi) const
    {
    Real ex_ = 0;
    for(size_t n = 1; n < lmpo_.size(); ++n)
        {
        ex_ += lmpo_.at(n).expect(phi);
        }
    return ex_;
    }

template <class Tensor>
Tensor inline LocalMPOSet<Tensor>::
deltaRho(const Tensor& AA,
         const CombinerT& comb, Direction dir) const
    {
    Tensor delta = lmpo_.at(1).deltaRho(AA,comb,dir);
    for(size_t n = 2; n < lmpo_.size(); ++n)
        {
        delta += lmpo_.at(n).deltaRho(AA,comb,dir);
        }
    return delta;
    }

template <class Tensor>
Tensor inline LocalMPOSet<Tensor>::
diag() const
    {
    Tensor D = lmpo_.at(1).diag();
    for(size_t n = 2; n < lmpo_.size(); ++n)
        {
        D += lmpo_.at(n).diag();
        }
    return D;
    }

template <class Tensor>
template <class MPSType> 
void inline LocalMPOSet<Tensor>::
position(int b, const MPSType& psi)
    {
    for(size_t n = 1; n < lmpo_.size(); ++n)
        lmpo_[n].position(b,psi);
    }

template <class Tensor>
void inline LocalMPOSet<Tensor>::
combineMPO(bool val)
    {
    for(size_t n = 1; n < lmpo_.size(); ++n)
        lmpo_[n].combineMPO(val);
    }

template <class Tensor>
void inline LocalMPOSet<Tensor>::
numCenter(int val)
    {
    for(size_t n = 1; n < lmpo_.size(); ++n)
        lmpo_[n].numCenter(val);
    }

#endif
