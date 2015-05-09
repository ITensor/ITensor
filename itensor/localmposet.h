//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_LOCALMPOSET
#define __ITENSOR_LOCALMPOSET
#include "localmpo.h"

namespace itensor {

template <class Tensor>
class LocalMPOSet
    {
    public:

    LocalMPOSet();

    LocalMPOSet(const std::vector<MPOt<Tensor> >& Op,
                const Args& args = Global::args());

    using IndexT = typename Tensor::IndexT;

    using LocalMPOT = LocalMPO<Tensor>;

    void
    product(const Tensor& phi, Tensor& phip) const;

    Real
    expect(const Tensor& phi) const;

    Tensor
    deltaRho(const Tensor& AA, 
             const Tensor& comb, Direction dir) const;

    Tensor
    diag() const;

    template <class MPSType>
    void
    position(int b, const MPSType& psi);

    bool
    combineMPO() const { return lmpo_.front().combineMPO(); }
    void
    combineMPO(bool val);

    int
    numCenter() const { return lmpo_.front().numCenter(); }
    void
    numCenter(int val);

    int
    size() const { return lmpo_.front().size(); }

    bool
    isNull() const { return Op_ == 0; }

    bool
    doWrite() const { return false; }
    void
    doWrite(bool val) 
        { 
        if(val) Error("Write to disk not yet supported LocalMPOSet");
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
LocalMPOSet(const std::vector<MPOt<Tensor> >& Op,
            const Args& args)
    : 
    Op_(&Op),
    lmpo_(Op.size())
    { 
    for(size_t n = 0; n < lmpo_.size(); ++n)
        {
        lmpo_[n] = LocalMPOT(Op.at(n));
        }
    }

template <class Tensor>
void inline LocalMPOSet<Tensor>::
product(const Tensor& phi, Tensor& phip) const
    {
    lmpo_.front().product(phi,phip);

    Tensor phi_n;
    for(size_t n = 1; n < lmpo_.size(); ++n)
        {
        lmpo_[n].product(phi,phi_n);
        phip += phi_n;
        }
    }

template <class Tensor>
Real inline LocalMPOSet<Tensor>::
expect(const Tensor& phi) const
    {
    Real ex_ = 0;
    for(size_t n = 0; n < lmpo_.size(); ++n)
        {
        ex_ += lmpo_[n].expect(phi);
        }
    return ex_;
    }

template <class Tensor>
Tensor inline LocalMPOSet<Tensor>::
deltaRho(const Tensor& AA,
         const Tensor& comb, Direction dir) const
    {
    Tensor delta = lmpo_.front().deltaRho(AA,comb,dir);
    for(size_t n = 1; n < lmpo_.size(); ++n)
        {
        delta += lmpo_[n].deltaRho(AA,comb,dir);
        }
    return delta;
    }

template <class Tensor>
Tensor inline LocalMPOSet<Tensor>::
diag() const
    {
    Tensor D = lmpo_.front().diag();
    for(size_t n = 1; n < lmpo_.size(); ++n)
        {
        D += lmpo_[n].diag();
        }
    return D;
    }

template <class Tensor>
template <class MPSType> 
void inline LocalMPOSet<Tensor>::
position(int b, const MPSType& psi)
    {
    for(size_t n = 0; n < lmpo_.size(); ++n)
        lmpo_[n].position(b,psi);
    }

template <class Tensor>
void inline LocalMPOSet<Tensor>::
combineMPO(bool val)
    {
    for(size_t n = 0; n < lmpo_.size(); ++n)
        lmpo_[n].combineMPO(val);
    }

template <class Tensor>
void inline LocalMPOSet<Tensor>::
numCenter(int val)
    {
    for(size_t n = 0; n < lmpo_.size(); ++n)
        lmpo_[n].numCenter(val);
    }

} //namespace itensor

#endif
