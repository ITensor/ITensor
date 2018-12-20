//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_LOCALMPO_MPS
#define __ITENSOR_LOCALMPO_MPS
#include "itensor/mps/localmpo.h"

namespace itensor {

template <class Tensor>
class LocalMPO_MPS
    {
    public:
    using LocalMPOType = LocalMPO<Tensor>;
    private:
    MPOt<Tensor> const* Op_ = nullptr;
    std::vector<MPSt<Tensor>> const* psis_ = nullptr;
    //LocalMPO object representing projected version
    //of the MPO Op_
    LocalMPOType lmpo_; 
    //LocalMPO objects representing projected version
    //of each MPS in psis_
    std::vector<LocalMPOType> lmps_;
    Real weight_ = 1;
    public:

    LocalMPO_MPS() { }

    LocalMPO_MPS(MPOt<Tensor> const& Op, 
                 std::vector<MPSt<Tensor> > const& psis,
                 Args const& args = Args::global());

    LocalMPO_MPS(MPOt<Tensor> const& Op, 
                 Tensor const& LOp,
                 Tensor const& ROp,
                 std::vector<MPSt<Tensor>> const& psis,
                 std::vector<Tensor> const& Lpsi,
                 std::vector<Tensor> const& Rpsi,
                 Args const& args = Args::global());


    //
    // Sparse matrix methods
    //

    void
    product(Tensor const& phi, 
            Tensor& phip) const;

    Real
    expect(Tensor const& phi) const { return lmpo_.expect(phi); }

    Tensor
    deltaRho(Tensor const& AA, 
             Tensor const& comb, 
             Direction dir) const
        { return lmpo_.deltaRho(AA,comb,dir); }

    Tensor
    diag() const { return lmpo_.diag(); }

    template <class MPSType>
    void
    position(int b, MPSType const& psi);

    size_t
    size() const { return lmpo_.size(); }

    explicit
    operator bool() const { return bool(Op_); }

    Real
    weight() const { return weight_; }
    void
    weight(Real val) { weight_ = val; }

    bool
    doWrite() const { return lmpo_.doWrite(); }
    void
    doWrite(bool val, Args const& args = Args::global()) { lmpo_.doWrite(val,args); }

    };

template <class Tensor>
inline LocalMPO_MPS<Tensor>::
LocalMPO_MPS(MPOt<Tensor> const& Op,
             std::vector<MPSt<Tensor>> const& psis,
             Args const& args)
  : Op_(&Op),
    psis_(&psis),
    lmps_(psis.size()),
    weight_(args.getReal("Weight",1))
    { 
    lmpo_ = LocalMPOType(Op);

    for(auto j : range(lmps_.size()))
        {
        lmps_[j] = LocalMPOType(psis[j]);
        }
    }

template <class Tensor>
inline LocalMPO_MPS<Tensor>::
LocalMPO_MPS(MPOt<Tensor> const& Op, 
             Tensor const& LOp,
             Tensor const& ROp,
             std::vector<MPSt<Tensor>> const& psis,
             std::vector<Tensor> const& Lpsi,
             std::vector<Tensor> const& Rpsi,
             Args const& args)
  : Op_(&Op),
    psis_(&psis),
    lmps_(psis.size()),
    weight_(args.getReal("Weight",1))
    { 
    lmpo_ = LocalMPOType(Op,LOp,ROp);
#ifdef DEBUG
    if(Lpsi.size() != psis.size()) Error("Lpsi must have same number of elements as psis");
    if(Rpsi.size() != psis.size()) Error("Rpsi must have same number of elements as psis");
#endif

    for(auto j : range(lmps_.size()))
        {
        lmps_[j] = LocalMPOType(psis[j],Lpsi[j],Rpsi[j]);
        }
    }

template <class Tensor>
void inline LocalMPO_MPS<Tensor>::
product(Tensor const& phi, 
        Tensor & phip) const
    {
    lmpo_.product(phi,phip);

    Tensor outer;
    for(auto& M : lmps_)
        {
        M.product(phi,outer);
        outer *= weight_;
        phip += outer;
        }
    }

template <class Tensor>
template <class MPSType> 
void inline LocalMPO_MPS<Tensor>::
position(int b, const MPSType& psi)
    {
    lmpo_.position(b,psi);
    for(auto& M : lmps_)
        {
        M.position(b,psi);
        }
    }

} //namespace itensor

#endif
