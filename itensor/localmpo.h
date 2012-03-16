//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_LOCALMPO
#define __ITENSOR_LOCALMPO
#include "mpo.h"
#include "localop.h"

//
// The LocalMPO class projects an MPO 
// into the reduced Hilbert space of
// some number of sites of an MPS.
// (The default is 2 sites.)
//
//   .----..---             ----...-.
//   |  |    |    |      |   |      | 
//   W1 W2.. Wj-1 Wj - Wj+1  Wj+2..WN
//   |  |    |    |      |   |      | 
//   '----..---             ----...-'
//
// 
//  Here the W's are the site tensors
//  of the MPO "Op" and the method position(j,psi)
//  has been called using the MPS 'psi' as a basis 
//  for the projection.
//
//  This results in an unprojected region of
//  num_center sites starting at site j.
//

template <class Tensor>
class LocalMPO
    {
    public:

    LocalMPO();

    LocalMPO(const MPOt<Tensor>& Op, int num_center = 2);

    typedef typename Tensor::IndexT
    IndexT;

    typedef typename Tensor::CombinerT
    CombinerT;

    void
    product(const Tensor& phi, Tensor& phip) const { lop_.product(phi,phip); }

    Real
    expect(const Tensor& phi) const { return lop_.expect(phi); }

    Tensor
    deltaRho(const Tensor& rho, 
             const CombinerT& comb, Direction dir) const
        { return lop_.deltaRho(rho,comb,dir); }

    Tensor
    deltaPhi(const Tensor& phi) const { return lop_.deltaPhi(phi); }

    void
    diag(Tensor& D) const { lop_.diag(D); }

    template <class MPSType>
    void
    position(int b, const MPSType& psi);

    const Tensor&
    L() const 
        { 
        return PH_[LHlim_];
        }

    const Tensor&
    R() const 
        { 
        return PH_[RHlim_];
        }

    const Tensor&
    bondTensor() const { return lop_.bondTensor(); }

    bool
    combineMPO() const { return lop_.combineMPO(); }
    void
    combineMPO(bool val) { lop_.combineMPO(val); }

    int
    numCenter() const { return nc_; }
    void
    numCenter(int val) 
        { 
        if(val < 1) Error("numCenter must be set >= 1");
        nc_ = val; 
        }

    int
    size() const { return lop_.size(); }

    bool
    isNull() const { return Op_ == 0; }
    bool
    isNotNull() const { return Op_ != 0; }

    static LocalMPO& Null()
        {
        static LocalMPO Null_;
        return Null_;
        }

    private:

    /////////////////
    //
    // Data Members
    //

    const MPOt<Tensor>* Op_;
    std::vector<Tensor> PH_;
    int LHlim_,RHlim_;
    int nc_;

    LocalOp<Tensor> lop_;

    //
    /////////////////

    template <class MPSType>
    void
    makeL(const MPSType& psi, int k);

    template <class MPSType>
    void
    makeR(const MPSType& psi, int k);


    };

template <class Tensor>
inline LocalMPO<Tensor>::
LocalMPO()
    : Op_(0),
      LHlim_(-1),
      RHlim_(-1),
      nc_(2)
    { }

template <class Tensor>
inline LocalMPO<Tensor>::
LocalMPO(const MPOt<Tensor>& Op, int num_center)
    : Op_(&Op),
      PH_(Op.NN()+2),
      LHlim_(0),
      RHlim_(Op.NN()+1),
      nc_(2)
    { 
    numCenter(num_center);
    }

template <class Tensor>
template <class MPSType> 
inline void LocalMPO<Tensor>::
position(int b, const MPSType& psi)
    {
    if(this->isNull()) Error("LocalMPO is null");

    makeL(psi,b-1);
    makeR(psi,b+nc_);

    LHlim_ = b-1; //not redundant since LHlim_ could be > b-1
    RHlim_ = b+nc_; //not redundant since RHlim_ could be < b+nc_

#ifdef DEBUG
    if(nc_ != 2)
        {
        Error("LocalOp only supports 2 center sites currently");
        }
#endif

    lop_.update(Op_->AA(b),Op_->AA(b+1),L(),R());
    }

template <class Tensor>
template <class MPSType> 
inline void LocalMPO<Tensor>::
makeL(const MPSType& psi, int k)
    {
    if(!PH_.empty())
    while(LHlim_ < k)
        {
        const int ll = LHlim_;
        psi.projectOp(ll+1,Fromleft,PH_.at(ll),Op_->AA(ll+1),PH_.at(ll+1));
        ++LHlim_;
        }
    }

template <class Tensor>
template <class MPSType> 
inline void LocalMPO<Tensor>::
makeR(const MPSType& psi, int k)
    {
    if(!PH_.empty())
    while(RHlim_ > k)
        {
        const int rl = RHlim_;
        psi.projectOp(rl-1,Fromright,PH_.at(rl),Op_->AA(rl-1),PH_.at(rl-1));
        --RHlim_;
        }
    }

#endif
