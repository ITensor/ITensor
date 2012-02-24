//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_LOCALMPO
#define __ITENSOR_LOCALMPO
#include "mpo.h"
#include "localop.h"

template <class Tensor>
class LocalMPO
    {
    public:

    LocalMPO();

    LocalMPO(const MPOt<Tensor>& Op, int num_center = 2);

    LocalMPO(const MPOt<Tensor>& Op, const Tensor& L, const Tensor& R, int num_center = 2);

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

    void
    diag(Tensor& D) const { lop_.diag(D); }

    template <class MPSType>
    void
    position(int b, const MPSType& psi);

    const Tensor&
    L() const 
        { 
        if(pL_ == 0) Error("pL_ not set");
        return *pL_; 
        }

    const Tensor&
    R() const 
        { 
        if(pR_ == 0) Error("pR_ not set");
        return *pR_; 
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
    numCenter(int val) { nc_ = val; }

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
    std::vector<Tensor> L_,R_;
    int LHlim_,RHlim_;
    const Tensor *pL_, *pR_;
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
      LHlim_(0),
      RHlim_(0),
      pL_(0),
      pR_(0),
      nc_(2)
    { }

template <class Tensor>
inline LocalMPO<Tensor>::
LocalMPO(const MPOt<Tensor>& Op, int num_center)
    : Op_(&Op),
      L_(Op.NN()+1),
      R_(Op.NN()+1),
      LHlim_(1),
      RHlim_(Op.NN()),
      pL_(&(L_[LHlim_])),
      pR_(&(R_[RHlim_])),
      nc_(num_center)
    { 
    }

template <class Tensor>
inline LocalMPO<Tensor>::
LocalMPO(const MPOt<Tensor>& Op, 
            const Tensor& L, const Tensor& R, int num_center)
    : Op_(&Op),
      LHlim_(0),
      RHlim_(0),
      nc_(num_center)
    { 
    pL_ = &L;
    pR_ = &R;
    }

template <class Tensor>
template <class MPSType> 
inline void LocalMPO<Tensor>::
position(int b, const MPSType& psi)
    {
    if(this->isNull()) Error("LocalMPO is null");

    makeL(psi,b);
    makeR(psi,b+nc_-1);

    LHlim_ = b; //not redundant since LHlim_ could be > b
    pL_ = &(L_.at(LHlim_));

    RHlim_ = b+nc_-1; //not redundant since RHlim_ could be < b+nc_-1
    pR_ = &(R_.at(RHlim_));

    lop_ = LocalOp<Tensor>(Op_->AA(b),Op_->AA(b+1),*pL_,*pR_);
    }

template <class Tensor>
template <class MPSType> 
inline void LocalMPO<Tensor>::
makeL(const MPSType& psi, int k)
    {
    if(!L_.empty())
    while(LHlim_ < k)
        {
        const int ll = LHlim_;
        //std::cout << boost::format("Shifting L from %d to %d") % ll % (ll+1) << std::endl;
        psi.projectOp(ll,Fromleft,L_.at(ll),Op_->AA(ll),L_.at(ll+1));
        ++LHlim_;
        pL_ = &(L_.at(LHlim_));
        }
    }

template <class Tensor>
template <class MPSType> 
inline void LocalMPO<Tensor>::
makeR(const MPSType& psi, int k)
    {
    if(!R_.empty())
    while(RHlim_ > k)
        {
        const int rl = RHlim_;
        psi.projectOp(rl,Fromright,R_.at(rl),Op_->AA(rl),R_.at(rl-1));
        --RHlim_;
        pR_ = &(R_.at(RHlim_));
        }
    }

#endif
