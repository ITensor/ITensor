//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_PROJECTED_OP
#define __ITENSOR_PROJECTED_OP
#include "mpo.h"

template <class Tensor>
class ProjectedOp
    {
    public:

    ProjectedOp();

    ProjectedOp(const MPOt<Tensor>& Op, int num_center = 2);

    ProjectedOp(const MPOt<Tensor>& Op, const Tensor& L, const Tensor& R, int num_center = 2);

    typedef typename Tensor::IndexT
    IndexT;

    typedef typename Tensor::CombinerT
    CombinerT;

    void
    product(const Tensor& phi, Tensor& phip) const;

    Real
    expect(const Tensor& phi) const;

    Tensor
    deltaRho(const Tensor& rho, 
             const CombinerT& comb, Direction dir) const;

    void
    diag(Tensor& D) const;

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

    bool
    combineMPO() const { return combine_mpo_; }
    void
    combineMPO(bool val) { combine_mpo_ = val; }

    int
    numCenter() const { return nc_; }
    void
    numCenter(int val) { nc_ = val; }

    int
    size() const { return size_; }

    bool
    isNull() const { return Op_ == 0; }
    bool
    isNotNull() const { return Op_ != 0; }

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
    int size_;
    bool combine_mpo_;
    Tensor mpoh_;

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
inline ProjectedOp<Tensor>::
ProjectedOp()
    : Op_(0),
      LHlim_(0),
      RHlim_(0),
      pL_(0),
      pR_(0),
      nc_(2),
      size_(-1),
      combine_mpo_(true)
    { }

template <class Tensor>
inline ProjectedOp<Tensor>::
ProjectedOp(const MPOt<Tensor>& Op, int num_center)
    : Op_(&Op),
      L_(Op.NN()+1),
      R_(Op.NN()+1),
      LHlim_(1),
      RHlim_(Op.NN()),
      pL_(&(L_[LHlim_])),
      pR_(&(R_[RHlim_])),
      nc_(num_center),
      size_(-1),
      combine_mpo_(true)
    { }

template <class Tensor>
inline ProjectedOp<Tensor>::
ProjectedOp(const MPOt<Tensor>& Op, 
            const Tensor& L, const Tensor& R, int num_center)
    : Op_(&Op),
      LHlim_(0),
      RHlim_(0),
      nc_(num_center),
      size_(-1),
      combine_mpo_(true)
    { 
    pL_ = &L;
    pR_ = &R;
    }

template <class Tensor>
inline void ProjectedOp<Tensor>::
product(const Tensor& phi, Tensor& phip) const
    {
    if(this->isNull()) Error("ProjectedOp is null");

    if(L().isNull())
        {
        phip = phi;
        if(R().isNotNull()) 
            phip *= R(); //m^3 k d

        if(combine_mpo_)
            {
            phip *= mpoh_;
            }
        else
            {
            for(int j = RHlim_; j >= LHlim_; --j)
                phip *= Op_->AA(j); //m^2 k^2
            }
        }
    else
        {
        phip = phi * L(); //m^3 k d
        if(combine_mpo_)
            {
            phip *= mpoh_;
            }
        else
            {
        for(int j = LHlim_; j <= RHlim_; ++j)
            phip *= Op_->AA(j); //m^2 k^2
            }
        if(R().isNotNull()) 
            phip *= R();
        }
    phip.mapprime(1,0);
    }

template <class Tensor>
inline Real ProjectedOp<Tensor>::
expect(const Tensor& phi) const
    {
    Tensor phip;
    product(phi,phip);
    return Dot(phip,phi);
    }

template <class Tensor>
inline Tensor ProjectedOp<Tensor>::
deltaRho(const Tensor& rho, const CombinerT& comb, Direction dir) const
    {
    Tensor delta(rho);

    Tensor A;
    IndexT hl;
    if(dir == Fromleft)
        {
        A = Op_->AA(LHlim_);
        if(L().isNotNull()) A *= L();
        hl = Op_->RightLinkInd(LHlim_);
        }
    else //dir == Fromright
        {
        A = Op_->AA(RHlim_);
        if(R().isNotNull()) A *= R();
        hl = Op_->LeftLinkInd(RHlim_);
        }

    A = conj(comb) * A;
    A = primed(comb) * A;

    A.mapprime(1,2);
    delta *= A;
    delta.mapprime(2,0);

    A.conj();
    A.mapprime(0,1);
    delta *= A;
    delta.mapprime(2,1);

    delta.trace(hl,primed(hl));

    return delta;
    }

template <class Tensor>
inline void ProjectedOp<Tensor>::
diag(Tensor& D) const
    {
    if(this->isNull()) Error("ProjectedOp is null");

    IndexT toTie;
    bool found = false;

    Tensor Diag = Op_->AA(LHlim_);
    for(int j = 1; j <= Diag.r(); ++j)
        {
        const IndexT& s = Diag.index(j);
        if(s.primeLevel() == 0 && s.type() == Site) 
            {
            toTie = s;
            found = true;
            break;
            }
        }
    if(!found) Error("Couldn't find Index");
    Diag.tieIndices(toTie,primed(toTie),toTie);

    const Tensor& Op2 = Op_->AA(RHlim_);
    found = false;
    for(int j = 1; j <= Op2.r(); ++j)
        {
        const IndexT& s = Op2.index(j);
        if(s.primeLevel() == 0 && s.type() == Site) 
            {
            toTie = s;
            found = true;
            break;
            }
        }
    if(!found) Error("Couldn't find Index");
    Diag *= tieIndices(toTie,primed(toTie),toTie,Op2);

    if(L().isNotNull())
        {
        found = false;
        for(int j = 1; j <= L().r(); ++j)
            {
            const IndexT& ll = L().index(j);
            if(ll.primeLevel() == 0 && !Diag.hasindex(ll))
                {
                toTie = ll;
                found = true;
                break;
                }
            }
        if(found)
            Diag *= tieIndices(toTie,primed(toTie),toTie,L());
        else
            Diag *= L();
        }

    if(R().isNotNull())
        {
        found = false;
        for(int j = 1; j <= R().r(); ++j)
            {
            const IndexT& ll = R().index(j);
            if(ll.primeLevel() == 0 && !Diag.hasindex(ll))
                {
                toTie = ll;
                found = true;
                break;
                }
            }
        if(found)
            Diag *= tieIndices(toTie,primed(toTie),toTie,R());
        else
            Diag *= R();
        }

    D.assignFrom(Diag);
    }

template <class Tensor>
template <class MPSType> 
inline void ProjectedOp<Tensor>::
position(int b, const MPSType& psi)
    {
    if(this->isNull()) Error("ProjectedOp is null");

    makeL(psi,b);
    makeR(psi,b+nc_-1);

    LHlim_ = b; //not redundant since LHlim_ could be > b
    pL_ = &(L_.at(LHlim_));

    RHlim_ = b+nc_-1; //not redundant since RHlim_ could be < b+nc_-1
    pR_ = &(R_.at(RHlim_));

    if(combine_mpo_)
        {
        if(nc_ != 2) Error("nc_ must be 2 for combine_mpo_");
        mpoh_ = Op_->AA(b)*Op_->AA(b+1);
        }

    //Calculate linear size of this projected
    //op as a square matrix
    size_ = 1;
    if(L().isNotNull()) 
        {
        size_ *= index_in_common(psi.AA(LHlim_),L(),Link).m();
        }
    if(R().isNotNull()) 
        {
        size_ *= index_in_common(psi.AA(RHlim_),R(),Link).m();
        }
    for(int j = LHlim_; j <= RHlim_; ++j)
        {
        size_ *= psi.AA(j).findtype(Site).m();
        }
    }

template <class Tensor>
template <class MPSType> 
inline void ProjectedOp<Tensor>::
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
inline void ProjectedOp<Tensor>::
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
