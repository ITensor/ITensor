#ifndef __ITENSOR_PROJECTED_OP
#define __ITENSOR_PROJECTED_OP

template <class Tensor>
class ProjectedOp
    {
    public:

    ProjectedOp(const MPOt<Tensor>& Op, int num_center = 2);

    void
    product(const Tensor& phi, Tensor& phip) const;

    void
    diag(Tensor& D) const;

    void
    setBond(int b, const MPSt<Tensor>& psi);

    const Tensor&
    L() const { return L_.at(LHlim_); }

    const Tensor&
    R() const { return R_.at(RHlim_); }

    int
    numCenter() const { return nc_; }
    void
    numCenter(int val) { nc_ = val; }

    int
    size() const { return size_; }

    typedef typename Tensor::IndexT
    IndexT;

    private:

    void
    makeL(const MPSt<Tensor>& psi, int k);

    void
    makeR(const MPSt<Tensor>& psi, int k);

    const MPOt<Tensor>& Op_;
    std::vector<Tensor> L_,R_;
    int LHlim_,RHlim_;
    int nc_;
    int size_;

    };

template <class Tensor>
inline ProjectedOp<Tensor>::
ProjectedOp(const MPOt<Tensor>& Op, int num_center)
    : Op_(Op),
      L_(Op.NN()+1),
      R_(Op.NN()+1),
      LHlim_(1),
      RHlim_(Op.NN()),
      nc_(num_center),
      size_(-1)
    { }

template <class Tensor>
inline void ProjectedOp<Tensor>::
product(const Tensor& phi, Tensor& phip) const
    {
    if(L().isNull())
        {
        phip = phi;
        if(R().isNotNull()) 
            phip *= R(); //m^3 k d
        for(int j = RHlim_; j >= LHlim_; --j)
            phip *= Op_.AA(j); //m^2 k^2
        }
    else
        {
        phip = phi * L(); //m^3 k d
        for(int j = LHlim_; j <= RHlim_; ++j)
            phip *= Op_.AA(j); //m^2 k^2
        if(R().isNotNull()) 
            phip *= R();
        }
    phip.mapprime(1,0);
    }

template <class Tensor>
inline void ProjectedOp<Tensor>::
diag(Tensor& D) const
    {
    IndexT toTie;
    bool found = false;

    Tensor Diag = Op_.AA(LHlim_);
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

    const Tensor& Op2 = Op_.AA(RHlim_);
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
inline void ProjectedOp<Tensor>::
setBond(int b, const MPSt<Tensor>& psi)
    {
    makeL(psi,b);
    makeR(psi,b+nc_-1);
    LHlim_ = b; //not redundant since LHlim_ could be > b
    RHlim_ = b+nc_-1; //not redundant since RHlim_ could be < b+nc_-1

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
inline void ProjectedOp<Tensor>::
makeL(const MPSt<Tensor>& psi, int k)
    {
    while(LHlim_ < k)
        {
        const int ll = LHlim_;
        psi.projectOp(ll,Fromleft,L_.at(ll),Op_.AA(ll),L_.at(ll+1));
        ++LHlim_;
        }
    }

template <class Tensor>
inline void ProjectedOp<Tensor>::
makeR(const MPSt<Tensor>& psi, int k)
    {
    while(RHlim_ > k)
        {
        const int rl = RHlim_;
        psi.projectOp(rl,Fromright,R_.at(rl),Op_.AA(rl),R_.at(rl-1));
        --RHlim_;
        }
    }


#endif
