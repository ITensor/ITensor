#ifndef __ITENSOR_PROJECTED_OP
#define __ITENSOR_PROJECTED_OP

class ProjectedOp
    {
    public:

    ProjectedOp(const MPO& Op, int b = 1, int num_center = 2)
        : Op_(Op),
          b_(b),
          L_(Op.NN()+1),
          R_(Op.NN()+1),
          LHlim_(1),
          RHlim_(Op.NN()),
          nc_(num_center)
        { }

    void
    product(const ITensor& phi, ITensor& phip) const;

    void
    diag(ITensor& D) const;

    void
    setBond(int b, const MPS& psi);

    const ITensor&
    L() const { return L_.at(b_); }

    const ITensor&
    R() const { return R_.at(b_+nc_-1); }

    int
    numCenter() const { return nc_; }
    void
    numCenter(int val) { nc_ = val; }

    private:

    void
    makeL(const MPS& psi, int k = -1);

    void
    makeR(const MPS& psi, int k = -1);

    const MPO& Op_;
    int b_;
    std::vector<ITensor> L_,R_;
    int LHlim_,RHlim_;
    int nc_;

    };

inline void ProjectedOp::
product(const ITensor& phi, ITensor& phip) const
    {
    phip = (L().is_null() ? phi : L() * phi); //m^3 k d
    phip *= Op_.AA(b_);   //m^2 k^2 d^2
    phip *= Op_.AA(b_+nc_-1); //m^2 k^2 d^3
    if(R().is_not_null()) 
        phip *= R();
    phip.mapprime(1,0);
    }

inline void ProjectedOp::
diag(ITensor& D) const
    {
    Index toTie;
    bool found = false;

    ITensor Diag = Op_.AA(b_);
    for(int j = 1; j <= Diag.r(); ++j)
        {
        const Index& s = Diag.index(j);
        if(s.primeLevel() == 0 && s.type() == Site) 
            {
            toTie = s;
            found = true;
            break;
            }
        }
    if(!found) Error("Couldn't find Index");
    Diag.tieIndices(toTie,primed(toTie),toTie);

    const ITensor& Op2 = Op_.AA(b_+nc_-1);
    found = false;
    for(int j = 1; j <= Op2.r(); ++j)
        {
        const Index& s = Op2.index(j);
        if(s.primeLevel() == 0 && s.type() == Site) 
            {
            toTie = s;
            found = true;
            break;
            }
        }
    if(!found) Error("Couldn't find Index");
    Diag *= tieIndices(toTie,primed(toTie),toTie,Op2);

    if(L().is_not_null())
        {
        found = false;
        for(int j = 1; j <= L().r(); ++j)
            {
            const Index& ll = L().index(j);
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

    if(R().is_not_null())
        {
        found = false;
        for(int j = 1; j <= R().r(); ++j)
            {
            const Index& ll = R().index(j);
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

inline void ProjectedOp::
setBond(int b, const MPS& psi)
    {
    makeL(psi,b);
    makeR(psi,b+nc_-1);
    b_ = b;
    LHlim_ = b_;
    RHlim_ = b_+nc_-1;
    }

inline void ProjectedOp::
makeL(const MPS& psi, int k)
    {
    if(k == -1) k = b_;
    while(LHlim_ < k)
        {
        const int j = LHlim_;
        psi.projectOp(j,Fromleft,L_[j],Op_.AA(j),L_[j+1]);
        ++LHlim_;
        }
    }

inline void ProjectedOp::
makeR(const MPS& psi, int k)
    {
    if(k == -1) k = b_+nc_-1;
    while(RHlim_ > k)
        {
        const int j = RHlim_;
        psi.projectOp(j,Fromright,R_[j],Op_.AA(j),R_[j-1]);
        --RHlim_;
        }
    }


#endif
