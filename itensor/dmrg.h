#ifndef ALGORITHM_H
#define ALGORITHM_H
#include "mpo.h"
#include "sparse.h"
#include "davidson.h"
#include "Sweeps.h"
#include "DMRGOpts.h"

template<class Tensor,class TensorSet, class OpTensorSet>
void applyProjOp(const Tensor& phi, const TensorSet& L, const TensorSet& R, const OpTensorSet& H, Tensor& Hphi)
{
    bool useL(L.size() == 0 ? false : L[0].is_not_null()),
         useR(R.size() == 0 ? false : R[0].is_not_null());
    Hphi = (useL ? L[0]*phi : phi); 
    Hphi *= H[0];
    if(useR) Hphi *= R[0];
    for(unsigned int j = 1; j < H.size(); ++j)
    {
        Tensor phij = (useL ? L[j]*phi : phi);
        phij *= H[j];
        if(useR) phij *= R[j];
        Hphi += phij;
    }
    Hphi.mapprime(1,0);
}

template<class Tensor, class OpTensor>
void applyProjOp(const Tensor& phi, const Tensor& L, const Tensor& R, const OpTensor& H, Tensor& Hphi)
{
    bool useL = L.is_not_null(),
         useR = R.is_not_null();
    Hphi = (useL ? L*phi : phi); 
    Hphi *= H;
    if(useR) Hphi *= R;
    Hphi.mapprime(1,0);
}

template<class Tensor>
class BaseLocalHam : public BigMatrix // to do DMRG using an MPO
{
public:
    virtual ~BaseLocalHam() { }
    virtual int Size() const = 0;
    virtual VectorRef DiagRef() const = 0;
    virtual Vector operator*(const VectorRef &A) const = 0;
    void product(const VectorRef &A , VectorRef & B) const = 0;
};

template<class Tensor, class TensorSet>
class LocalHam : public BaseLocalHam<Tensor>
{
    typedef BaseLocalHam<Tensor> Parent;
    Tensor& psi;
    Vector diag;
    const TensorSet &LeftTerm, &RightTerm, &MPOTerm;
public:
    LocalHam(const TensorSet& le, const TensorSet& ri, const TensorSet& mpo, Tensor& psi_) 
	: psi(psi_), LeftTerm(le), RightTerm(ri), MPOTerm(mpo)
    { diag.ReDimension(psi.vec_size()); diag = 1; }

    int Size() const { return psi.vec_size(); }
    VectorRef DiagRef() const { return diag; }

    Vector operator*(const VectorRef &A) const
	{ Vector res(Size()); product(A,res); return res; }

    void product(const VectorRef& A, VectorRef& B) const
	{
        psi.assignFromVec(A);
        Tensor psip; 
        applyProjOp(psi,LeftTerm,RightTerm,MPOTerm,psip);
        psi.assignFrom(psip);
        psi.assignToVec(B);
	}
};

template<class Tensor>
class LocalHamOrth : public BaseLocalHam<Tensor> // to do DMRG using an MPO, ortho to other vecs
{
    typedef BaseLocalHam<Tensor> Parent;

    Tensor& psi;
    Vector diag;
    const Tensor &LeftTerm, &RightTerm, &MPOTerm;
    bool useleft, useright;
    Real weight;
public:
    std::vector<Tensor> other;

    LocalHamOrth(const Tensor& le, const Tensor& ri, const Tensor& mpo, Tensor& psi_, Real weight_) 
	: psi(psi_), LeftTerm(le), RightTerm(ri), MPOTerm(mpo), 
      useleft(le.is_not_null()), useright(ri.is_not_null()), weight(weight_)
    { diag.ReDimension(psi.vec_size()); diag = 1; }

    int Size() const { return psi.vec_size(); }
    VectorRef DiagRef() const { return diag; }

    Vector operator*(const VectorRef &A) const
	{ Vector res(Size()); product(A,res); return res; }

    void product(const VectorRef &A , VectorRef & B) const
	{
        psi.assignFromVec(A);
        Tensor psip;
        applyProjOp(psi,LeftTerm,RightTerm,MPOTerm,psip);
        foreach(const ITensor& phi, other)
        {
            Real re,im; Dot(phi,psi,re,im);
            if(fabs(im) < 1E-10)
            { psip += (weight*re) * phi;}
            else
            { psip += weight*(re*Complex_1 + im*Complex_i) * phi; }
        }
        psi.assignFrom(psip);
        psi.assignToVec(B);
	}
};

template<class Tensor, class TensorSet>
void putInQNs(Tensor& phi, const TensorSet& mpoh, const TensorSet& LH, const TensorSet& RH)
{
    Tensor phip;
    for(int cnt = 1; cnt <= 1E5; ++cnt)
    {
        applyProjOp(phi,LH,RH,mpoh,phip);
        phip *= -0.00232341; //arbitrary small number
        phip += phi; //evolve by (1-tau*H)
        int phisize = phi.vec_size();
        phi = phip;
        if(cnt > 10) std::cerr << "Warning: large number of time evolution steps in putInQNs." << std::endl;
        if(phisize == 0) { if(cnt > 9) Error("phi has zero size in putInQNs."); else continue; }
        else if(phip.vec_size() == phisize) break;
    }
}
template<class Tensor, class TensorSet>
void putInQNs(std::vector<Tensor>& phi, const TensorSet& mpoh, const TensorSet& LH, const TensorSet& RH)
{
    for(size_t n = 0; n < phi.size(); ++n)
    {
        Tensor phip;
        if(phi[n].is_null() || phi[n].vec_size() == 0)
        {
            Print(n); Print(phi[n]);
            Error("Null or zero size tensor in putInQNs.");
        }
        for(int cnt = 1; cnt <= 1E5; ++cnt)
        {
            applyProjOp(phi[n],LH,RH,mpoh,phip);
            phip *= -0.00232341; //arbitrary small number
            phip += phi[n]; //evolve by (1-tau*H)
            int phisize = phi[n].vec_size();
            phi[n] = phip;
            if(cnt > 10) std::cerr << "Warning: large number of time evolution steps in putInQNs." << std::endl;
            if(phisize == 0) { if(cnt > 9) Error("phi has zero size in putInQNs."); else continue; }
            else if(phip.vec_size() == phisize) break;
        }
    }
}
template<class TensorSet>
void putInQNs(ITensor& phi, const TensorSet& mpoh, const TensorSet& LH, const TensorSet& RH) { }

template<class Tensor, class TensorSet>
Real doDavidson(Tensor& phi, const TensorSet& mpoh, const TensorSet& LH, const TensorSet& RH, int niter, int debuglevel, Real errgoal)
{
    putInQNs(phi,mpoh,LH,RH);
    LocalHam<Tensor,TensorSet> lham(LH,RH,mpoh,phi);
    if(niter < 1)
    {
        //Just return the current energy (no optimization via Davidson)
        Vector Phi(phi.vec_size()),HPhi(phi.vec_size()); 
        phi.assignToVec(Phi);
        Phi /= Norm(Phi);
        lham.product(Phi,HPhi);
        Tensor Hphi(phi); Hphi.assignFromVec(HPhi);
        phi.assignFromVec(Phi);
        return Dot(phi,Hphi);
    }
    else
    {
        Matrix evecs(niter,phi.vec_size()); Vector evals;
        phi.assignToVec(evecs.Row(1));
        evecs.Row(1) /= Norm(evecs.Row(1));
        David(lham,1,errgoal,evals,evecs,1,1,debuglevel);
        phi.assignFromVec(evecs.Row(1));
        return evals(1); //energy
    }
    return 1000;
}

template<class Tensor, class TensorSet>
Vector doDavidson(std::vector<Tensor>& phi, const TensorSet& mpoh, const TensorSet& LH, const TensorSet& RH, int niter, int debuglevel, Real errgoal)
{
    const int ntarget = phi.size();
    assert(ntarget != 0);

    putInQNs(phi,mpoh,LH,RH);
    LocalHam<Tensor,TensorSet> lham(LH,RH,mpoh,phi[0]);

    Matrix evecs(max(ntarget,niter),phi[0].vec_size()); Vector evals;
    for(int n = 0; n < ntarget; ++n)
    { 
        phi[n].assignToVec(evecs.Row(1+n)); 
        evecs.Row(1+n) /= Norm(evecs.Row(1+n));
    }
    David(lham,1,errgoal,evals,evecs,1,1,debuglevel);
    Vector energies(ntarget);
    for(int n = 0; n < ntarget; ++n)
    { 
        phi[n].assignFromVec(evecs.Row(1+n));
        energies(1+n) = evals(1+n);
    }
    return energies;
}

void onesite_sweepnext(int &l, int &ha, int N)
{
    if(ha == 1)
    {
        if(++l == N) ha = 2;
        return;
    }
    if(--l == 1) ha = 3;
}


template <class MPSType, class MPOType, class DMRGOptions>
Real onesitedmrg(MPSType& psi, const MPOType& H, const Sweeps& sweeps, DMRGOptions& opts)
{
    typedef typename MPSType::TensorT Tensor;
    typedef typename MPOType::TensorT MPOTensor;
    const Real orig_cutoff = psi.cutoff(); 
    const int orig_minm = psi.minm(), orig_maxm = psi.maxm();
    int debuglevel = (opts.quiet ? 0 : 1);
    int N = psi.NN();
    Real energy;

    psi.position(1);
    //if(H.is_complex()) psi.AAnc(1) *= Complex_1;

    std::vector<MPOTensor> LH(N+1);
    std::vector<MPOTensor> RH(N+1);
    for(int l = N-1; l >= 1; --l) psi.projectOp(l+1,Fromright,RH.at(l+1),H.AA(l+1),RH.at(l));

    for(int sw = 1; sw <= sweeps.nsweep(); ++sw)
    {
    psi.cutoff(sweeps.cutoff(sw)); psi.minm(sweeps.minm(sw)); psi.maxm(sweeps.maxm(sw));
    for(int b = 1, ha = 1; ha != 3; onesite_sweepnext(b,ha,N))
    {
        if(!opts.quiet) 
        {
            std::cout << boost::format("Sweep=%d, HS=%d, Bond=(%d,%d)\n") 
                                % sw   % ha     % b % (b+1);
        }

	    Direction dir = (ha==1?Fromleft:Fromright);
	    Tensor phi = psi.AA(b);

	    const Real errgoal = 1E-4;
	    
	    energy = doDavidson(phi, H.AA(b), LH.at(b), RH.at(b), sweeps.niter(sw), debuglevel, errgoal);

        if(ha == 1)
        {
            phi *= psi.AA(b+1);
            psi.doSVD(b,phi,dir);
        }
        else
        {
            phi *= psi.AA(b-1);
            psi.doSVD(b-1,phi,dir);
        }


        if(!opts.quiet) { std::cout << boost::format("    Truncated to Cutoff=%.1E, Max_m=%d, %s\n") 
                                  % sweeps.cutoff(sw) % sweeps.maxm(sw) 
                                  % (ha == 1 ? psi.LinkInd(b) : psi.LinkInd(b-1)).showm(); }

        opts.measure(sw,ha,(ha==1 ? b : b-1),psi,energy);

        if(ha == 1 && b != N) psi.projectOp(b,Fromleft,LH.at(b),H.AA(b),LH.at(b+1));
        if(ha == 2 && b >= 1)   psi.projectOp(b,Fromright,RH.at(b),H.AA(b),RH.at(b-1));
    } //for loop over b

        if(opts.checkDone(sw,psi,energy))
        {
            psi.cutoff(orig_cutoff); 
            psi.minm(orig_minm); 
            psi.maxm(orig_maxm);
            return energy;
        }

    } //for loop over sw

    psi.cutoff(orig_cutoff); 
    psi.minm(orig_minm); 
    psi.maxm(orig_maxm);
    return energy;
}

template <class MPSType, class MPOType>
Real onesitedmrg(MPSType& psi, const MPOType& H, const Sweeps& sweeps)
{
    DMRGOpts opts; 
    return onesitedmrg(psi,H,sweeps,opts);
}

//Orthogonalizing DMRG. Puts in an energy penalty if psi has an overlap with any MPS in 'other'.
Real dmrg(MPS& psi, const MPO& finalham, const Sweeps& sweeps, 
          const std::vector<MPS>& other, DMRGOpts& opts);

// Deprecated, use MPOSet to work with a set of MPOs
//Real dmrg(MPS& psi, const std::vector<MPO>& H, const Sweeps& sweeps, DMRGOpts& opts);

Real ucdmrg(MPS& psi, const ITensor& LB, const ITensor& RB, const MPO& H, 
            const Sweeps& sweeps, DMRGOpts& opts, bool preserve_edgelink);

#endif
