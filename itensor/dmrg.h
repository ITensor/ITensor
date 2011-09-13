#ifndef ALGORITHM_H
#define ALGORITHM_H
#include "mpo.h"
#include "sparse.h"
#include "davidson.h"

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
    vector<Tensor> other;

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
        if(cnt > 10) cerr << "Warning: large number of time evolution steps in putInQNs." << endl;
        if(phisize == 0) { if(cnt > 9) Error("phi has zero size in putInQNs."); else continue; }
        else if(phip.vec_size() == phisize) break;
    }
}
template<class Tensor, class TensorSet>
void putInQNs(vector<Tensor>& phi, const TensorSet& mpoh, const TensorSet& LH, const TensorSet& RH)
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
            if(cnt > 10) cerr << "Warning: large number of time evolution steps in putInQNs." << endl;
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
Vector doDavidson(vector<Tensor>& phi, const TensorSet& mpoh, const TensorSet& LH, const TensorSet& RH, int niter, int debuglevel, Real errgoal)
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

enum SweepScheme {ramp_m, fixed_m, fixed_cutoff};

inline void sweepnext(int &l, int &ha, int N, int min_l = 1)
{
    if(ha == 1)
	{
        if(++l == N) 
            l = N-1, ha = 2;
        return;
	}
    if(l-- == min_l) ha = 3;
}

class Sweeps
{
public:
    SweepScheme scheme;
    int Minm;
    vector<int>  Maxm, Niter;
    vector<Real> Cutoff;
    int Nsweep;
    int num_site_center;		// May not be implemented in some cases
    Sweeps(SweepScheme sch, int nsw, int _minm, int _maxm, Real _cut)
	    : scheme(sch), Minm(_minm), Maxm(nsw+1), Niter(nsw+1,4), Cutoff(nsw+1), Nsweep(nsw), num_site_center(2)
	{
        if(scheme == ramp_m)
        {
            for(int s = 1; s <= Nsweep; s++)
            { Cutoff.at(s) = _cut; Maxm.at(s) = (int)(_minm + (s-1.0)/nsw * (_maxm - _minm)); }
        }
        else if(scheme == fixed_m || scheme == fixed_cutoff)
        {
            for(int s = 1; s <= Nsweep; s++)
            { Cutoff.at(s) = _cut; Maxm.at(s) = _maxm; }
        }

        for(int s = 1; s <= min(Nsweep,4); s++)
        { Niter.at(s) = 10 - s; }
	}
    Real cutoff(int sw) const { return Cutoff.at(sw); }
    int minm(int sw) const { return Minm; }
    int maxm(int sw) const { return Maxm.at(sw); }
    int nsweep() const { return Nsweep; }
    int niter(int sw) const { return Niter.at(sw); }
};

//Struct of options for fine-tuning DMRG algorithms
struct DMRGOpts
{
private:
    int largest_m,max_eigs_bond;
    Vector max_eigs,center_eigs;
public:
    Real bulk_entanglement_gap; //Difference between leading denmat evals at center bond
    Real bulk_ent_errgoal;
    Real energy_errgoal;    //Stop DMRG once energy has converged to this precision
    int ntarget;            //Number of states to target
    Real orth_weight;       //How much to penalize non-orthogonality in multiple-state DMRG
    bool printeigs;         //Print slowest decaying eigenvalues after every sweep
    bool quiet;             //Show/don't show info after every step

    DMRGOpts() : 
    largest_m(-1),
    max_eigs_bond(-1),
    max_eigs(1),
    center_eigs(1),
    bulk_entanglement_gap(-1),
    bulk_ent_errgoal(-1),
    energy_errgoal(-1), 
    ntarget(1),
    orth_weight(1),
    printeigs(true), 
    quiet(true)
    { }

    template<class MPSType>
    void measure(int sw, int ha, int b, const MPSType& psi, Real energy)
    {
        if(printeigs)
        {
            Index bi = psi.LinkInd(b);
            if(b == 1 && ha == 1)
            {
                largest_m = -1;
                max_eigs_bond = -1;
                max_eigs = Vector(1); max_eigs = 2;
                center_eigs = Vector(1); center_eigs = 2;
            }

            largest_m = max(largest_m,bi.m());
            assert(lastd.Length() > 0);
            assert(max_eigs.Length() > 0);
            if(lastd(1) < max_eigs(1) && b != 1 && b != (psi.NN()-1)) { max_eigs = lastd; max_eigs_bond = b; }
            if(b == psi.NN()/2) 
            {
                center_eigs = lastd;
                bulk_entanglement_gap = (lastd.Length() >= 2 ? lastd(1)-lastd(2) : 1);
            }

            if(b == 1 && ha == 2) 
            {
                cout << "\n    Largest m during sweep " << sw << " was " << largest_m << "\n";
                cout << boost::format("    Eigs at bond %d: ") % max_eigs_bond;
                for(int j = 1; j <= min(max_eigs.Length(),10); ++j) 
                {
                    cout << boost::format(max_eigs(j) > 1E-2 ? ("%.2f") : ("%.2E")) % max_eigs(j);
                    cout << ((j != min(max_eigs.Length(),10)) ? ", " : "\n");
                }
                cout << "    Eigs at center bond: ";
                for(int j = 1; j <= min(center_eigs.Length(),10); ++j) 
                {
                    cout << boost::format(center_eigs(j) > 1E-2 ? ("%.2f") : ("%.2E")) % center_eigs(j);
                    cout << ((j != min(center_eigs.Length(),10)) ? ", " : "\n");
                }
                cout << boost::format("    Bulk entanglement gap = %f\n") % bulk_entanglement_gap;

                cout << boost::format("    Energy after sweep %d is %f\n") % sw % energy;
            }
        }
    }

    template<class MPSType>
    bool checkDone(int sw, const MPSType& psi, Real energy)
    {
        static Real last_energy;

        if(sw == 1) last_energy = 1000;
        if(energy_errgoal > 0 && sw%2 == 0)
        {
            Real dE = fabs(energy-last_energy);
            if(dE < energy_errgoal)
            {
                cout << boost::format("    Energy error goal met (dE = %E); returning after %d sweeps.\n") % dE % sw;
                return true;
            }
        }
        last_energy = energy;

        return false;
    }
};

template <class Tensor, class DMRGOptions>
Real dmrg(MPSt<Tensor>& psi, const MPOt<Tensor>& H, const Sweeps& sweeps, DMRGOptions& opts)
{
    const Real orig_cutoff = psi.cutoff; const int orig_minm = psi.minm, orig_maxm = psi.maxm;
    int debuglevel = (opts.quiet ? 0 : 1);
    int N = psi.NN();
    Real energy;

    psi.position(1);
    //if(H.is_complex()) psi.AAnc(1) *= Complex_1;

    vector<Tensor> PH(N+1);
    for(int l = N-1; l >= 2; --l) psi.projectOp(l+1,Fromright,PH[l+1],H.AA(l+1),PH[l]);

    for(int sw = 1; sw <= sweeps.nsweep(); ++sw)
    {
        psi.cutoff = sweeps.cutoff(sw); psi.minm = sweeps.minm(sw); psi.maxm = sweeps.maxm(sw);
        for(int b = 1, ha = 1; ha != 3; sweepnext(b,ha,N))
        {
            if(!opts.quiet) cout << boost::format("Sweep=%d, HS=%d, Bond=(%d,%d)\n") % sw % ha % b % (b+1);

            energy = psi.bondDavidson(b,H.bondTensor(b),PH[b],PH[b+1],
                     sweeps.niter(sw),debuglevel,(ha==1?Fromleft:Fromright));

            if(!opts.quiet) { cout << boost::format("    Truncated to Cutoff=%.1E, Max_m=%d, %s\n") 
                                      % sweeps.cutoff(sw) % sweeps.maxm(sw) % psi.LinkInd(b).showm(); }

            opts.measure(sw,ha,b,psi,energy);

            if(ha == 1 && b != N-1) psi.projectOp(b,Fromleft,PH[b],H.AA(b),PH[b+1]);
            if(ha == 2 && b != 1)   psi.projectOp(b+1,Fromright,PH[b+1],H.AA(b+1),PH[b]);
        } //for loop over b

        if(opts.checkDone(sw,psi,energy))
        {
            psi.cutoff = orig_cutoff; psi.minm = orig_minm; psi.maxm = orig_maxm;
            return energy;
        }

    } //for loop over sw

    psi.cutoff = orig_cutoff; psi.minm = orig_minm; psi.maxm = orig_maxm;
    return energy;
}
template <class MPSType, class MPOType>
Real dmrg(MPSType& psi, const MPOType& H, const Sweeps& sweeps)
{
    DMRGOpts opts; 
    return dmrg(psi,H,sweeps,opts);
}

//Orthogonalizing DMRG. Puts in an energy penalty if psi has an overlap with any MPS in 'other'.
Real dmrg(MPS& psi, const MPO& finalham, const Sweeps& sweeps, 
          const vector<MPS>& other, DMRGOpts& opts);

// Deprecated, use MPOSet to work with a set of MPOs
//Real dmrg(MPS& psi, const vector<MPO>& H, const Sweeps& sweeps, DMRGOpts& opts);

Real ucdmrg(MPS& psi, const ITensor& LB, const ITensor& RB, const MPO& H, 
            const Sweeps& sweeps, DMRGOpts& opts, bool preserve_edgelink);

#endif
