//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MPO_H
#define __ITENSOR_MPO_H
#include "mps.h"


namespace itensor {

static const Real DefaultLogRefScale(2.0255);

//
// class MPOt
//
// (defines MPO and IQMPO via typedefs)
//
template<class Tensor>
class MPOt : private MPSt<Tensor>, public safe_bool<MPOt<Tensor> >
    {
    public:

    typedef MPSt<Tensor> 
    Parent;

    typedef Tensor 
    TensorT;

    typedef typename Tensor::IndexT 
    IndexT;

    typedef typename Tensor::IndexValT 
    IndexValT;

    typedef typename Tensor::CombinerT 
    CombinerT;

    //MPOt: Constructors -----------------------------------------

    MPOt();

    MPOt(const SiteSet& sites, 
         Real _refNorm = DefaultLogRefScale);

    //Accessor Methods ------------------------------

    using Parent::N;

    using Parent::sites;
    using Parent::valid;

    using Parent::rightLim;
    using Parent::leftLim;

    using Parent::A;
    using Parent::Anc;

    using Parent::doWrite;

    using Parent::read;
    using Parent::write;

    //MPOt: operators ------------------------------------------------------

    MPOt& 
    operator*=(Real a) { Parent::operator*=(a); return *this; }
    MPOt
    operator*(Real r) const { MPOt res(*this); res *= r; return res; }
    friend MPOt inline
    operator*(Real r, MPOt res) { res *= r; return res; }

    MPOt& 
    operator*=(Complex z) { Parent::operator*=(z); return *this; }
    MPOt
    operator*(Complex z) const { MPOt res(*this); res *= z; return res; }
    friend MPOt inline
    operator*(Complex z, MPOt res) { res *= z; return res; }

    MPOt&
    plusEq(const MPOt& R,
           const OptSet& opts = Global::opts());


    MPOt<ITensor>
    toMPO() const;

    MPOt<IQTensor>
    toIQMPO() const;

    //MPOt: index methods --------------------------------------------------

    using Parent::mapprime;
    using Parent::primelinks;
    using Parent::noprimelink;

    void 
    primeall()	// sites i,i' -> i',i'';  link:  l -> l'
        {
        for(int i = 1; i <= this->N(); i++)
            {
            Anc(i).mapprime(0,1,Link);
            Anc(i).mapprime(1,2,Site);
            Anc(i).mapprime(0,1,Site);
            }
        }

    void 
    svdBond(int b, const Tensor& AA, Direction dir, const OptSet& opts = Global::opts())
        { 
        Parent::svdBond(b,AA,dir,opts & Opt("UseSVD") & Opt("LogRefNorm",logrefNorm_)); 
        }

    //Move the orthogonality center to site i 
    //(l_orth_lim_ = i-1, r_orth_lim_ = i+1)
    void 
    position(int i, const OptSet& opts = Global::opts()) { Parent::position(i,opts & Opt("UseSVD")); }

    void 
    orthogonalize(const OptSet& opts = Global::opts()) { Parent::orthogonalize(opts & Opt("UseSVD")); }

    using Parent::isOrtho;
    using Parent::orthoCenter;

    using Parent::isComplex;

    void 
    toIQ(QN totalq, MPOt<IQTensor>& res, Real cut = 1E-12) const
        {
        res = MPOt<IQTensor>(*sites_,logrefNorm_);
        convertToIQ(*sites_,A_,res.A_,totalq,cut);
        }

    private:

    ///////////
    using Parent::N_;
    using Parent::A_;
    using Parent::l_orth_lim_;
    using Parent::r_orth_lim_;
    using Parent::sites_;
    Real logrefNorm_;
    ///////////

    MPOt&
    addAssumeOrth(const MPOt& oth, const OptSet& opts = Global::opts()) 
        { 
        Parent::addAssumeOrth(oth,opts & Opt("UseSVD") & Opt("LogRefNorm",logrefNorm_)); 
        return *this; 
        }

    friend class MPOt<ITensor>;
    friend class MPOt<IQTensor>;

    }; //class MPOt<Tensor>

typedef MPOt<ITensor> MPO;
typedef MPOt<IQTensor> IQMPO;

template <> inline
MPO IQMPO::
toMPO() const
    {
    MPO res(*sites_,logrefNorm_);
    for(int j = 0; j <= N()+1; ++j)
        {
        res.A_.at(j) = A(j).toITensor();
        }
    return res;
    }

//toMPO method fails unless template class 
//Tensor is set to IQTensor (object is an IQMPO)
template<class Tensor>
MPO MPOt<Tensor>::
toMPO() const
    {
    Error("toMPO only implemented for class IQMPO");
    return MPO();
    }

template <> inline
IQMPO MPO::
toIQMPO() const
    {
    IQMPO res(*sites_,logrefNorm_);
    convertToIQ(*sites_,A_,res.A_);
    return res;
    }

//toMPO method fails unless template class 
//Tensor is set to IQTensor (object is an IQMPO)
template<class Tensor>
IQMPO MPOt<Tensor>::
toIQMPO() const
    {
    Error("toIQMPO only implemented for class MPO");
    return IQMPO();
    }


int
findCenter(const IQMPO& psi);

void inline 
checkQNs(const MPO& psi) { }

void
checkQNs(const IQMPO& psi);

template <class Tensor>
MPOt<Tensor>
sum(const MPOt<Tensor>& L, 
    const MPOt<Tensor>& R, 
    const OptSet& opts = Global::opts())
    {
    MPOt<Tensor> res(L);
    res.plusEq(R,opts);
    return res;
    }


template <class Tensor>
void 
psiHphi(const MPSt<Tensor>& psi, const MPOt<Tensor>& H, const MPSt<Tensor>& phi, Real& re, Real& im) //<psi|H|phi>
    {
    const int N = H.N();
    if(phi.N() != N || psi.N() != N) Error("psiHphi: mismatched N");

    Tensor L = phi.A(1); 
    //Some Hamiltonians may store edge tensors in H.A(0) and H.A(N+1)
    L *= (H.A(0) ? H.A(0)*H.A(1) : H.A(1));
    L *= dag(prime(psi.A(1)));
    for(int i = 2; i < N; ++i) 
        { 
        L *= phi.A(i); 
        L *= H.A(i); 
        L *= dag(prime(psi.A(i))); 
        }
    L *= phi.A(N); 
    L *= H.A(N);
    if(H.A(N+1)) L *= H.A(N+1);

    Complex z = BraKet(prime(psi.A(N)),L);
    re = z.real();
    im = z.imag();
    }
template <class Tensor>
Real 
psiHphi(const MPSt<Tensor>& psi, const MPOt<Tensor>& H, const MPSt<Tensor>& phi) //Re[<psi|H|phi>]
    {
    Real re, im;
    psiHphi(psi,H,phi,re,im);
    if(std::fabs(im) > 1.0e-12 * std::fabs(re))
        printfln("\nReal psiHphi: WARNING, dropping non-zero (=%.5E) imaginary part of expectation value.",im);
    return re;
    }
template <class Tensor>
Complex 
psiHphiC(const MPSt<Tensor>& psi, const MPOt<Tensor>& H, const MPSt<Tensor>& phi) //Re[<psi|H|phi>]
    {
    Real re, im;
    psiHphi(psi,H,phi,re,im);
    return Complex(re,im);
    }

template<class Tensor>
void
psiHphi(const MPSt<Tensor>& psi, 
        const MPOt<Tensor>& H, 
        const Tensor& LB, 
        const Tensor& RB, 
        const MPSt<Tensor>& phi, 
        Real& re, 
        Real& im) //<psi|H|phi>
    {
    int N = psi.N();
    if(N != phi.N() || H.N() < N) Error("mismatched N in psiHphi");

    Tensor L = (LB ? LB*phi.A(1) : phi.A(1));
    L *= H.A(1); 
    L *= dag(prime(psi.A(1)));
    for(int i = 2; i <= N; ++i)
        { 
        L *= phi.A(i); 
        L *= H.A(i); 
        L *= dag(prime(psi.A(i))); 
        }

    if(RB) L *= RB;

    Complex z = L.toComplex();
    re = z.real();
    im = z.imag();
    }

template <class Tensor>
Real
psiHphi(const MPSt<Tensor>& psi, const MPOt<Tensor>& H, const Tensor& LB, const Tensor& RB, const MPSt<Tensor>& phi) //Re[<psi|H|phi>]
    {
    Real re,im; psiHphi(psi,H,LB,RB,phi,re,im);
    if(std::fabs(im) > 1.0e-12 * std::fabs(re))
        printfln("Real psiHphi: WARNING, dropping non-zero imaginary part (=%.5E) of expectation value.",im);
    return re;
    }

template <class Tensor>
void
psiHKphi(const MPSt<Tensor>& psi, const MPOt<Tensor>& H, const MPOt<Tensor>& K,const MPSt<Tensor>& phi, Real& re, Real& im) //<psi|H K|phi>
    {
    if(psi.N() != phi.N() || psi.N() != H.N() || psi.N() != K.N()) Error("Mismatched N in psiHKphi");
    int N = psi.N();
    MPSt<Tensor> psidag(psi);
    for(int i = 1; i <= N; i++)
        {
        psidag.Anc(i) = dag(psi.A(i));
        psidag.Anc(i).mapprime(0,2);
        }
    MPOt<Tensor> Kp(K);
    Kp.mapprime(1,2);
    Kp.mapprime(0,1);

    //scales as m^2 k^2 d
    Tensor L = (((phi.A(1) * H.A(1)) * Kp.A(1)) * psidag.A(1));
    for(int i = 2; i < N; i++)
        {
        //scales as m^3 k^2 d + m^2 k^3 d^2
        L = ((((L * phi.A(i)) * H.A(i)) * Kp.A(i)) * psidag.A(i));
        }
    //scales as m^2 k^2 d
    L = ((((L * phi.A(N)) * H.A(N)) * Kp.A(N)) * psidag.A(N));
    //cout << "in psiHKpsi, L is "; PrintData(L);
    Complex z = L.toComplex();
    re = z.real();
    im = z.imag();
    }

template <class Tensor>
Real
psiHKphi(const MPSt<Tensor>& psi, const MPOt<Tensor>& H, const MPOt<Tensor>& K,const MPSt<Tensor>& phi) //<psi|H K|phi>
    {
    Real re,im;
    psiHKphi(psi,H,K,phi,re,im);
    if(std::fabs(im) > 1.0e-12 * std::fabs(re))
	Error("Non-zero imaginary part in psiHKphi");
    return re;
    }

template <class Tensor>
Complex
psiHKphiC(const MPSt<Tensor>& psi, const MPOt<Tensor>& H, const MPOt<Tensor>& K,const MPSt<Tensor>& phi) //<psi|H K|phi>
    {
    Real re,im;
    psiHKphi(psi,H,K,phi,re,im);
    return Complex(re,im);
    }

template <class MPOType>
void 
nmultMPO(const MPOType& Aorig, const MPOType& Borig, MPOType& res,
         const OptSet& opts = Global::opts());

//
// Applies an MPO to an MPS using the zip-up method described
// more fully in Stoudenmire and White, New. J. Phys. 12, 055026 (2010).
//
// This method applies the MPO to an MPS one site at a time,
// with the new MPS being calculated at each step via an
// SVD of the MPO-MPS product.
//
// Uses cutoff and max of MPS psi unless specified.
//
template<class Tensor>
void 
zipUpApplyMPO(const MPSt<Tensor>& psi, 
              const MPOt<Tensor>& K, 
              MPSt<Tensor>& res, 
              const OptSet& opts = Global::opts());

//Applies an MPO K to an MPS x with no approximation (|res>=K|x>)
//The bond dimension of res will be the product of bond dimensions
//of x and K.
template<class Tensor>
void 
exactApplyMPO(const MPSt<Tensor>& x, 
              const MPOt<Tensor>& K, 
              MPSt<Tensor>& res,
              const OptSet& opts = Global::opts());

//Applies an MPO K to an MPS psi (|res>=K|psi>) using a sweeping/DMRG-like
//fitting approach. Warning: this method can get stuck i.e. fail to converge
//if the initial value of res is too different from the product K|psi>.
//List of options recognized:
//   Nsweep (default: 1) - number of sweeps to use
//   Maxm (default: res.maxm()) - maximum number of states to keep
//   Minm (default: res.minm()) - minimum number of states to keep
//   Cutoff (default: res.cutoff()) - maximum truncation error goal
template<class Tensor>
void
fitApplyMPO(const MPSt<Tensor>& psi,
            const MPOt<Tensor>& K,
            MPSt<Tensor>& res,
            const OptSet& opts = Global::opts());

//Applies an MPO K to an MPS psi including an overall scalar factor (|res>=fac*K|psi>) 
//using a sweeping/DMRG-like fitting approach. 
//Warning: this method can get stuck i.e. fail to converge
//if the initial value of res is too different from the product fac*K|psi>.
//   Nsweep (default: 1) - number of sweeps to use
//   Maxm (default: res.maxm()) - maximum number of states to keep
//   Minm (default: res.minm()) - minimum number of states to keep
//   Cutoff (default: res.cutoff()) - maximum truncation error goal
template<class Tensor>
void
fitApplyMPO(Real fac,
            const MPSt<Tensor>& psi,
            const MPOt<Tensor>& K,
            MPSt<Tensor>& res,
            const OptSet& opts = Global::opts());

//Computes |res> = |psiA> + mpofac*H*|psiB>
//using a sweeping/DMRG-like fitting approach. 
//Warning: this method can get stuck i.e. fail to converge
//if the initial value of res is too different from desired exact result.
//   Nsweep (default: 1) - number of sweeps to use
//   Maxm (default: res.maxm()) - maximum number of states to keep
//   Minm (default: res.minm()) - minimum number of states to keep
//   Cutoff (default: res.cutoff()) - maximum truncation error goal
template<class Tensor>
Real
fitApplyMPO(const MPSt<Tensor>& psiA, 
            Real mpofac,
            const MPSt<Tensor>& psiB,
            const MPOt<Tensor>& H,
            MPSt<Tensor>& res,
            const OptSet& opts = Global::opts());

//Computes |res> = mpsfac*|psiA> + mpofac*H*|psiB>
//using a sweeping/DMRG-like fitting approach. 
//Warning: this method can get stuck i.e. fail to converge
//if the initial value of res is too different from desired exact result.
//   Nsweep (default: 1) - number of sweeps to use
//   Maxm (default: res.maxm()) - maximum number of states to keep
//   Minm (default: res.minm()) - minimum number of states to keep
//   Cutoff (default: res.cutoff()) - maximum truncation error goal
template<class Tensor>
Real
fitApplyMPO(Real mpsfac,
            const MPSt<Tensor>& psiA, 
            Real mpofac,
            const MPSt<Tensor>& psiB,
            const MPOt<Tensor>& H,
            MPSt<Tensor>& res,
            const OptSet& opts = Global::opts());

//Computes the exponential of the MPO H: K=exp(-tau*(H-Etot))
template<class Tensor>
void 
expH(const MPOt<Tensor>& H, MPOt<Tensor>& K, Real tau, Real Etot,
     Real Kcutoff, int ndoub, const OptSet& opts = Global::opts());

//
//Approximately computes |res> = exp(-tau*H)|psi>.
//Works by expanding the exponential to order n.
//The default order is n=4 but can be increased via the "Order" Opt.
//E.g. applyExpH(H,tau,psi,res,opts&Opt("Order",n));
//List of options recognized:
//   Order  - order of Taylor series expansion of exp(-tau*H)
//   Cutoff - maximum truncation error allowed
//   Maxm   - maximum number of states after truncation
//   Minm   - minimum number of states after truncation
//   Nsweep - number of sweeps used to apply H MPO to intermediate MPS
//
template<class Tensor>
void
applyExpH(const MPSt<Tensor>& psi, 
          const MPOt<Tensor>& H, 
          Real tau, 
          MPSt<Tensor>& res, 
          const OptSet& opts = Global::opts());

//Given an MPO with no Link indices between site operators,
//put in links (of bond dimension 1).
//In the IQMPO case ensure that links carry the proper QNs.
void
putMPOLinks(MPO& W, const OptSet& opts = Global::opts());
void
putMPOLinks(IQMPO& W, const OptSet& opts = Global::opts());

template <class Tensor>
std::ostream& 
operator<<(std::ostream& s, const MPOt<Tensor>& M);

}; //namespace itensor


#endif
