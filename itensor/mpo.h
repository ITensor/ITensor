//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MPO_H
#define __ITENSOR_MPO_H
#include "mps.h"

#define Cout std::cout
#define Endl std::endl
#define Format boost::format

//
// class MPOt
//
// (defines MPO and IQMPO via typedefs)
//
template<class Tensor>
class MPOt : private MPSt<Tensor>
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

    MPOt(const Model& model, 
         int maxm_ = MAX_M, 
         Real cutoff_ = MIN_CUT, 
         bool _doRelCutoff = true, 
         LogNumber _refNorm = DefaultRefScale);


    //Accessor Methods ------------------------------

    using Parent::N;

    using Parent::model;
    using Parent::isNull;
    using Parent::isNotNull;

    using Parent::si;
    using Parent::siP;

    using Parent::rightLim;
    using Parent::leftLim;

    using Parent::A;
    using Parent::Anc;
    using Parent::bondTensor;

    using Parent::doWrite;
    using Parent::doRelCutoff;
    using Parent::refNorm;
    using Parent::cutoff;
    using Parent::minm;
    using Parent::maxm;
    using Parent::truncerr;
    using Parent::eigsKept;
    using Parent::spectrum;


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
    addAssumeOrth(const MPOt& oth, const OptSet& opts = Global::opts()) 
        { Parent::addAssumeOrth(oth,opts & Opt("UseSVD")); return *this; }

    MPOt& 
    operator+=(const MPOt& oth);


    MPOt 
    operator+(MPOt res) const { res += *this; return res; }

    MPOt 
    operator-(MPOt res) const { res *= -1; res += *this; return res; }

    operator MPOt<IQTensor>()
        { 
        MPOt<IQTensor> res(*model_,maxm(),cutoff(),doRelCutoff(),refNorm()); 
        res.spectrum_ = spectrum_;
        convertToIQ(*model_,A_,res.A_);
        return res; 
        }

    MPOt<ITensor>
    toMPO() const;

    //MPOt: index methods --------------------------------------------------

    using Parent::mapprime;
    using Parent::primelinks;
    using Parent::noprimelink;

    using Parent::LinkInd;
    using Parent::RightLinkInd;
    using Parent::LeftLinkInd;

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
        { Parent::svdBond(b,AA,dir,opts & Opt("UseSVD")); }

    //Move the orthogonality center to site i 
    //(l_orth_lim_ = i-1, r_orth_lim_ = i+1)
    void 
    position(int i, const OptSet& opts = Global::opts()) { Parent::position(i,opts & Opt("UseSVD")); }

    void 
    orthogonalize(const OptSet& opts = Global::opts()) { Parent::orthogonalize(opts & Opt("UseSVD")); }

    using Parent::isOrtho;
    using Parent::orthoCenter;

    using Parent::isComplex;

    using Parent::averageM;

    using Parent::applygate;

    friend inline std::ostream& 
    operator<<(std::ostream& s, const MPOt& M)
        {
        s << "\n";
        for(int i = 1; i <= M.N(); ++i) s << M.A(i) << "\n";
        return s;
        }

    void 
    toIQ(QN totalq, MPOt<IQTensor>& res, Real cut = 1E-12) const
        {
        res = MPOt<IQTensor>(*model_,maxm(),cutoff());
        res.spectrum_ = spectrum_;
        convertToIQ(*model_,A_,res.A_,totalq,cut);
        }

    private:

    ///////////
    using Parent::N_;
    using Parent::A_;
    using Parent::l_orth_lim_;
    using Parent::r_orth_lim_;
    using Parent::model_;
    using Parent::spectrum_;
    ///////////

    friend class MPOt<ITensor>;
    friend class MPOt<IQTensor>;

    }; //class MPOt<Tensor>

typedef MPOt<ITensor> MPO;
typedef MPOt<IQTensor> IQMPO;

template <> inline
MPO MPOt<IQTensor>::
toMPO() const
    {
    MPO res(*model_,maxm(),cutoff(),doRelCutoff(),refNorm());
    res.spectrum_ = spectrum_;
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


int
findCenter(const IQMPO& psi);

void inline 
checkQNs(const MPO& psi) { }

void
checkQNs(const IQMPO& psi);


template <class MPSType, class MPOType>
void 
psiHphi(const MPSType& psi, const MPOType& H, const MPSType& phi, Real& re, Real& im) //<psi|H|phi>
    {
    typedef typename MPSType::TensorT Tensor;
    const int N = H.N();
    if(phi.N() != N || psi.N() != N) Error("psiHphi: mismatched N");

    Tensor L = phi.A(1); 
    //Some Hamiltonians may store edge tensors in H.A(0) and H.A(N+1)
    L *= (H.A(0).isNull() ? H.A(1) : H.A(0)*H.A(1));
    L *= conj(primed(psi.A(1)));
    for(int i = 2; i < N; ++i) 
        { 
        L *= phi.A(i); 
        L *= H.A(i); 
        L *= conj(primed(psi.A(i))); 
        }
    L *= phi.A(N); 
    L *= H.A(N);
    if(!H.A(N+1).isNull()) L *= H.A(N+1);

    Complex z = BraKet(primed(psi.A(N)),L);
    re = z.real();
    im = z.imag();
    }
template <class MPSType, class MPOType>
Real 
psiHphi(const MPSType& psi, const MPOType& H, const MPSType& phi) //Re[<psi|H|phi>]
    {
    Real re, im;
    psiHphi(psi,H,phi,re,im);
    if(fabs(im) > 1.0e-12 * fabs(re))
	std::cerr << boost::format("\nReal psiHphi: WARNING, dropping non-zero (im = %.5f) imaginary part of expectation value.\n")%im;
    return re;
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

    Tensor L = (LB.isNull() ? phi.A(1) : LB * phi.A(1));
    L *= H.A(1); 
    L *= conj(primed(psi.A(1)));
    for(int i = 2; i <= N; ++i)
        { 
        L *= phi.A(i); 
        L *= H.A(i); 
        L *= conj(primed(psi.A(i))); 
        }

    if(!RB.isNull()) L *= RB;

    Complex z = L.toComplex();
    re = z.real();
    im = z.imag();
    }

template <class Tensor>
Real
psiHphi(const MPSt<Tensor>& psi, const MPOt<Tensor>& H, const Tensor& LB, const Tensor& RB, const MPSt<Tensor>& phi) //Re[<psi|H|phi>]
    {
    Real re,im; psiHphi(psi,H,LB,RB,phi,re,im);
    if(fabs(im) > 1.0e-12 * fabs(re))
	std::cout << "Real psiHphi: WARNING, dropping non-zero imaginary part of expectation value.\n";
    return re;
    }

template <class Tensor>
void
psiHKphi(const MPSt<Tensor>& psi, const MPOt<Tensor>& H, const MPOt<Tensor>& K,const MPSt<Tensor>& phi, Real& re, Real& im) //<psi|H K|phi>
    {
    if(psi.N() != phi.N() || psi.N() != H.N() || psi.N() != K.N()) Error("Mismatched N in psiHKphi");
    int N = psi.N();
    MPSt<Tensor> psiconj(psi);
    for(int i = 1; i <= N; i++)
        {
        psiconj.Anc(i) = conj(psi.A(i));
        psiconj.Anc(i).mapprime(0,2);
        }
    MPOt<Tensor> Kp(K);
    Kp.mapprime(1,2);
    Kp.mapprime(0,1);

    //scales as m^2 k^2 d
    Tensor L = (((phi.A(1) * H.A(1)) * Kp.A(1)) * psiconj.A(1));
    for(int i = 2; i < N; i++)
        {
        //scales as m^3 k^2 d + m^2 k^3 d^2
        L = ((((L * phi.A(i)) * H.A(i)) * Kp.A(i)) * psiconj.A(i));
        }
    //scales as m^2 k^2 d
    L = ((((L * phi.A(N)) * H.A(N)) * Kp.A(N)) * psiconj.A(N));
    //cout << "in psiHKpsi, L is "; PrintDat(L);
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
    if(fabs(im) > 1.0e-12 * fabs(re))
	Error("Non-zero imaginary part in psiHKphi");
    return re;
    }

template <class MPOType>
void 
nmultMPO(const MPOType& Aorig, const MPOType& Borig, MPOType& res,Real cut, int maxm);

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
zipUpApplyMPO(const MPSt<Tensor>& psi, const MPOt<Tensor>& K, MPSt<Tensor>& res, Real cutoff = -1, int maxm = -1,
              const OptSet& opts = Global::opts());

//Applies an MPO K to an MPS x with no approximation (|res>=K|x>)
//The bond dimension of res will be the product of bond dimensions
//of x and K.
template<class Tensor>
void 
exactApplyMPO(const MPSt<Tensor>& x, const MPOt<Tensor>& K, MPSt<Tensor>& res);

template<class Tensor>
void
fitApplyMPO(const MPOt<Tensor>& K,
            const MPSt<Tensor>& psi,
            MPSt<Tensor>& res,
            const OptSet& opts = Global::opts());

template<class Tensor>
void
fitApplyMPO(Real fac,
            const MPOt<Tensor>& K,
            const MPSt<Tensor>& psi,
            MPSt<Tensor>& res,
            const OptSet& opts = Global::opts());

template<class Tensor>
Real
fitApplyMPO(const MPSt<Tensor>& psiA, 
            Real mpofac,
            const MPOt<Tensor>& H,
            const MPSt<Tensor>& psiB,
            MPSt<Tensor>& res,
            const OptSet& opts = Global::opts());

template<class Tensor>
Real
fitApplyMPO(Real mpsfac,
            const MPSt<Tensor>& psiA, 
            Real mpofac,
            const MPOt<Tensor>& H,
            const MPSt<Tensor>& psiB,
            MPSt<Tensor>& res,
            const OptSet& opts = Global::opts());

//Computes the exponential of the MPO H: K=exp(-tau*(H-Etot))
template<class Tensor>
void 
expH(const MPOt<Tensor>& H, MPOt<Tensor>& K, Real tau, Real Etot,
     Real Kcutoff, int ndoub);

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
applyExpH(const MPOt<Tensor>& H, 
          Real tau, 
          const MPSt<Tensor>& psi, 
          MPSt<Tensor>& res, 
          const OptSet& opts = Global::opts());

//Given an MPO with no Link indices between site operators,
//put in links (of bond dimension 1).
//In the IQMPO case ensure that links carry the proper QNs.
void
putMPOLinks(MPO& W, const OptSet& opts = Global::opts());
void
putMPOLinks(IQMPO& W, const OptSet& opts = Global::opts());

#undef Cout
#undef Endl
#undef Format

#endif
