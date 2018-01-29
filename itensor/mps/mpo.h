//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MPO_H
#define __ITENSOR_MPO_H
#include "itensor/mps/mps.h"
#include "itensor/mps/sweeps.h"


namespace itensor {

const Real DefaultLogRefScale = 2.0255;

template<class Tensor> 
class MPOt;

using MPO = MPOt<ITensor>;
using IQMPO = MPOt<IQTensor>;

//
// class MPOt
//
// (defines MPO and IQMPO via above typedefs)
//
template<class Tensor>
class MPOt : private MPSt<Tensor>
    {
    using Parent = MPSt<Tensor>;
    using Parent::N_;
    using Parent::A_;
    using Parent::l_orth_lim_;
    using Parent::r_orth_lim_;
    using Parent::sites_;
    Real logrefNorm_;
    public:
    using TensorT = Tensor;
    using IndexT = typename Tensor::index_type;
    using IndexValT = typename Tensor::indexval_type;

    MPOt();

    MPOt(int N);

    MPOt(SiteSet const& sites, 
         Real _refNorm = DefaultLogRefScale);

    explicit operator bool() const { return Parent::operator bool(); }

    using Parent::N;
    using Parent::sites;

    using Parent::rightLim;
    using Parent::leftLim;

    using Parent::A;
    using Parent::Aref;
    using Parent::Anc;
    using Parent::setA;

    using Parent::doWrite;

    using Parent::read;
    using Parent::write;

    Real
    logRefNorm() const { return logrefNorm_; }
    void
    logRefNorm(Real lrn) { logrefNorm_ = lrn; }

    MPOt&
    plusEq(const MPOt& R,
           const Args& args = Args::global());


    MPOt<ITensor>
    toMPO() const;

    //MPOt<IQTensor>
    //toIQMPO() const;

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
    svdBond(int b, const Tensor& AA, Direction dir, const Args& args = Args::global())
        { 
        Parent::svdBond(b,AA,dir,args + Args("UseSVD",true,"LogRefNorm",logrefNorm_)); 
        }

    //Move the orthogonality center to site i 
    //(l_orth_lim_ = i-1, r_orth_lim_ = i+1)
    void 
    position(int i, const Args& args = Args::global()) { Parent::position(i,args + Args("UseSVD")); }

    void 
    orthogonalize(Args const& args = Args::global()) 
        { 
        Parent::orthogonalize(args + Args("UseSVD")); 
        }


    private:


    friend class MPOt<ITensor>;
    friend class MPOt<IQTensor>;
    
    public:

    //
    // Deprecated methods
    //
    //use isOrtho(W) instead
    using Parent::isOrtho;
    //use orthoCenter(W) instead
    using Parent::orthoCenter;
    //use isComplex(W) instead
    using Parent::isComplex;

    //void 
    //toIQ(QN totalq, MPOt<IQTensor>& res, Real cut = 1E-12) const
    //    {
    //    res = MPOt<IQTensor>(*sites_,logrefNorm_);
    //    convertToIQ(*sites_,A_,res.A_,totalq,cut);
    //    }

    }; //class MPOt<Tensor>

//template<typename T>
//MPOt<T>&
//addAssumeOrth(MPOt<T> & L, MPOt<T> const& R, Args const& args = Args::global()) 
//    { 
//    MPSt<T>::addAssumeOrth(L,R,{args,"UseSVD",true,"LogRefNorm",L.logRefNorm()}); 
//    return L;
//    }


template<class T>
MPOt<T>& 
operator*=(MPOt<T> & W, Real a) { W.Anc(W.leftLim()+1) *= a; return W; }

template<class T>
MPOt<T>& 
operator*=(MPOt<T> & W, Cplx a) { W.Anc(W.leftLim()+1) *= a; return W; }

template<class T>
MPOt<T>& 
operator/=(MPOt<T> & W, Real a) { W.Anc(W.leftLim()+1) /= a; return W; }

template<class T>
MPOt<T>& 
operator/=(MPOt<T> & W, Cplx a) { W.Anc(W.leftLim()+1) /= a; return W; }

template<typename T>
MPOt<T>
operator*(MPOt<T> W, Real r) { return W *= r; }

template<typename T>
MPOt<T>
operator*(Real r, MPOt<T> W) { return W *= r; }

template<typename T>
MPOt<T>
operator*(MPOt<T> W, Cplx z) { return W *= z; }

template<typename T>
MPOt<T>
operator*(Cplx z, MPOt<T> W) { return W *= z; }

//Convert an IQMPO to an MPO
MPO
toMPO(IQMPO const& K);

template<typename T>
bool
isComplex(MPOt<T> const& W);

template<typename T>
bool
isOrtho(MPOt<T> const& W);

template<typename T>
int
orthoCenter(MPOt<T> const& W);

int
findCenter(IQMPO const& psi);

void inline 
checkQNs(MPO const& psi) { }

void
checkQNs(IQMPO const& psi);

template <class Tensor>
MPOt<Tensor>
sum(MPOt<Tensor> L, 
    MPOt<Tensor> const& R, 
    Args const& args = Args::global());

//<psi|H|phi>
template <class Tensor>
void 
overlap(MPSt<Tensor> const& psi, 
        MPOt<Tensor> const& H, 
        MPSt<Tensor> const& phi, 
        Real& re, 
        Real& im);

template <class Tensor>
Real 
overlap(MPSt<Tensor> const& psi, 
        MPOt<Tensor> const& H, 
        MPSt<Tensor> const& phi);

template <class Tensor>
Complex 
overlapC(MPSt<Tensor> const& psi, 
         MPOt<Tensor> const& H, 
         MPSt<Tensor> const& phi);

template<class Tensor>
void
overlap(MPSt<Tensor> const& psi, 
        MPOt<Tensor> const& H, 
        Tensor const& LB, 
        Tensor const& RB, 
        MPSt<Tensor> const& phi, 
        Real& re, 
        Real& im);

template <class Tensor>
Real
overlap(MPSt<Tensor> const& psi, 
        MPOt<Tensor> const& H, 
        Tensor const& LB, 
        Tensor const& RB, 
        MPSt<Tensor> const& phi);

template <class Tensor>
void
overlap(MPSt<Tensor> const& psi, 
         MPOt<Tensor> const& H, 
         MPOt<Tensor> const& K,
         MPSt<Tensor> const& phi, 
         Real& re, 
         Real& im);

template <class Tensor>
Real
overlap(MPSt<Tensor> const& psi, 
        MPOt<Tensor> const& H, 
        MPOt<Tensor> const& K,
        MPSt<Tensor> const& phi);

template <class Tensor>
Complex
overlapC(MPSt<Tensor> const& psi, 
         MPOt<Tensor> const& H, 
         MPOt<Tensor> const& K,
         MPSt<Tensor> const& phi);

template<class MPOType>
void 
nmultMPO(MPOType const& Aorig, 
         MPOType const& Borig, 
         MPOType& res,
         Args args = Args::global());

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
zipUpApplyMPO(MPSt<Tensor> const& psi, 
              MPOt<Tensor> const& K, 
              MPSt<Tensor>& res, 
              Args const& args = Args::global());

//
//Applies an MPO K to an MPS x (K|x>) with no approximation
//made in the application of K to x. Compresses
//the result back into an MPS whose bond dimension
//is at most the product of the bond dimension of K
//and the bond dimension of x. The result can 
//be controllably truncated further by providing
//optional truncation args "Cutoff" and "Maxm"
//
template<class Tensor>
MPSt<Tensor>
exactApplyMPO(MPOt<Tensor> const& K,
              MPSt<Tensor> const& x,
              Args const& args = Args::global());



//Applies an MPO K to an MPS psi (|res>=K|psi>) using a sweeping/DMRG-like
//fitting approach. Warning: this method can get stuck i.e. fail to converge
//if the initial value of res is too different from the product K|psi>.
//List of options recognized:
//   Normalize (default: true) - normalize state to 1 after applying MPO
//   Nsweep (default: 1) - number of sweeps to use
//   Maxm (default: res.maxm()) - maximum number of states to keep
//   Minm (default: res.minm()) - minimum number of states to keep
//   Cutoff (default: res.cutoff()) - maximum truncation error goal
template<class Tensor>
MPSt<Tensor>
fitApplyMPO(MPSt<Tensor> const& psi,
            MPOt<Tensor> const& K,
            Args const& args = Args::global());

template<class Tensor>
void
fitApplyMPO(MPSt<Tensor> const& psi,
            MPOt<Tensor> const& K,
            MPSt<Tensor>& res,
            Args const& args = Args::global());

//Applies an MPO K to an MPS psi including an overall scalar factor (|res>=fac*K|psi>) 
//using a sweeping/DMRG-like fitting approach. 
//Warning: this method can get stuck i.e. fail to converge
//if the initial value of res is too different from the product fac*K|psi>.
//   Normalize (default: true) - normalize state to 1 after applying MPO
//   Nsweep (default: 1) - number of sweeps to use
//   Maxm (default: res.maxm()) - maximum number of states to keep
//   Minm (default: res.minm()) - minimum number of states to keep
//   Cutoff (default: res.cutoff()) - maximum truncation error goal
template<class Tensor>
void
fitApplyMPO(Real fac,
            MPSt<Tensor> const& psi,
            MPOt<Tensor> const& K,
            MPSt<Tensor>& res,
            Args const& args = Args::global());

//Applies an MPO K to an MPS psi including an overall scalar factor (|res>=fac*K|psi>) 
//using a sweeping/DMRG-like fitting approach. 
//Warning: this method can get stuck i.e. fail to converge
//if the initial value of res is too different from the product fac*K|psi>.
//Try setting noise > 0 in the Sweeps argument to overcome this.
//Arguments recognized:
//   Verbose (default: false): print out extra information
//   Normalize (default: true): normalize the state to 1 after applying MPO
//
template<class Tensor>
void
fitApplyMPO(Real fac,
            MPSt<Tensor> const& psi,
            MPOt<Tensor> const& K,
            MPSt<Tensor>& res,
            Sweeps const& sweeps,
            Args args);

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
fitApplyMPO(MPSt<Tensor> const& psiA, 
            Real mpofac,
            MPSt<Tensor> const& psiB,
            MPOt<Tensor> const& H,
            MPSt<Tensor>& res,
            Args const& args = Args::global());

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
            MPSt<Tensor> const& psiA, 
            Real mpofac,
            MPSt<Tensor> const& psiB,
            MPOt<Tensor> const& H,
            MPSt<Tensor>& res,
            Args const& args = Args::global());

//Computes the exponential of the MPO H: K=exp(-tau*(H-Etot))
template<class Tensor>
void 
expH(MPOt<Tensor> const& H, 
     MPOt<Tensor>& K, 
     Real tau, 
     Real Etot,
     Real Kcutoff, 
     int ndoub, 
     Args args = Args::global());

//
//Approximately computes |res> = exp(-tau*H)|psi>.
//Works by expanding the exponential to order n.
//The default order is n=4 but can be increased via the "Order" argument.
//E.g. applyExpH(H,tau,psi,res,{"Order",n});
//List of named arguments recognized:
//   "Order" : order of Taylor series expansion of exp(-tau*H)
//   "Cutoff": maximum truncation error allowed
//   "Maxm"  : maximum number of states after truncation
//   "Minm"  : minimum number of states after truncation
//   "Nsweep": number of sweeps used to apply H MPO to intermediate MPS
//
template<class Tensor>
void
applyExpH(MPSt<Tensor> const& psi, 
          MPOt<Tensor> const& H, 
          Real tau, 
          MPSt<Tensor>& res, 
          Args const& args = Args::global());

//Given an MPO with no Link indices between site operators,
//put in links (of bond dimension 1).
//In the IQMPO case ensure that links carry the proper QNs.
void
putMPOLinks(MPO& W, Args const& args = Args::global());
void
putMPOLinks(IQMPO& W, Args const& args = Args::global());

template <class Tensor>
std::ostream& 
operator<<(std::ostream& s, MPOt<Tensor> const& M);

template<class Tensor>
Real
checkMPOProd(MPSt<Tensor> const& psi2,
             MPOt<Tensor> const& K, 
             MPSt<Tensor> const& psi1);

//
// Deprecated interfaces - kept for backwards compatibility
//

//Older interface for exactApplyMPO with different ordering
template<class Tensor>
MPSt<Tensor>
exactApplyMPO(MPSt<Tensor> const& x,
              MPOt<Tensor> const& K,
              Args const& args = Args::global());

//Older interface for exactApplyMPO with reference argument for result
template<class Tensor>
void 
exactApplyMPO(MPSt<Tensor> const& x, 
              MPOt<Tensor> const& K, 
              MPSt<Tensor>      & res,
              Args const& args = Args::global());

} //namespace itensor

#include "mpo.ih"

#endif
