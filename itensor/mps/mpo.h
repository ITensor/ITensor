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

class MPO : private MPS
    {
    using Parent = MPS;
    using Parent::N_;
    using Parent::A_;
    using Parent::l_orth_lim_;
    using Parent::r_orth_lim_;
    using Parent::sites_;
    Real logrefNorm_;
    public:

    MPO();

    MPO(int N);

    MPO(SiteSet const& sites, 
         Real _refNorm = DefaultLogRefScale);

    explicit operator bool() const { return Parent::operator bool(); }

    using Parent::length;
    using Parent::sites;

    using Parent::rightLim;
    using Parent::leftLim;

    using Parent::operator();
    using Parent::ref;
    using Parent::set;

    using Parent::doWrite;

    using Parent::read;
    using Parent::write;

    Real
    logRefNorm() const { return logrefNorm_; }
    void
    logRefNorm(Real lrn) { logrefNorm_ = lrn; }

    MPO&
    plusEq(MPO const& R,
           Args const& args = Args::global());

    using Parent::addTags;
    using Parent::removeTags;
    using Parent::replaceTags;
    using Parent::swapTags;
    using Parent::prime;
    using Parent::setPrime;
    using Parent::noPrime;

    void 
    svdBond(int b, 
            ITensor const& AA, 
            Direction dir, 
            Args const& args = Args::global())
        { 
        Parent::svdBond(b,AA,dir,args + Args("UseSVD",true,"LogRefNorm",logrefNorm_)); 
        }

    //Move the orthogonality center to site i 
    //(l_orth_lim_ = i-1, r_orth_lim_ = i+1)
    void 
    position(int i, 
             Args const& args = Args::global()) 
        { 
        Parent::position(i,args + Args("UseSVD")); 
        }

    void 
    orthogonalize(Args const& args = Args::global()) 
        { 
        Parent::orthogonalize(args + Args("UseSVD")); 
        }

    //
    // Deprecations
    //

    using Parent::N;
    using Parent::A;
    using Parent::Aref;
    using Parent::setA;

    }; //class MPO<Tensor>

//template<typename T>
//MPO<T>&
//addAssumeOrth(MPO<T> & L, MPO<T> const& R, Args const& args = Args::global()) 
//    { 
//    MPS<T>::addAssumeOrth(L,R,{args,"UseSVD",true,"LogRefNorm",L.logRefNorm()}); 
//    return L;
//    }


inline MPO& 
operator*=(MPO & W, Real a) { W.ref(W.leftLim()+1) *= a; return W; }

inline MPO& 
operator*=(MPO & W, Cplx a) { W.ref(W.leftLim()+1) *= a; return W; }

inline MPO& 
operator/=(MPO & W, Real a) { W.ref(W.leftLim()+1) /= a; return W; }

inline MPO& 
operator/=(MPO & W, Cplx a) { W.ref(W.leftLim()+1) /= a; return W; }

MPO inline
operator*(MPO W, Real r) { return W *= r; }

MPO inline
operator*(Real r, MPO W) { return W *= r; }

MPO inline
operator*(MPO W, Cplx z) { return W *= z; }

MPO inline
operator*(Cplx z, MPO W) { return W *= z; }

////Convert an IQMPO to an MPO
//MPO
//toMPO(IQMPO const& K);

template<typename... VarArgs>
MPO
addTags(MPO A,
        VarArgs&&... vargs)
    {
    A.addTags(std::forward<VarArgs>(vargs)...);
    return A;
    }
    
template<typename... VarArgs>
MPO
removeTags(MPO A,
           VarArgs&&... vargs)
    {
    A.removeTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
MPO
replaceTags(MPO A,
            VarArgs&&... vargs)
    {
    A.replaceTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
MPO
swapTags(MPO A,
         VarArgs&&... vargs)
    {
    A.swapTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
MPO
prime(MPO A,
      VarArgs&&... vargs)
    {
    A.prime(std::forward<VarArgs>(vargs)...);
    return A; 
    }

template<typename... VarArgs>
MPO
setPrime(MPO A,
         VarArgs&&... vargs)
    {
    A.setPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
MPO
noPrime(MPO A,
        VarArgs&&... vargs)
    {
    A.noPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

Index inline
siteIndex(MPO const& W, int b, TagSet const& tsmatch = TagSet("0"))
    {
    return uniqueIndex(W(b),{W(b-1),W(b+1)},tsmatch);
    }

inline int
length(MPO const& W) { return W.length(); }

bool
isComplex(MPO const& W);

bool
isOrtho(MPO const& W);

int
orthoCenter(MPO const& W);

int
findCenter(MPO const& psi);

void
checkQNs(MPO const& psi);

//MPO
//sum(MPO L, 
//    MPO const& R, 
//    Args const& args = Args::global());

//<psi|H|phi>
void 
overlap(MPS const& psi, 
        MPO const& H, 
        MPS const& phi, 
        Real& re, 
        Real& im);

Real 
overlap(MPS const& psi, 
        MPO const& H, 
        MPS const& phi);

Complex 
overlapC(MPS const& psi, 
         MPO const& H, 
         MPS const& phi);

void
overlap(MPS const& psi, 
        MPO const& H, 
        ITensor const& LB, 
        ITensor const& RB, 
        MPS const& phi, 
        Real& re, 
        Real& im);

Real
overlap(MPS const& psi, 
        MPO const& H, 
        ITensor const& LB, 
        ITensor const& RB, 
        MPS const& phi);

void
overlap(MPS const& psi, 
        MPO const& H, 
        MPO const& K,
        MPS const& phi, 
        Real& re, 
        Real& im);

Real
overlap(MPS const& psi, 
        MPO const& H, 
        MPO const& K,
        MPS const& phi);

Complex
overlapC(MPS const& psi, 
         MPO const& H, 
         MPO const& K,
         MPS const& phi);

void 
nmultMPO(MPO const& Aorig, 
         MPO const& Borig, 
         MPO & res,
         Args args = Args::global());

MPS
applyMPO(MPO const& K,
         MPS const& x,
         Args const& args = Args::global());

MPS
applyMPO(MPO const& K,
         MPS const& x,
         MPS const& x0,
         Args const& args = Args::global());

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
void 
zipUpApplyMPO(MPS const& psi, 
              MPO const& K, 
              MPS & res, 
              Args const& args = Args::global());

//
//Applies an MPO K to an MPS x (K|x>) with no approximation
//made in the application of K to x. Compresses
//the result back into an MPS whose bond dimension
//is at most the product of the bond dimension of K
//and the bond dimension of x. The result can 
//be controllably truncated further by providing
//optional truncation args "Cutoff" and "MaxDim"
//
MPS
exactApplyMPO(MPO const& K,
              MPS const& x,
              Args args = Args::global());



//Applies an MPO K to an MPS psi (|res>=K|psi>) using a sweeping/DMRG-like
//fitting approach. Warning: this method can get stuck i.e. fail to converge
//if the initial value of res is too different from the product K|psi>.
//List of options recognized:
//   Normalize (default: true) - normalize state to 1 after applying MPO
//   Nsweep (default: 1) - number of sweeps to use
//   MaxDim (default: res.maxdim()) - maximum number of states to keep
//   MinDim (default: res.mindim()) - minimum number of states to keep
//   Cutoff (default: res.cutoff()) - maximum truncation error goal
MPS
fitApplyMPO(MPS const& psi,
            MPO const& K,
            Args const& args = Args::global());

void
fitApplyMPO(MPS const& psi,
            MPO const& K,
            MPS & res,
            Args const& args = Args::global());

//Applies an MPO K to an MPS psi including an overall scalar factor (|res>=fac*K|psi>) 
//using a sweeping/DMRG-like fitting approach. 
//Warning: this method can get stuck i.e. fail to converge
//if the initial value of res is too different from the product fac*K|psi>.
//   Normalize (default: true) - normalize state to 1 after applying MPO
//   Nsweep (default: 1) - number of sweeps to use
//   MaxDim (default: res.maxdim()) - maximum number of states to keep
//   MinDim (default: res.mindim()) - minimum number of states to keep
//   Cutoff (default: res.cutoff()) - maximum truncation error goal
void
fitApplyMPO(Real fac,
            MPS const& psi,
            MPO const& K,
            MPS & res,
            Args args = Args::global());

//Applies an MPO K to an MPS psi including an overall scalar factor (|res>=fac*K|psi>) 
//using a sweeping/DMRG-like fitting approach. 
//Warning: this method can get stuck i.e. fail to converge
//if the initial value of res is too different from the product fac*K|psi>.
//Try setting noise > 0 in the Sweeps argument to overcome this.
//Arguments recognized:
//   Verbose (default: false): print out extra information
//   Normalize (default: true): normalize the state to 1 after applying MPO
//
void
fitApplyMPO(Real fac,
            MPS const& psi,
            MPO const& K,
            MPS & res,
            Sweeps const& sweeps,
            Args args);

//Computes |res> = |psiA> + mpofac*H*|psiB>
//using a sweeping/DMRG-like fitting approach. 
//Warning: this method can get stuck i.e. fail to converge
//if the initial value of res is too different from desired exact result.
//   Nsweep (default: 1) - number of sweeps to use
//   MaxDim (default: res.maxdim()) - maximum number of states to keep
//   MinDim (default: res.mindim()) - minimum number of states to keep
//   Cutoff (default: res.cutoff()) - maximum truncation error goal
Real
fitApplyMPO(MPS const& psiA, 
            Real mpofac,
            MPS const& psiB,
            MPO const& H,
            MPS & res,
            Args const& args = Args::global());

//Computes |res> = mpsfac*|psiA> + mpofac*H*|psiB>
//using a sweeping/DMRG-like fitting approach. 
//Warning: this method can get stuck i.e. fail to converge
//if the initial value of res is too different from desired exact result.
//   Nsweep (default: 1) - number of sweeps to use
//   MaxDim (default: res.maxdim()) - maximum number of states to keep
//   MinDim (default: res.mindim()) - minimum number of states to keep
//   Cutoff (default: res.cutoff()) - maximum truncation error goal
Real
fitApplyMPO(Real mpsfac,
            MPS const& psiA, 
            Real mpofac,
            MPS const& psiB,
            MPO const& H,
            MPS & res,
            Args const& args = Args::global());

//Computes the exponential of the MPO H: K=exp(-tau*(H-Etot))
void 
expH(MPO const& H, 
     MPO & K, 
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
//   "MaxDim"  : maximum number of states after truncation
//   "MinDim"  : minimum number of states after truncation
//   "Nsweep": number of sweeps used to apply H MPO to intermediate MPS
//
void
applyExpH(MPS const& psi, 
          MPO const& H, 
          Real tau, 
          MPS & res, 
          Args const& args = Args::global());

//Given an MPO with no Link indices between site operators,
//put in links (of bond dimension 1).
//In the QN conserving case ensure that links carry the proper QNs.
void
putMPOLinks(MPO& W, Args const& args = Args::global());

std::ostream& 
operator<<(std::ostream& s, MPO const& M);

Real
errorMPOProd(MPS const& psi2,
             MPO const& K, 
             MPS const& psi1);

Real
checkMPOProd(MPS const& psi2,
             MPO const& K, 
             MPS const& psi1);

bool
checkMPOProd(MPS const& psi2,
             MPO const& K, 
             MPS const& psi1,
             Real threshold);

//
// Deprecated interfaces - kept for backwards compatibility
//

//Older interface for exactApplyMPO with different ordering
MPS
exactApplyMPO(MPS const& x,
              MPO const& K,
              Args const& args = Args::global());

//Older interface for exactApplyMPO with reference argument for result
void 
exactApplyMPO(MPS const& x, 
              MPO const& K, 
              MPS      & res,
              Args const& args = Args::global());

} //namespace itensor

#endif
