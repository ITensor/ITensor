//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef __ITENSOR_MPO_H
#define __ITENSOR_MPO_H
#include "itensor/mps/mps.h"
#include "itensor/mps/sweeps.h"


namespace itensor {

const Real DefaultLogRefScale = 2.0255;

class MPO : protected MPS
    {
    using Parent = MPS;
    using Parent::N_;
    using Parent::A_;
    using Parent::l_orth_lim_;
    using Parent::r_orth_lim_;
    Real logrefNorm_;
    public:

    MPO();

    MPO(int N);

    MPO(SiteSet const& sites, 
         Real _refNorm = DefaultLogRefScale);

    explicit operator bool() const { return Parent::operator bool(); }

    using Parent::length;

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

    using Parent::dag;
    using Parent::replaceLinkInds;
    using Parent::setTags;
    using Parent::noTags;
    using Parent::addTags;
    using Parent::removeTags;
    using Parent::replaceTags;
    using Parent::swapTags;
    using Parent::prime;
    using Parent::setPrime;
    using Parent::mapPrime;
    using Parent::swapPrime;
    using Parent::noPrime;

    // Replace the site indices of MPO A from sites_old
    // to sites_new (replaces one site index per MPO tensor)
    MPO&
    replaceSiteInds(IndexSet const& sites_old, IndexSet const& sites_new);

    // Swap the bra and ket indices of the MPO (i.e. transpose the MPO)
    MPO&
    swapSiteInds();

    Spectrum 
    svdBond(int b, 
            ITensor const& AA, 
            Direction dir, 
            Args const& args = Args::global())
        { 
        return Parent::svdBond(b,AA,dir,args + Args("UseSVD",true,"LogRefNorm",logrefNorm_)); 
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
    // Developer level methods
    //

    using Parent::uref;

    //
    // Deprecations
    //

    using Parent::N;
    using Parent::A;
    using Parent::Aref;
    using Parent::setA;

    }; //class MPO<Tensor>

MPO& 
operator*=(MPO & W, Real a);

MPO& 
operator*=(MPO & W, Cplx a);

MPO& 
operator/=(MPO & W, Real a);

MPO& 
operator/=(MPO & W, Cplx a);

MPO
operator*(MPO W, Real r);

MPO
operator*(Real r, MPO W);

MPO
operator*(MPO W, Cplx z);

MPO
operator*(Cplx z, MPO W);

int
length(MPO const& W);

MPO
dag(MPO W);

//
// MPO Index functions
//

// Check if the MPO A has the site indices
// sites
bool
hasSiteInds(MPO const& A, IndexSet const& sites);

//
// Find site indices
//

// Find a site index of the MPO by matching the tag 'tsmatch'
Index
siteIndex(MPO const& A, int b, TagSet const& tsmatch = TagSet("0"));

// Get the site Indices of the MPO A*B at site b
// as if MPO A and MPO B were contracted
IndexSet
siteInds(MPO const& A, MPO const& B, int b);

// Get the site Index of the MPS A|x>
// as if MPO A was applied to MPS x
Index
uniqueSiteIndex(MPO const& A, MPS const& x, int b);

// Find the site Index of the bth MPO tensor of W
// that is not the input site Index s
Index
uniqueSiteIndex(MPO const& W, IndexSet const& s, int b);

// Get the site Indices that are unique to A
IndexSet
uniqueSiteInds(MPO const& A, MPS const& x);

// Get the site Indices that are unique to A
// (not in the IndexSet of site indices sites)
IndexSet
uniqueSiteInds(MPO const& A, IndexSet const& sites);

// Get the site Index that is unique to A
Index
uniqueSiteIndex(MPO const& A, MPO const& B, int b);

// Get the site Indices that are unique to the MPO A
// If both indices are shared by A and B, they will 
// be default valued indices (Index())
IndexSet
uniqueSiteInds(MPO const& A, MPO const& B);

//
// Modify site indices
//

// Replace the site indices of MPO A from sites_old
// to sites_new (replaces one site index per MPO tensor)
MPO
replaceSiteInds(MPO A, IndexSet const& sites_old, IndexSet const& sites_new);

// Swap the bra and ket indices of the MPO (i.e. transpose the MPO)
MPO
swapSiteInds(MPO A);

//
// MPO tag functions
//

MPO
setTags(MPO A, TagSet const& ts, IndexSet const& is);

template <typename... VarArgs>
MPO
setTags(MPO A,
        VarArgs&&... vargs)
    {
    A.setTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

MPO
noTags(MPO A, IndexSet const& is);

template <typename... VarArgs>
MPO
noTags(MPO A,
       VarArgs&&... vargs)
    {
    A.noTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

MPO
addTags(MPO A, TagSet const& ts, IndexSet const& is);

template <typename... VarArgs>
MPO
addTags(MPO A,
        VarArgs&&... vargs)
    {
    A.addTags(std::forward<VarArgs>(vargs)...);
    return A;
    }
   
MPO
removeTags(MPO A, TagSet const& ts, IndexSet const& is);

template <typename... VarArgs>
MPO
removeTags(MPO A,
           VarArgs&&... vargs)
    {
    A.removeTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

MPO
replaceTags(MPO A, TagSet const& ts1, TagSet const& ts2, IndexSet const& is);

template <typename... VarArgs>
MPO
replaceTags(MPO A,
            VarArgs&&... vargs)
    {
    A.replaceTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

MPO
swapTags(MPO A, TagSet const& ts1, TagSet const& ts2, IndexSet const& is);

template <typename... VarArgs>
MPO
swapTags(MPO A,
         VarArgs&&... vargs)
    {
    A.swapTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

MPO
prime(MPO A, int plev, IndexSet const& is);

MPO
prime(MPO A, IndexSet const& is);

template <typename... VarArgs>
MPO
prime(MPO A,
      VarArgs&&... vargs)
    {
    A.prime(std::forward<VarArgs>(vargs)...);
    return A;
    }

MPO
setPrime(MPO A, int plev, IndexSet const& is);

template <typename... VarArgs>
MPO
setPrime(MPO A,
         VarArgs&&... vargs)
    {
    A.setPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

MPO
mapPrime(MPO A, int plevold, int plevnew, IndexSet const& is);

template <typename... VarArgs>
MPO
mapPrime(MPO A,
         VarArgs&&... vargs)
    {
    A.mapPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

MPO
swapPrime(MPO A, int plevold, int plevnew, IndexSet const& is);

template <typename... VarArgs>
MPO
swapPrime(MPO A,
         VarArgs&&... vargs)
    {
    A.swapPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

MPO
noPrime(MPO A, IndexSet const& is);

template <typename... VarArgs>
MPO
noPrime(MPO A,
        VarArgs&&... vargs)
    {
    A.noPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

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

// Re[Tr(A)]
Real
trace(MPO const& A);

// Tr(A)
Cplx
traceC(MPO const& A);

// Tr(A)
void
trace(MPO const& A,
      Real& re, Real& im);

// Re[Tr(AB)]
Real
trace(MPO const& A,
      MPO const& B);

// Tr(AB)
Cplx
traceC(MPO const& A,
       MPO const& B);

// Tr(AB)
void
trace(MPO const& A,
      MPO const& B,
      Real& re, Real& im);

// Calculate <x|A|y>
void 
inner(MPS const& x,
      MPO const& A, 
      MPS const& y, 
      Real& re, 
      Real& im);

// Calculate <x|A|y>
Real 
inner(MPS const& x, 
      MPO const& A, 
      MPS const& y);

// Calculate <x|A|y>
Complex 
innerC(MPS const& x,
       MPO const& A, 
       MPS const& y);

// Generic calculate <x|A|y> 
template <typename T> 
T 
innerT(MPS const& x,
       MPO const& A, 
       MPS const& y);

template <> inline
Real
innerT<Real>(MPS const& x,
             MPO const& A, 
             MPS const& y)
{
    return inner(x,A,y);
}

template <> inline
Complex
innerT<Complex>(MPS const& x,
                MPO const& A, 
                MPS const& y)
{
    return innerC(x,A,y);
}

      
// Calculate <Ax|By>
void
inner(MPO const& A, 
      MPS const& x,
      MPO const& B,
      MPS const& y, 
      Real& re, 
      Real& im);

// Calculate <Ax|By>
Real
inner(MPO const& A, 
      MPS const& x, 
      MPO const& B,
      MPS const& y);

// Calculate <Ax|By>
Complex
innerC(MPO const& A,
       MPS const& x, 
       MPO const& B,
       MPS const& y);

// Calculate <x|AB|y>
void
inner(MPS const& x,
      MPO const& A,
      MPO const& B,
      MPS const& y,
      Real& re,
      Real& im);

// Calculate <x|AB|y>
Real
inner(MPS const& x,
      MPO const& A,
      MPO const& B,
      MPS const& y);

// Calculate <x|AB|y>
Complex
innerC(MPS const& x,
       MPO const& A,
       MPO const& B,
       MPS const& y);

// Calculate AB
void 
nmultMPO(MPO const& Aorig, 
         MPO const& Borig, 
         MPO & res,
         Args args = Args::global());

// Calculate AB
MPO
nmultMPO(MPO const& A,
         MPO const& B,
         Args args = Args::global());

//
//{"Method=","DensityMatrix"}:
//Applies an MPO K to an MPS x (K|x>) with no approximation
//made in the application of K to x. Compresses
//the result back into an MPS whose bond dimension
//is at most the product of the bond dimension of K
//and the bond dimension of x. The result can 
//be controllably truncated further by providing
//optional truncation args "Cutoff" and "MaxDim"
//
//{"Method=","Fit"}
//Applies an MPO K to an MPS psi (|res>=K|psi>) using a sweeping/DMRG-like
//fitting approach. Warning: this method can get stuck i.e. fail to converge
//if the initial value of res is too different from the product K|psi>.
//List of options recognized:
//   Normalize (default: true) - normalize state to 1 after applying MPO
//   Nsweep (default: 1) - number of sweeps to use
//   MaxDim (default: res.maxdim()) - maximum number of states to keep
//   MinDim (default: res.mindim()) - minimum number of states to keep
//   Cutoff (default: res.cutoff()) - maximum truncation error goal
//
MPS
applyMPO(MPO const& K,
         MPS const& x,
         Args args = Args::global());

//Takes a starting guess wavefunction (only for {"Method=","Fit"})
MPS
applyMPO(MPO const& K,
         MPS const& x,
         MPS const& x0,
         Args args = Args::global());

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

//
// Deprecated
//


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

#ifdef ITENSOR_USE_HDF5
void
h5_write(h5::group parent, std::string const& name, MPO const& M);
void
h5_read(h5::group parent, std::string const& name, MPO & M);
#endif

} //namespace itensor

#endif
