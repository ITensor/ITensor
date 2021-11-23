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
#include "itensor/util/print_macro.h"
#include "itensor/util/str.h"
#include "itensor/mps/mpo.h"
#include "itensor/mps/localop.h"

namespace itensor {

using std::istream;
using std::ostream;
using std::cout;
using std::endl;
using std::vector;
using std::find;
using std::pair;
using std::make_pair;
using std::string;

MPO::
MPO() 
  : Parent(),
    logrefNorm_(DefaultLogRefScale)
    { 
    }

MPO::
MPO(int N)
  : Parent(N),
    logrefNorm_(DefaultLogRefScale)
    { 
    }

MPO::
MPO(const SiteSet& sites,
     Real _logrefNorm) 
    : 
    Parent(sites)
    { 
    // Norm of psi^2 = 1 = norm = sum of denmat evals. 
    // This translates to Tr{Adag A} = norm.  
    // Ref. norm is Tr{1} = d^N, d = 2 S=1/2, d = 4 for Electron, etc
    if(_logrefNorm == DefaultLogRefScale) logrefNorm_ = sites.length();

    //Set all tensors to identity ops
    for(int j = 1; j <= length(); ++j)
        {
        ref(j) = sites.op("Id",j);
        }
    putMPOLinks(*this);
    }


MPO& MPO::
plusEq(MPO const& other_,
       Args const& args)
    {
    if(doWrite())
        Error("operator+= not supported if doWrite(true)");

    //cout << "calling new orthog in sum" << endl;
    if(!itensor::isOrtho(*this))
        {
        try { 
            orthogonalize(); 
            }
        catch(const ResultIsZero& rz) 
            { 
            *this = other_;
            return *this;
            }
        }

    if(!itensor::isOrtho(other_))
        {
        auto other = other_;
        try { 
            other.orthogonalize(); 
            }
        catch(const ResultIsZero& rz) 
            { 
            return *this;
            }
        return addAssumeOrth(*this,other,args);
        }

    return addAssumeOrth(*this,other_,args);
    }

MPO&
operator*=(MPO & W, Real a) { W.ref(W.leftLim()+1) *= a; return W; }

MPO&
operator*=(MPO & W, Cplx a) { W.ref(W.leftLim()+1) *= a; return W; }

MPO&
operator/=(MPO & W, Real a) { W.ref(W.leftLim()+1) /= a; return W; }

MPO&
operator/=(MPO & W, Cplx a) { W.ref(W.leftLim()+1) /= a; return W; }

MPO
operator*(MPO W, Real r) { return W *= r; }

MPO
operator*(Real r, MPO W) { return W *= r; }

MPO
operator*(MPO W, Cplx z) { return W *= z; }

MPO
operator*(Cplx z, MPO W) { return W *= z; }

MPO
dag(MPO W)
    {
    W.dag();
    return W;
    }

int
length(MPO const& W)
    {
    return W.length();
    }

bool
isOrtho(MPO const& W)
    {
    return W.leftLim()+1 == W.rightLim()-1;
    }

int
orthoCenter(MPO const& W)
    {
    if(!isOrtho(W)) Error("orthogonality center not well defined.");
    return (W.leftLim() + 1);
    }

bool
hasSiteInds(MPO const& A, IndexSet const& sites)
    {
    auto N = length(A);
    if( N!=length(sites) ) Error("In hasSiteInds(MPO,IndexSet), lengths of MPO and IndexSet of site indices don't match");
    for( auto n : range1(N) )
      {
      if( !hasIndex(A(n),sites(n)) ) return false;
      }
    return true;
    }

// Find the site Index of the bth MPO tensor of W
// having the tags tsmatch
Index
siteIndex(MPO const& W, int b, TagSet const& tsmatch)
    {
    return findIndex(siteInds(W,b),tsmatch);
    }

// Find the site Index of the bth MPO tensor of W
// that is not in the IndexSet is
Index
uniqueSiteIndex(MPO const& W, IndexSet const& is, int b)
    {
    return findIndex(uniqueInds(siteInds(W,b),is));
    }

// Get the site Index of the MPS W*A 
// as if MPO W was applied to MPS A
Index
uniqueSiteIndex(MPO const& W, MPS const& A, int b)
    {
    return uniqueIndex(W(b),{W(b-1),W(b+1),A(b)});
    }

// Get the site Index that is unique to A
Index
uniqueSiteIndex(MPO const& A, MPO const& B, int b)
    {
    return uniqueIndex(A(b),{A(b-1),A(b+1),B(b)});
    }

IndexSet
uniqueSiteInds(MPO const& A, MPS const& x)
    {
    auto N = length(x);
    if( N!=length(x) ) Error("In uniqueSiteInds(MPO,MPS), lengths of MPO and MPS do not match");
    auto inds = IndexSetBuilder(N);
    for( auto n : range1(N) )
      {
      auto s = uniqueSiteIndex(A,x,n);
      inds.nextIndex(std::move(s));
      }
    return inds.build();
    }

// Get the site Indices of the MPO A*B 
// as if MPO A and MPO B were contracted
IndexSet
siteInds(MPO const& A, MPO const& B, int b)
    {
    auto sA = uniqueSiteIndex(A,B,b);
    auto sB = uniqueSiteIndex(B,A,b);
    return IndexSet(sA,sB);
    }

// Get the site Indices that are unique to A
IndexSet
uniqueSiteInds(MPO const& A, MPO const& B)
    {
    auto N = length(A);
    if( N!=length(B) ) Error("In uniqueSiteInds(MPO,MPO), lengths of MPO and MPS do not match");
    auto inds = IndexSetBuilder(N);
    for( auto n : range1(N) )
      {
      auto s = uniqueSiteIndex(A,B,n);
      inds.nextIndex(std::move(s));
      }
    return inds.build();
    }

// Get the site Indices that are unique to A
// (on A but not in the input IndexSet of site indices)
IndexSet
uniqueSiteInds(MPO const& A, IndexSet const& sites)
    {
    auto N = length(A);
    if( N!=length(sites) ) Error("In uniqueSiteInds(MPO,IndexSet), lengths of MPO and IndexSet do not match");
    auto inds = IndexSetBuilder(N);
    for( auto n : range1(N) )
      {
      auto sn = uniqueSiteIndex(A,{sites(n)},n);
      inds.nextIndex(std::move(sn));
      }
    return inds.build();
    }

MPO& MPO::
replaceSiteInds(IndexSet const& sites_old, IndexSet const& sites_new)
    {
    auto& A = *this;
    auto N = itensor::length(A);
    if( itensor::length(sites_new)!=N ) Error("In replaceSiteInds(MPO,IndexSet,IndexSet), number of new sites must be equal length of MPO");
    if( itensor::hasSiteInds(A,sites_new) ) return A;
    for( auto n : range1(N) )
        A_[n].replaceInds({sites_old(n)},{sites_new(n)});
    return A;
    }

MPO
replaceSiteInds(MPO A, IndexSet const& sites_old, IndexSet const& sites_new)
    {
    A.replaceSiteInds(sites_old,sites_new);
    return A;
    }

MPO& MPO::
swapSiteInds()
    {
    auto& A = *this;
    auto N = itensor::length(A);
    for( auto n : range1(N) )
        {
        auto s = itensor::siteInds(A,n);
        A_[n].swapInds({s(1)},{s(2)});
        }
    return A;
    }

MPO
swapSiteInds(MPO A)
    {
    A.swapSiteInds();
    return A;
    }

int 
findCenter(MPO const& psi)
    {
    for(int j = 1; j <= length(psi); ++j) 
        {
        const auto& A = psi(j);
        if(A.order() == 0) Error("Zero order tensor in MPO");
        bool allOut = true;
        for(const auto& I : A.inds())
            {
            //Only look at Link IQIndices
            if(!hasTags(I,"Link")) continue;

            if(I.dir() != Out)
                {
                allOut = false;
                break;
                }
            }

        //Found the ortho. center
        if(allOut) return j;
        }
    return -1;
    }

void
checkQNs(MPO const& H)
    {
    const int N = length(H);

    const QN Zero;

    int center = findCenter(H);
    if(center == -1)
        {
        Error("Did not find an ortho. center");
        }

    //Check that all ITensors have zero div
    //including the ortho. center
    for(int i = 1; i <= N; ++i) 
        {
        if(!H(i))
            {
            println("A(",i,") null, QNs not well defined");
            Error("QNs not well defined");
            }
        if(div(H(i)) != Zero)
            {
            cout << "At i = " << i << endl;
            Print(H(i));
            Error("Non-zero div ITensor in MPO");
            }
        }

    //Check arrows from left edge
    for(int i = 1; i < center; ++i)
        {
        if(dir(linkIndex(H,i)) != In) 
            {
            println("checkQNs: At site ",i," to the left of the OC, Right side Link not pointing In");
            Error("Incorrect Arrow in MPO");
            }
        if(i > 1)
            {
            if(dir(linkIndex(H,i-1)) != Out) 
                {
                println("checkQNs: At site ",i," to the left of the OC, Left side Link not pointing Out");
                Error("Incorrect Arrow in MPO");
                }
            }
        }

    //Check arrows from right edge
    for(int i = N; i > center; --i)
        {
        if(i < N)
        if(dir(linkIndex(H,i)) != Out) 
            {
            println("checkQNs: At site ",i," to the right of the OC, Right side Link not pointing Out");
            Error("Incorrect Arrow in MPO");
            }
        if(dir(linkIndex(H,i-1)) != In) 
            {
            println("checkQNs: At site ",i," to the right of the OC, Left side Link not pointing In");
            Error("Incorrect Arrow in MPO");
            }
        }
    }



std::ostream& 
operator<<(std::ostream& s, MPO const& M)
    {
    s << "\n";
    for(int i = 1; i <= length(M); ++i) s << M(i) << "\n";
    return s;
    }

void
putMPOLinks(MPO& W, Args const& args)
    {
    string pfix = args.getString("Prefix","l=");
    if(not hasQNs(W))
        {
        auto links = vector<Index>(length(W));
        for(int b = 1; b < length(W); ++b)
            {
            string ts = "Link,"+pfix+str(b);
            links.at(b) = Index(1,ts);
            }
        W.ref(1) *= setElt(links.at(1)(1));
        for(int b = 2; b < length(W); ++b)
            {
            W.ref(b) *= setElt(links.at(b-1)(1));
            W.ref(b) *= setElt(links.at(b)(1));
            }
        W.ref(length(W)) *= setElt(links.at(length(W)-1)(1));
        }
    else
        {
        QN q;
        auto N = length(W);

        auto links = vector<Index>(N);
        for(int b = 1; b < N; ++b)
            {
            q += div(W(b));
            string ts = "Link,"+pfix+str(b);
            links.at(b) = Index(q,1,Out,ts);
            }

        W.ref(1) *= setElt(links.at(1)(1));
        for(int b = 2; b < N; ++b)
            {
            W.ref(b) *= setElt(dag(links.at(b-1)(1)));
            W.ref(b) *= setElt(links.at(b)(1));
            }
        W.ref(N) *= setElt(dag(links.at(N-1)(1)));
        }
    }


//MPO
//toMPO(MPO const& K)
//    {
//    int N = length(K);
//    MPO res;
//    if(K.sites()) res = MPO(K.sites());
//    else          res = MPO(N);
//    res.logRefNorm(K.logRefNorm());
//    for(int j = 0; j <= N+1; ++j)
//        {
//        res.ref(j) = ITensor(K(j));
//        }
//    res.leftLim(K.leftLim());
//    res.rightLim(K.rightLim());
//    return res;
//    }

MPO
setTags(MPO A, TagSet const& ts, IndexSet const& is)
    {
    A.setTags(ts,is);
    return A;
    }

MPO
noTags(MPO A, IndexSet const& is)
    {
    A.noTags(is);
    return A;
    }

MPO
addTags(MPO A, TagSet const& ts, IndexSet const& is)
    {
    A.addTags(ts,is);
    return A;
    }

MPO
removeTags(MPO A, TagSet const& ts, IndexSet const& is)
    {
    A.removeTags(ts,is);
    return A;
    }

MPO
replaceTags(MPO A, TagSet const& ts1, TagSet const& ts2, IndexSet const& is)
    {
    A.replaceTags(ts1,ts2,is);
    return A;
    }

MPO
swapTags(MPO A, TagSet const& ts1, TagSet const& ts2, IndexSet const& is)
    {
    A.swapTags(ts1,ts2,is);
    return A;
    }

MPO
prime(MPO A, int plev, IndexSet const& is)
    {
    A.prime(plev,is);
    return A;
    }

MPO
prime(MPO A, IndexSet const& is)
    {
    A.prime(is);
    return A;
    }

MPO
setPrime(MPO A, int plev, IndexSet const& is)
    {
    A.setPrime(plev,is);
    return A;
    }

MPO
noPrime(MPO A, IndexSet const& is)
    {
    A.noPrime(is);
    return A;
    }

bool
isComplex(MPO const& W)
    {
    for(auto j : range1(length(W)))
        {
        if(itensor::isComplex(W(j))) return true;
        }
    return false;
    }

Cplx
traceC(MPO const& A)
    {
    auto N = length(A);
    auto trA_n = A(1) * delta(dag(siteInds(A,1)));
    auto L = trA_n;
    if(N == 1) return eltC(L);
    for(auto n : range1(2,N) )
        {
        trA_n = A(n) * delta(dag(siteInds(A,n)));
        L *= trA_n;
        }
    return eltC(L);
    }

void
trace(MPO const& A,
      Real& re, Real& im)
    {
    auto z = traceC(A);
    re = real(z);
    im = imag(z);
    }

Real
trace(MPO const& A)
    {
    Real re, im;
    trace(A,re,im);
    if(std::fabs(im) > (1E-12 * std::fabs(re)) )
        printfln("Real inner: WARNING, dropping non-zero imaginary part (=%.5E) of expectation value.",im);
    return re;
    }

Cplx
traceC(MPO const& A,
       MPO const& B)
    {
    auto N = length(A);
    if(N != length(B)) Error("traceC(MPO,MPO): mismatched N");

    // Make the site indices of the MPOs match
    // and the links not match
    auto sA = uniqueSiteInds(A,B);
    auto sB = uniqueSiteInds(B,A);
    auto Bp = replaceSiteInds(B,sB,sA);
    Bp.replaceLinkInds(sim(linkInds(Bp)));

    auto L = A(1) * Bp(1);
    if(N == 1) return eltC(L);
    for(auto i : range1(2,N) )
        L = L * A(i) * Bp(i);
    return eltC(L);
    }

void
trace(MPO const& A,
      MPO const& B,
      Real& re, Real& im)
    {
    auto z = traceC(A,B);
    re = real(z);
    im = imag(z);
    }

Real
trace(MPO const& A, MPO const& B) //Re[<psi|phi>]
    {
    Real re, im;
    trace(A,B,re,im);
    if(std::fabs(im) > (1E-12 * std::fabs(re)) )
        printfln("Real inner: WARNING, dropping non-zero imaginary part (=%.5E) of expectation value.",im);
    return re;
    }

//<x|A|y>
void 
inner(MPS const& x, 
      MPO const& A, 
      MPS const& y, 
      Real& re, 
      Real& im)
    {
    auto N = length(A);
    if( length(y) != N || length(x) != N ) Error("inner: mismatched N");

    // Make the indices of |x> and A|y> match
    auto sAy = uniqueSiteInds(A,y);
    auto xp = replaceSiteInds(x,sAy);

    // Dagger x, since it is the ket
    auto xdag = dag(xp);
    xdag.replaceLinkInds(sim(linkInds(xdag)));

    auto L = y(1) * A(1) * xdag(1);

    // TODO: some MPOs may store edge tensors
    // in A(0) and A(N+1). Add this back?
    //L *= (A(0) ? A(0)*A(1) : A(1));

    for( auto n : range1(2,N) ) 
        L = L * y(n) * A(n) * xdag(n); 

    // TODO: some MPOs may store edge tensors
    // in A(0) and A(N+1). Add this back?
    //if(A(N+1)) L *= A(N+1);

    auto z = eltC(L);
    re = real(z);
    im = imag(z);
    }

Real 
inner(MPS const& psi, 
      MPO const& H, 
      MPS const& phi) //Re[<psi|H|phi>]
    {
    if(isComplex(psi) || isComplex(H) || isComplex(phi)) Error("Cannot use inner(...) with complex MPS/MPO, use innerC(...) instead");
    Real re, im;
    inner(psi,H,phi,re,im);
    return re;
    }


Cplx 
innerC(MPS const& psi, 
       MPO const& H, 
       MPS const& phi) //Re[<psi|H|phi>]
    {
    Real re, im;
    inner(psi,H,phi,re,im);
    return Cplx(re,im);
    }

// Calculate <Ax|By>
void
inner(MPO const& A,
      MPS const& x,
      MPO const& B,
      MPS const& y,
      Real& re,
      Real& im)
  {
  if(length(x) != length(y) || length(x) != length(A) || length(y) != length(B)) Error("Mismatched N in inner");
  auto N = length(y);

  // Automatically match site indices
  auto Ap = replaceSiteInds(A,uniqueSiteInds(A,x),uniqueSiteInds(B,y));

  // Prime the links to avoid clashes
  auto Adag = dag(A);
  Adag.replaceLinkInds(sim(linkInds(Adag)));
  auto xdag = dag(x);
  xdag.replaceLinkInds(sim(linkInds(xdag)));

  //scales as m^2 k^2 d
  auto L = y(1) * B(1) * Adag(1) * xdag(1);
  for(int i = 2; i < N; i++)
      {
      //scales as m^3 k^2 d + m^2 k^3 d^2
      L = L * y(i) * B(i) * Adag(i) * xdag(i);
      }
  //scales as m^2 k^2 d
  L = L * y(N) * B(N) * Adag(N) * xdag(N);
  auto z = eltC(L);
  re = real(z);
  im = imag(z);
  }

// Calculate <Ax|By>
Real
inner(MPO const& A,
      MPS const& x,
      MPO const& B,
      MPS const& y)
    {
    if(isComplex(A) || isComplex(x) || isComplex(B) || isComplex(y)) Error("Cannot use inner(...) with complex MPS/MPO, use innerC(...) instead");
    Real re,im;
    inner(A,x,B,y,re,im);
    return re;
    }

// Calculate <Ax|By>
Complex
innerC(MPO const& A,
       MPS const& x,
       MPO const& B,
       MPS const& y)
    {
    Real re,im;
    inner(A,x,B,y,re,im);
    return Cplx(re,im);
    }

void
inner(MPS const& x, 
      MPO const& A, 
      MPO const& B,
      MPS const& y, 
      Real& re, 
      Real& im)
    {
    if(length(x) != length(y) || length(x) != length(A) || length(x) != length(B)) Error("Mismatched N in inner");
    auto N = length(x);

    // Assume order of operations A(B|y>), use replaceInds
    // to handle the case where A and B share all indices
    auto sABy = uniqueSiteInds(A,uniqueSiteInds(B,y)); 
    auto sAByp = sim(sABy);
    auto Ap = replaceSiteInds(A,sABy,sAByp);
    Ap.replaceLinkInds(sim(linkInds(Ap)));
    auto xp = replaceSiteInds(x,sAByp);
    auto xdag = dag(xp);
    xdag.replaceLinkInds(sim(linkInds(xdag)));

    //scales as m^2 k^2 d
    auto L = y(1) * B(1) * Ap(1) * xdag(1);
    for(int i = 2; i < N; i++)
        {
        //scales as m^3 k^2 d + m^2 k^3 d^2
        L = L * y(i) * B(i) * Ap(i) * xdag(i);
        }
    //scales as m^2 k^2 d
    L = L * y(N) * B(N) * Ap(N) * xdag(N);
    auto z = eltC(L);
    re = real(z);
    im = imag(z);
    }

Real
inner(MPS const& psi, 
      MPO const& H, 
      MPO const& K,
      MPS const& phi) //<psi|H K|phi>
    {
    if(isComplex(psi) || isComplex(H) || isComplex(K) || isComplex(phi)) Error("Cannot use inner(...) with complex MPS/MPO, use innerC(...) instead");
    Real re,im;
    inner(psi,H,K,phi,re,im);
    return re;
    }

Cplx
innerC(MPS const& psi, 
       MPO const& H, 
       MPO const& K,
       MPS const& phi) //<psi|H K|phi>
    {
    Real re,im;
    inner(psi,H,K,phi,re,im);
    return Cplx(re,im);
    }

// Check how close and approximation to A|x>, called |y>,
// is to the exact A|x>
//||y> - A|x>| / || A|x> || = sqrt{(<y|-<x|Ad)(|y>-A|x>) / <x|AdA|x>}
//                             = sqrt{1 + (<y|y>-2*Re[<y|A|x>]) / <x|AdA|x>}
Real
errorMPOProd(MPS const& y,
             MPO const& A, 
             MPS const& x)
    {
    auto err = real(innerC(y,y));
    err += -2.*real(innerC(y,A,x));
    err /= real(innerC(A,x,A,x));
    err = std::sqrt(std::abs(1.0+err));
    return err;
    }

//
// Deprecated
//

void 
overlap(MPS const& psi, 
        MPO const& H, 
        MPS const& phi, 
        Real& re, 
        Real& im)
    {
    Global::warnDeprecated("overlap is deprecated in favor of inner/trace");
    auto N = length(H);
    if(length(phi) != N || length(psi) != N) Error("psiHphi: mismatched N");

    auto L = phi(1); 
    //Some Hamiltonians may store edge tensors in H(0) and H(N+1)
    L *= (H(0) ? H(0)*H(1) : H(1));
    L *= dag(prime(psi(1)));
    for(int i = 2; i < N; ++i) 
        { 
        L *= phi(i); 
        L *= H(i); 
        L *= dag(prime(psi(i))); 
        }
    L *= phi(N); 
    L *= H(N);
    if(H(N+1)) L *= H(N+1);

    auto z = (dag(prime(psi(N)))*L).eltC();
    re = z.real();
    im = z.imag();
    }

Real 
overlap(MPS const& psi, 
        MPO const& H, 
        MPS const& phi) //Re[<psi|H|phi>]
    {
    Global::warnDeprecated("overlap is deprecated in favor of inner/trace");
    Real re, im;
    overlap(psi,H,phi,re,im);
    if(std::fabs(im) > 1E-5 * std::fabs(re) || std::fabs(im) > 1E-9)
        {
        printfln("\nReal psiHphi: WARNING, dropping non-zero (=%.5E) imaginary part of expectation value.",im);
        }
    return re;
    }


Cplx 
overlapC(MPS const& psi, 
         MPO const& H, 
         MPS const& phi) //Re[<psi|H|phi>]
    {
    Global::warnDeprecated("overlap is deprecated in favor of inner/trace");
    Real re, im;
    overlap(psi,H,phi,re,im);
    return Cplx(re,im);
    }

void
overlap(MPS const& psi, 
        MPO const& H, 
        ITensor const& LB, 
        ITensor const& RB, 
        MPS const& phi, 
        Real& re, 
        Real& im) //<psi|H|phi>
    {
    Global::warnDeprecated("overlap is deprecated in favor of inner/trace");
    auto N = length(psi);
    if(N != length(phi) || length(H) < N) Error("mismatched N in psiHphi");

    auto L = (LB ? LB*phi(1) : phi(1));
    L *= H(1); 
    L *= dag(prime(psi(1)));
    for(int i = 2; i <= N; ++i)
        { 
        L *= phi(i); 
        L *= H(i); 
        L *= dag(prime(psi(i))); 
        }

    if(RB) L *= RB;

    auto z = eltC(L);
    re = z.real();
    im = z.imag();
    }

Real
overlap(MPS const& psi, 
        MPO const& H, 
        ITensor const& LB, 
        ITensor const& RB, 
        MPS const& phi) //Re[<psi|H|phi>]
    {
    Global::warnDeprecated("overlap is deprecated in favor of inner/trace");
    Real re,im; 
    overlap(psi,H,LB,RB,phi,re,im);
    if(std::fabs(im) > 1.0e-12 * std::fabs(re))
        printfln("Real overlap: WARNING, dropping non-zero imaginary part (=%.5E) of expectation value.",im);
    return re;
    }

void
overlap(MPS const& psi, 
        MPO const& H, 
        MPO const& K,
        MPS const& phi, 
        Real& re, 
        Real& im)
    {
    Global::warnDeprecated("overlap is deprecated in favor of inner/trace");
    //println("Running psiHKphi");
    if(length(psi) != length(phi) || length(psi) != length(H) || length(psi) != length(K)) Error("Mismatched N in overlap");
    auto N = length(psi);
    auto psidag = psi;
    psidag.dag().prime(2);
    auto Hp = H;
    Hp.replaceTags("1","2").replaceTags("0","1");

    //scales as m^2 k^2 d
    auto L = phi(1) * K(1) * Hp(1) * psidag(1);
    //printfln("L%02d = %s",1,L);
    for(int i = 2; i < N; i++)
        {
        //scales as m^3 k^2 d + m^2 k^3 d^2
        L = L * phi(i) * K(i) * Hp(i) * psidag(i);
        //printfln("L%02d = %s",i,L);
        }
    //scales as m^2 k^2 d
    L = L * phi(N) * K(N) * Hp(N) * psidag(N);
    //PrintData(L);
    auto z = L.eltC();
    re = z.real();
    im = z.imag();
    }

Real
overlap(MPS const& psi, 
        MPO const& H, 
        MPO const& K,
        MPS const& phi) //<psi|H K|phi>
    {
    Global::warnDeprecated("overlap is deprecated in favor of inner/trace");
    Real re,im;
    overlap(psi,H,K,phi,re,im);
    if(std::fabs(im) > 1.0e-12 * std::fabs(re))
        Error("Non-zero imaginary part in overlap, use overlapC instead.");
    return re;
    }

Cplx
overlapC(MPS const& psi, 
         MPO const& H, 
         MPO const& K,
         MPS const& phi) //<psi|H K|phi>
    {
    Global::warnDeprecated("overlap is deprecated in favor of inner/trace");
    Real re,im;
    overlap(psi,H,K,phi,re,im);
    return Cplx(re,im);
    }

#ifdef ITENSOR_USE_HDF5

void
h5_write(h5::group parent, string const& name, MPO const& M)
    {
    auto g = parent.create_group(name);
    h5_write_attribute(g,"type","MPO",true);
    h5_write_attribute(g,"version",long(1));
    h5_write(g,"length",long(M.length()));
    h5_write(g,"rlim",long(M.rightLim()));
    h5_write(g,"llim",long(M.leftLim()));
    for(auto n : range1(M.length()))
        {
        h5_write(g,format("MPO[%d]",n),M(n));
        }
    }

void
h5_read(h5::group parent, string const& name, MPO & M)
    {
    auto g = parent.open_group(name);
    auto type = h5_read_attribute<string>(g,"type");
    if(type != "MPO") Error("Group does not contain MPO data in HDF5 file");
    auto N = h5_read<long>(g,"length");
    auto rlim = h5_read<long>(g,"rlim");
    auto llim = h5_read<long>(g,"llim");
    M = MPO(N);
    for(auto n : range1(N))
        {
        M.ref(n) = h5_read<ITensor>(g,format("MPO[%d]",n));
        }
    M.leftLim(llim);
    M.rightLim(rlim);
    }

#endif //ITENSOR_USE_HDF5

} //namespace itensor
