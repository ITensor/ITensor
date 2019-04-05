//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/util/print_macro.h"
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
    // Ref. norm is Tr{1} = d^N, d = 2 S=1/2, d = 4 for Hubbard, etc
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

Index
siteIndex(MPO const& W, int b, TagSet const& tsmatch)
    {
    return findIndex(siteInds(W,b),tsmatch);
    }

// Get the site Index of the MPS W*A 
// as if MPO W was applied to MPS A
Index
siteIndex(MPO const& W, MPS const& A, int b)
    {
    return uniqueIndex(W(b),{W(b-1),W(b+1),A(b)});
    }

IndexSet
siteInds(MPO const& A, MPS const& x)
    {
    auto N = length(x);
    auto inds = IndexSetBuilder(N);
    for( auto n : range1(N) )
      {
      auto s = siteIndex(A,x,n);
      inds.nextIndex(std::move(s));
      }
    return inds.build();
    }

// Get the site Indices of the MPO A*B 
// as if MPO A and MPO B were contracted
// TODO: implement
//Index inline
//siteInds(MPO const& A, MPO const& B, int b)
//    {
//    auto sA = uniqueIndex(A(b),{A(b-1),A(b+1),B(b)});
//    auto sB = uniqueIndex(B(b),{B(b-1),B(b+1),A(b)});
//    return IndexSet(sA,sB);
//    }

//Index
//siteIndex(MPO const& A, MPO const& B, int b, TagSet const& tsmatch = TagSet("0"))
//    {
//    return findIndex(siteInds(A,B,b),tsmatch);
//    }



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
    if(not hasQNs(W(1)))
        {
        string pfix = args.getString("Prefix","MPO");
        auto links = vector<Index>(length(W));
        for(int b = 1; b < length(W); ++b)
            {
            links.at(b) = Index(1,format("Link,%s,%d",pfix,b));
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
        auto pfix = args.getString("Prefix","l");

        auto links = vector<Index>(N);
        for(int b = 1; b < N; ++b)
            {
            string ts = format("%s,%d",pfix,b);
            q += div(W(b));
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

//<psi|H|phi>
void 
overlap(MPS const& psi, 
        MPO const& H, 
        MPS const& phi, 
        Real& re, 
        Real& im)
    {
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
    Real re,im;
    overlap(psi,H,K,phi,re,im);
    return Cplx(re,im);
    }

Real
errorMPOProd(MPS const& psi2,
             MPO const& K, 
             MPS const& psi1)
    {
    //||p2> - K|p1>| / || K|p1> || = sqrt{(<p2|-<p1|Kd)(|p2>-K|p1>) / <p1|KdK|p1>}
    //                             = sqrt{1+ (<p2|p2>-2*Re[<p2|K|p1>]) / <p1|KdK|p1>}
    Real err = overlap(psi2,psi2);
    err += -2.*real(overlapC(psi2,K,psi1));
    //Compute Kd, Hermitian conjugate of K
    auto Kd = K;
    for(auto j : range1(length(K)))
        {
        auto s = siteInds(Kd,j);
        Kd.ref(j) = dag(swapInds(Kd(j),{s(1)},{s(2)}));
        }
    err /= overlap(psi1,Kd,K,psi1);
    err = std::sqrt(std::abs(1.0+err));
    return err;
    }

bool
checkMPOProd(MPS const& psi2,
             MPO const& K, 
             MPS const& psi1,
             Real threshold)
    {
    Real err = errorMPOProd(psi2,K,psi1);
    return (std::norm(err) < threshold);
    }

//
// Deprecated
//

Real
checkMPOProd(MPS const& psi2,
             MPO const& K, 
             MPS const& psi1)
    {
    Global::warnDeprecated("checkMPOProd is deprecated in favor of errorMPOProd");
    //||p2> - K|p1>|^2 = (<p2|-<p1|Kd)(|p2>-K|p1>) = <p2|p2>+<p1|Kd*K|p1>-2*Re[<p2|K|p1>]
    Real res = overlap(psi2,psi2);
    res += -2.*real(overlapC(psi2,K,psi1));
    //Compute Kd, Hermitian conjugate of K
    auto Kd = K;
    for(auto j : range1(length(K)))
        {
        auto s = siteInds(K,j);
        Kd.ref(j) = dag(swapInds(K(j),{s(1)},{s(2)}));
        }
    res += overlap(psi1,Kd,K,psi1);
    return res;
    }

void 
psiHphi(MPS const& psi, 
        MPO const& H, 
        MPS const& phi, 
        Real& re, 
        Real& im)
    {
    Global::warnDeprecated("psiHphi deprecated in favor of overlap");
    overlap(psi,H,phi,re,im);
    }

Real 
psiHphi(MPS const& psi, 
        MPO const& H, 
        MPS const& phi) //Re[<psi|H|phi>]
    {
    Global::warnDeprecated("psiHphi deprecated in favor of overlap");
    return overlap(psi,H,phi);
    }

Complex 
psiHphiC(MPS const& psi, 
         MPO const& H, 
         MPS const& phi) //Re[<psi|H|phi>]
    {
    Global::warnDeprecated("psiHphiC deprecated in favor of overlapC");
    return overlapC(psi,H,phi);
    }

void
psiHphi(MPS const& psi, 
        MPO const& H, 
        ITensor const& LB, 
        ITensor const& RB, 
        MPS const& phi, 
        Real& re, 
        Real& im) //<psi|H|phi>
    {
    Global::warnDeprecated("psiHphi deprecated in favor of overlap");
    overlap(psi,H,LB,RB,phi,re,im);
    }

Real
psiHphi(MPS const& psi, 
        MPO const& H, 
        ITensor const& LB, 
        ITensor const& RB, 
        MPS const& phi) //Re[<psi|H|phi>]
    {
    Global::warnDeprecated("psiHphi deprecated in favor of overlap");
    return overlap(psi,H,LB,RB,phi);
    }

void
psiHKphi(MPS const& psi, 
         MPO const& H, 
         MPO const& K,
         MPS const& phi, 
         Real& re, 
         Real& im) //<psi|H K|phi>
    {
    Global::warnDeprecated("psiHKphi deprecated in favor of overlap");
    overlap(psi,H,K,phi,re,im);
    }

Real
psiHKphi(MPS const& psi, 
         MPO const& H, 
         MPO const& K,
         MPS const& phi) //<psi|H K|phi>
    {
    Global::warnDeprecated("psiHKphi deprecated in favor of overlap");
    return overlap(psi,H,K,phi);
    }

Complex
psiHKphiC(MPS const& psi, 
          MPO const& H, 
          MPO const& K,
          MPS const& phi) //<psi|H K|phi>
    {
    Global::warnDeprecated("psiHKphiC deprecated in favor of overlapC");
    return overlapC(psi,H,K,phi);
    }

} //namespace itensor
