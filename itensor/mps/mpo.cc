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
        Aref(j) = sites.op("Id",j);
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

//MPO
//sum(MPO L, 
//    MPO const& R, 
//    Args const& args)
//    {
//    L.plusEq(R,args);
//    return L;
//    }

void 
psiHphi(MPS const& psi, 
        MPO const& H, 
        MPS const& phi, 
        Real& re, 
        Real& im)
    {
    overlap(psi,H,phi,re,im);
    }

Real 
psiHphi(MPS const& psi, 
        MPO const& H, 
        MPS const& phi) //Re[<psi|H|phi>]
    {
    return overlap(psi,H,phi);
    }

Complex 
psiHphiC(MPS const& psi, 
         MPO const& H, 
         MPS const& phi) //Re[<psi|H|phi>]
    {
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
    overlap(psi,H,LB,RB,phi,re,im);
    }

Real
psiHphi(MPS const& psi, 
        MPO const& H, 
        ITensor const& LB, 
        ITensor const& RB, 
        MPS const& phi) //Re[<psi|H|phi>]
    {
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
    overlap(psi,H,K,phi,re,im);
    }

Real
psiHKphi(MPS const& psi, 
         MPO const& H, 
         MPO const& K,
         MPS const& phi) //<psi|H K|phi>
    {
    return overlap(psi,H,K,phi);
    }

Complex
psiHKphiC(MPS const& psi, 
          MPO const& H, 
          MPO const& K,
          MPS const& phi) //<psi|H K|phi>
    {
    return overlapC(psi,H,K,phi);
    }



int 
findCenter(MPO const& psi)
    {
    for(int j = 1; j <= length(psi); ++j) 
        {
        const auto& A = psi.A(j);
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
        if(!H.A(i))
            {
            println("A(",i,") null, QNs not well defined");
            Error("QNs not well defined");
            }
        if(div(H.A(i)) != Zero)
            {
            cout << "At i = " << i << endl;
            Print(H.A(i));
            Error("Non-zero div ITensor in MPO");
            }
        }

    //Check arrows from left edge
    for(int i = 1; i < center; ++i)
        {
        if(rightLinkInd(H,i).dir() != In) 
            {
            println("checkQNs: At site ",i," to the left of the OC, Right side Link not pointing In");
            Error("Incorrect Arrow in MPO");
            }
        if(i > 1)
            {
            if(leftLinkInd(H,i).dir() != Out) 
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
        if(rightLinkInd(H,i).dir() != Out) 
            {
            println("checkQNs: At site ",i," to the right of the OC, Right side Link not pointing Out");
            Error("Incorrect Arrow in MPO");
            }
        if(leftLinkInd(H,i).dir() != In) 
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
    for(int i = 1; i <= length(M); ++i) s << M.A(i) << "\n";
    return s;
    }

void
putMPOLinks(MPO& W, Args const& args)
    {
    if(not hasQNs(W.A(1)))
        {
        string pfix = args.getString("Prefix","MPO");
        auto links = vector<Index>(length(W));
        for(int b = 1; b < length(W); ++b)
            {
            links.at(b) = Index(1,format("Link,%s,%d",pfix,b));
            }
        W.Aref(1) *= setElt(links.at(1)(1));
        for(int b = 2; b < length(W); ++b)
            {
            W.Aref(b) *= setElt(links.at(b-1)(1));
            W.Aref(b) *= setElt(links.at(b)(1));
            }
        W.Aref(length(W)) *= setElt(links.at(length(W)-1)(1));
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
            q += div(W.A(b));
            links.at(b) = Index(q,1,Out,ts);
            }

        W.Aref(1) *= setElt(links.at(1)(1));
        for(int b = 2; b < N; ++b)
            {
            W.Aref(b) *= setElt(dag(links.at(b-1)(1)));
            W.Aref(b) *= setElt(links.at(b)(1));
            }
        W.Aref(N) *= setElt(dag(links.at(N-1)(1)));
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
//        res.Aref(j) = ITensor(K.A(j));
//        }
//    res.leftLim(K.leftLim());
//    res.rightLim(K.rightLim());
//    return res;
//    }

bool
isComplex(MPO const& W)
    {
    for(auto j : range1(length(W)))
        {
        if(itensor::isComplex(W.A(j))) return true;
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

    auto L = phi.A(1); 
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

    auto z = (dag(prime(psi.A(N)))*L).eltC();
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

    auto L = (LB ? LB*phi.A(1) : phi.A(1));
    L *= H.A(1); 
    L *= dag(prime(psi.A(1)));
    for(int i = 2; i <= N; ++i)
        { 
        L *= phi.A(i); 
        L *= H.A(i); 
        L *= dag(prime(psi.A(i))); 
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
        printfln("Real psiHphi: WARNING, dropping non-zero imaginary part (=%.5E) of expectation value.",im);
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
    if(length(psi) != length(phi) || length(psi) != length(H) || length(psi) != length(K)) Error("Mismatched N in psiHKphi");
    auto N = length(psi);
    auto psidag = psi;
    for(int i = 1; i <= N; i++)
        {
        psidag.Aref(i) = dag(prime(psi.A(i),2));
        }
    auto Hp = H;
    Hp.mapPrime(1,2);
    Hp.mapPrime(0,1);

    //scales as m^2 k^2 d
    auto L = phi.A(1) * K.A(1) * Hp.A(1) * psidag.A(1);
    //printfln("L%02d = %s",1,L);
    for(int i = 2; i < N; i++)
        {
        //scales as m^3 k^2 d + m^2 k^3 d^2
        L = L * phi.A(i) * K.A(i) * Hp.A(i) * psidag.A(i);
        //printfln("L%02d = %s",i,L);
        }
    //scales as m^2 k^2 d
    L = L * phi.A(N) * K.A(N) * Hp.A(N) * psidag.A(N);
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
	Error("Non-zero imaginary part in psiHKphi");
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
    err += -2.*overlapC(psi2,K,psi1).real();
    //Compute Kd, Hermitian conjugate of K
    auto Kd = K;
    for(auto j : range1(length(K)))
        {
        Kd.Aref(j) = dag(swapPrime(K.A(j),0,1,"Site"));
        }
    err /= overlap(psi1,Kd,K,psi1);
    err = std::sqrt(std::abs(1.0+err));
    return err;
    }

Real
checkMPOProd(MPS const& psi2,
             MPO const& K, 
             MPS const& psi1)
    {
    Global::warnDeprecated("checkMPOProd is deprecated in favor of errorMPOProd");
    //||p2> - K|p1>|^2 = (<p2|-<p1|Kd)(|p2>-K|p1>) = <p2|p2>+<p1|Kd*K|p1>-2*Re[<p2|K|p1>]
    Real res = overlap(psi2,psi2);
    res += -2.*overlapC(psi2,K,psi1).real();
    //Compute Kd, Hermitian conjugate of K
    auto Kd = K;
    for(auto j : range1(length(K)))
        {
        Kd.Aref(j) = dag(swapPrime(K.A(j),0,1,"Site"));
        }
    res += overlap(psi1,Kd,K,psi1);
    return res;
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

//template<class Tensor> 
//void MPO::
//position(int i, const Args& args)
//    {
//    if(isNull()) Error("position: MPS is null");
//
//    while(l_orth_lim_ < i-1)
//        {
//        if(l_orth_lim_ < 0) l_orth_lim_ = 0;
//        Tensor WF = A(l_orth_lim_+1) * A(l_orth_lim_+2);
//        svdBond(l_orth_lim_+1,WF,Fromleft,args);
//        }
//    while(r_orth_lim_ > i+1)
//        {
//        if(r_orth_lim_ > N_+1) r_orth_lim_ = N_+1;
//        Tensor WF = A(r_orth_lim_-2) * A(r_orth_lim_-1);
//        svdBond(r_orth_lim_-2,WF,Fromright,args);
//        }
//
//    is_ortho_ = true;
//    }
//template void MPO<ITensor>::
//position(int b, const Args& args);
//template void MPO<IQTensor>::
//position(int b, const Args& args);

//template <class Tensor>
//void MPO::
//orthogonalize(const Args& args)
//    {
//    //Do a half-sweep to the right, orthogonalizing each bond
//    //but do not truncate since the basis to the right might not
//    //be ortho (i.e. use the current m).
//    //svd_.useOrigM(true);
//    int orig_maxdim = maxdim();
//    Real orig_cutoff = cutoff();
//    for(Spectrum& spec : spectrum_)
//        {
//        spec.maxdim(MAX_DIM);
//        spec.cutoff(MIN_CUT);
//        }
//
//    position(1);
//    position(N_);
//
//    //Now basis is ortho, ok to truncate
//    for(Spectrum& spec : spectrum_)
//        {
//        spec.useOrigM(false);
//        spec.maxdim(orig_maxdim);
//        spec.cutoff(orig_cutoff);
//        }
//    position(1);
//
//    is_ortho_ = true;
//    }
//template
//void MPO<ITensor>::orthogonalize(const Args& args);
//template
//void MPO<IQTensor>::orthogonalize(const Args& args);

} //namespace itensor
