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

template <class Tensor>
MPOt<Tensor>::
MPOt() 
    : 
    Parent(),
    logrefNorm_(DefaultLogRefScale)
    { 
    }
template MPOt<ITensor>::MPOt();
template MPOt<IQTensor>::MPOt();

template <class Tensor>
MPOt<Tensor>::
MPOt(const SiteSet& sites,
     Real _logrefNorm) 
    : 
    Parent(sites)
    { 
    // Norm of psi^2 = 1 = norm = sum of denmat evals. 
    // This translates to Tr{Adag A} = norm.  
    // Ref. norm is Tr{1} = d^N, d = 2 S=1/2, d = 4 for Hubbard, etc
    if(_logrefNorm == DefaultLogRefScale) logrefNorm_ = sites.N();

    //Set all tensors to identity ops
    for(int j = 1; j <= N(); ++j)
        {
        Anc(j) = sites.op("Id",j);
        }
    putMPOLinks(*this);
    }
template
MPOt<ITensor>::
MPOt(const SiteSet& sites, Real _logrefNorm);
template
MPOt<IQTensor>::
MPOt(const SiteSet& sites, Real _logrefNorm);

/*
template<class Tensor> 
void MPOt<Tensor>::
position(int i, const Args& args)
    {
    if(isNull()) Error("position: MPS is null");

    while(l_orth_lim_ < i-1)
        {
        if(l_orth_lim_ < 0) l_orth_lim_ = 0;
        Tensor WF = A(l_orth_lim_+1) * A(l_orth_lim_+2);
        svdBond(l_orth_lim_+1,WF,Fromleft,args);
        }
    while(r_orth_lim_ > i+1)
        {
        if(r_orth_lim_ > N_+1) r_orth_lim_ = N_+1;
        Tensor WF = A(r_orth_lim_-2) * A(r_orth_lim_-1);
        svdBond(r_orth_lim_-2,WF,Fromright,args);
        }

    is_ortho_ = true;
    }
template void MPOt<ITensor>::
position(int b, const Args& args);
template void MPOt<IQTensor>::
position(int b, const Args& args);
*/

/*
template <class Tensor>
void MPOt<Tensor>::
orthogonalize(const Args& args)
    {
    //Do a half-sweep to the right, orthogonalizing each bond
    //but do not truncate since the basis to the right might not
    //be ortho (i.e. use the current m).
    //svd_.useOrigM(true);
    int orig_maxm = maxm();
    Real orig_cutoff = cutoff();
    for(Spectrum& spec : spectrum_)
        {
        spec.maxm(MAX_M);
        spec.cutoff(MIN_CUT);
        }

    position(1);
    position(N_);

    //Now basis is ortho, ok to truncate
    for(Spectrum& spec : spectrum_)
        {
        spec.useOrigM(false);
        spec.maxm(orig_maxm);
        spec.cutoff(orig_cutoff);
        }
    position(1);

    is_ortho_ = true;
    }
template
void MPOt<ITensor>::orthogonalize(const Args& args);
template
void MPOt<IQTensor>::orthogonalize(const Args& args);
*/


template <class Tensor>
MPOt<Tensor>& MPOt<Tensor>::
plusEq(const MPOt<Tensor>& other_,
       const Args& args)
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
template
MPOt<ITensor>& MPOt<ITensor>::plusEq(const MPOt<ITensor>& other, const Args&);
template
MPOt<IQTensor>& MPOt<IQTensor>::plusEq(const MPOt<IQTensor>& other, const Args&);

int 
findCenter(const IQMPO& psi)
    {
    for(int j = 1; j <= psi.N(); ++j) 
        {
        const IQTensor& A = psi.A(j);
        if(A.r() == 0) Error("Zero rank tensor in IQMPO");
        bool allOut = true;
        for(const IQIndex& I : A.inds())
            {
            //Only look at Link IQIndices
            if(I.type() != Link) continue;

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
checkQNs(const IQMPO& H)
    {
    const int N = H.N();

    const QN Zero;

    int center = findCenter(H);
    if(center == -1)
        {
        Error("Did not find an ortho. center");
        }

    //Check that all IQTensors have zero div
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
            Error("Non-zero div IQTensor in IQMPO");
            }
        }

    //Check arrows from left edge
    for(int i = 1; i < center; ++i)
        {
        if(rightLinkInd(H,i).dir() != In) 
            {
            println("checkQNs: At site ",i," to the left of the OC, Right side Link not pointing In");
            Error("Incorrect Arrow in IQMPO");
            }
        if(i > 1)
            {
            if(leftLinkInd(H,i).dir() != Out) 
                {
                println("checkQNs: At site ",i," to the left of the OC, Left side Link not pointing Out");
                Error("Incorrect Arrow in IQMPO");
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
            Error("Incorrect Arrow in IQMPO");
            }
        if(leftLinkInd(H,i).dir() != In) 
            {
            println("checkQNs: At site ",i," to the right of the OC, Left side Link not pointing In");
            Error("Incorrect Arrow in IQMPO");
            }
        }
    }



template <class Tensor>
std::ostream& 
operator<<(std::ostream& s, MPOt<Tensor> const& M)
    {
    s << "\n";
    for(int i = 1; i <= M.N(); ++i) s << M.A(i) << "\n";
    return s;
    }
template
std::ostream& 
operator<<(std::ostream& s, MPOt<ITensor> const& M);
template
std::ostream& 
operator<<(std::ostream& s, MPOt<IQTensor> const& M);

void
putMPOLinks(MPO& W, Args const& args)
    {
    const string pfix = args.getString("Prefix","l");
    vector<Index> links(W.N());
    for(int b = 1; b < W.N(); ++b)
        {
        links.at(b) = Index(format("%s%d",pfix,b));
        }
    W.Anc(1) *= links.at(1)(1);
    for(int b = 2; b < W.N(); ++b)
        {
        W.Anc(b) *= links.at(b-1)(1);
        W.Anc(b) *= links.at(b)(1);
        }
    W.Anc(W.N()) *= links.at(W.N()-1)(1);
    }

void
putMPOLinks(IQMPO& W, Args const& args)
    {
    QN q;
    const int N = W.N();
    const string pfix = args.getString("Prefix","l");

    vector<IQIndex> links(N);
    for(int b = 1; b < N; ++b)
        {
        string nm = format("%s%d",pfix,b);
               
        q += div(W.A(b));
        links.at(b) = IQIndex(nm,Index(nm),q);
        }

    W.Anc(1) *= links.at(1)(1);
    for(int b = 2; b < N; ++b)
        {
        W.Anc(b) *= dag(links.at(b-1)(1));
        W.Anc(b) *= links.at(b)(1);
        }
    W.Anc(N) *= dag(links.at(N-1)(1));
    }

} //namespace itensor
