//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MPS_IMPL_H_
#define __ITENSOR_MPS_IMPL_H_

namespace itensor {


template<class T>
SiteSet const& MPSt<T>::
sites() const 
    { 
    if(not sites_) Error("MPS SiteSet is default-initialized");
    return sites_; 
    }

template <class Tensor>
template <class BigMatrixT>
Spectrum MPSt<Tensor>::
svdBond(int b, const Tensor& AA, Direction dir, 
        const BigMatrixT& PH, const Args& args)
    {
    setBond(b);
    if(dir == Fromleft && b-1 > leftLim())
        {
        printfln("b=%d, l_orth_lim_=%d",b,leftLim());
        Error("b-1 > l_orth_lim_");
        }
    if(dir == Fromright && b+2 < rightLim())
        {
        printfln("b=%d, r_orth_lim_=%d",b,rightLim());
        Error("b+2 < r_orth_lim_");
        }

    auto noise = args.getReal("Noise",0.);
    auto cutoff = args.getReal("Cutoff",MIN_CUT);
    auto usesvd = args.getBool("UseSVD",false);

    Spectrum res;

    if(usesvd || (noise == 0 && cutoff < 1E-12))
        {
        //Need high accuracy, use svd which calls the
        //accurate SVD method in the MatrixRef library
        Tensor D;
        res = svd(AA,A_[b],D,A_[b+1],args);

        //Normalize the ortho center if requested
        if(args.getBool("DoNormalize",false))
            {
            D *= 1./itensor::norm(D);
            }

        //Push the singular values into the appropriate site tensor
        if(dir == Fromleft) A_[b+1] *= D;
        else                A_[b]   *= D;
        }
    else
        {
        //If we don't need extreme accuracy
        //or need to use noise term
        //use density matrix approach
        res = denmatDecomp(AA,A_[b],A_[b+1],dir,PH,args);

        //Normalize the ortho center if requested
        if(args.getBool("DoNormalize",false))
            {
            Tensor& oc = (dir == Fromleft ? A_[b+1] : A_[b]);
            auto nrm = itensor::norm(oc);
            if(nrm > 1E-16) oc *= 1./nrm;
            }
        }

    if(dir == Fromleft)
        {
        l_orth_lim_ = b;
        if(r_orth_lim_ < b+2) r_orth_lim_ = b+2;
        }
    else //dir == Fromright
        {
        if(l_orth_lim_ > b-1) l_orth_lim_ = b-1;
        r_orth_lim_ = b+1;
        }

    return res;
    }

template<typename T>
bool
isComplex(MPSt<T> const& psi)
    {
    for(auto j : range1(psi.N()))
        {
        if(itensor::isComplex(psi.A(j))) return true;
        }
    return false;
    }

template<typename T>
bool
isOrtho(MPSt<T> const& psi)
    {
    return psi.leftLim()+1 == psi.rightLim()-1;
    }

template<typename T>
int
orthoCenter(MPSt<T> const& psi)
    {
    if(!isOrtho(psi)) Error("orthogonality center not well defined.");
    return (psi.leftLim() + 1);
    }

template<typename T>
Real
norm(MPSt<T> const& psi)
    {
    if(not isOrtho(psi)) Error("\
MPS must have well-defined ortho center to compute norm; \
call .position(j) or .orthogonalize() to set ortho center");
    return itensor::norm(psi.A(orthoCenter(psi)));
    }

template<typename T>
Real
normalize(MPSt<T> & psi)
    {
    auto nrm = norm(psi);
    if(std::fabs(nrm) < 1E-20) Error("Zero norm");
    psi /= nrm;
    return nrm;
    }

template <typename MPST>
typename MPST::IndexT 
linkInd(MPST const& psi, int b)
    { 
    return commonIndex(psi.A(b),psi.A(b+1),Link); 
    }

template <typename MPST>
typename MPST::IndexT 
rightLinkInd(MPST const& psi, int i)
    { 
    return commonIndex(psi.A(i),psi.A(i+1),Link); 
    }

template <typename MPST>
typename MPST::IndexT 
leftLinkInd(MPST const& psi, int i)
    { 
    return commonIndex(psi.A(i),psi.A(i-1),Link); 
    }

template <typename MPST>
Real
averageM(MPST const& psi)
    {
    Real avgm = 0;
    for(int b = 1; b < psi.N(); ++b) 
        {
        avgm += linkInd(psi,b).m();
        }
    avgm /= (psi.N()-1);
    return avgm;
    }

template <typename MPST>
int
maxM(MPST const& psi)
    {
    int maxM_ = 0;
    for(int b = 1; b < psi.N(); ++b) 
        {
        int mb = commonIndex(psi.A(b),psi.A(b+1)).m();
        maxM_ = std::max(mb,maxM_);
        }
    return maxM_;
    }

template <class Tensor>
void 
applyGate(Tensor const& gate, 
          MPSt<Tensor>& psi,
          Args const& args)
    {
    auto fromleft = args.getBool("Fromleft",true);
    const int c = psi.orthoCenter();
    Tensor AA = psi.A(c) * psi.A(c+1) * gate;
    AA.noprime();
    if(fromleft) psi.svdBond(c,AA,Fromleft,args);
    else         psi.svdBond(c,AA,Fromright,args);
    }

template <class Tensor>
bool 
checkOrtho(const MPSt<Tensor>& psi,
           int i, 
           bool left)
    {
    using IndexT = typename Tensor::index_type;

    IndexT link = (left ? rightLinkInd(psi,i) : leftLinkInd(psi,i));
    Tensor rho = psi.A(i) * dag(prime(psi.A(i),link,4));
    Tensor Delta = delta(link, prime(link, 4));
    Tensor Diff = rho - Delta;

    const
    Real threshold = 1E-13;
    if(norm(Diff) < threshold) 
        {
        return true;
        }

    //Print any helpful debugging info here:
    println("checkOrtho: on line ",__LINE__," of mps.h,");
    println("checkOrtho: Tensor at position ",i," failed to be ",left?"left":"right"," ortho.");
    printfln("checkOrtho: norm(Diff) = %E",norm(Diff));
    printfln("checkOrtho: Error threshold set to %E",threshold);
    //-----------------------------

    return false;
    }

template <class Tensor>
bool 
checkOrtho(MPSt<Tensor> const& psi)
    {
    for(int i = 1; i <= psi.leftLim(); ++i)
    if(!checkOrtho(psi,i,true))
        {
        std::cout << "checkOrtho: A_[i] not left orthogonal at site i=" 
                  << i << std::endl;
        return false;
        }

    for(int i = psi.N(); i >= psi.rightLim(); --i)
    if(!checkOrtho(psi,i,false))
        {
        std::cout << "checkOrtho: A_[i] not right orthogonal at site i=" 
                  << i << std::endl;
        return false;
        }
    return true;
    }

template <class MPSType>
Cplx 
overlapC(MPSType const& psi, 
         MPSType const& phi)
    {
    auto N = psi.N();
    if(N != phi.N()) Error("overlap: mismatched N");

    auto l1 = linkInd(psi,1);
    auto L = phi.A(1);
    if(l1) L *= dag(prime(psi.A(1),l1)); 
    else   L *= dag(psi.A(1));

    if(N == 1) return L.cplx();

    for(decltype(N) i = 2; i < N; ++i) 
        { 
        L = L * phi.A(i) * dag(prime(psi.A(i),Link)); 
        }
    L = L * phi.A(N);

    auto lNm = linkInd(psi,N-1);
    if(lNm) return (dag(prime(psi.A(N),lNm))*L).cplx();
    return (dag(psi.A(N))*L).cplx();
    }

template <class MPSType>
void 
overlap(MPSType const& psi,MPSType const& phi, Real& re, Real& im)
    {
    auto z = overlapC(psi,phi);
    re = z.real();
    im = z.imag();
    }

template <class MPSType>
Real 
overlap(MPSType const& psi, MPSType const& phi) //Re[<psi|phi>]
    {
    Real re, im;
    overlap(psi,phi,re,im);
    if(std::fabs(im) > (1E-12 * std::fabs(re)) )
        printfln("Real overlap: WARNING, dropping non-zero imaginary part (=%.5E) of expectation value.",im);
    return re;
    }

template <class MPSType>
Complex 
psiphiC(MPSType const& psi, 
        MPSType const& phi)
    {
    return overlapC(psi,phi);
    }

template <class MPSType>
void 
psiphi(MPSType const& psi,MPSType const& phi, Real& re, Real& im)
    {
    overlap(psi,phi,re,im);
    }

template <class MPSType>
Real 
psiphi(MPSType const& psi, MPSType const& phi) //Re[<psi|phi>]
    {
    return overlap(psi,phi);
    }

template <class Tensor>
MPSt<Tensor>
sum(MPSt<Tensor> const& L, 
    MPSt<Tensor> const& R, 
    Args const& args)
    {
    MPSt<Tensor> res(L);
    res.plusEq(R,args);
    return res;
    }

template <typename MPSType>
MPSType 
sum(std::vector<MPSType> const& terms, 
    Args const& args)
    {
    auto Nt = terms.size();
    if(Nt == 2)
        { 
        return sum(terms.at(0),terms.at(1),args);
        }
    else 
    if(Nt == 1) 
        {
        return terms.at(0);
        }
    else 
    if(Nt > 2)
        {
        //Add all MPS in pairs
        auto nsize = (Nt%2==0 ? Nt/2 : (Nt-1)/2+1);
        std::vector<MPSType> newterms(nsize); 
        for(decltype(Nt) n = 0, np = 0; n < Nt-1; n += 2, ++np)
            {
            newterms.at(np) = sum(terms.at(n),terms.at(n+1),args);
            }
        if(Nt%2 == 1) newterms.at(nsize-1) = terms.back();

        //Recursively call sum again
        return sum(newterms,args);
        }
    return MPSType();
    }

} //namespace itensor

#endif
