//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MPS_IMPL_H_
#define __ITENSOR_MPS_IMPL_H_

namespace itensor {

template <class BigMatrixT>
Spectrum MPS::
svdBond(int b, ITensor const& AA, Direction dir, 
        BigMatrixT const& PH, Args const& args)
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
        ITensor D;
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
            ITensor& oc = (dir == Fromleft ? A_[b+1] : A_[b]);
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

bool inline
isComplex(MPS const& psi)
    {
    for(auto j : range1(psi.N()))
        {
        if(itensor::isComplex(psi.A(j))) return true;
        }
    return false;
    }

bool inline
isOrtho(MPS const& psi)
    {
    return psi.leftLim()+1 == psi.rightLim()-1;
    }

int inline
orthoCenter(MPS const& psi)
    {
    if(!isOrtho(psi)) Error("orthogonality center not well defined.");
    return (psi.leftLim() + 1);
    }

Real inline
norm(MPS const& psi)
    {
    if(not isOrtho(psi)) Error("\
MPS must have well-defined ortho center to compute norm; \
call .position(j) or .orthogonalize() to set ortho center");
    return itensor::norm(psi.A(orthoCenter(psi)));
    }

Real inline
normalize(MPS & psi)
    {
    auto nrm = norm(psi);
    if(std::fabs(nrm) < 1E-20) Error("Zero norm");
    psi /= nrm;
    return nrm;
    }

template<typename MPSType>
Index 
linkInd(MPSType const& psi, int b)
    { 
    return commonIndex(psi.A(b),psi.A(b+1),"Link"); 
    }

template<typename MPSType>
Index
rightLinkInd(MPSType const& psi, int i)
    { 
    return commonIndex(psi.A(i),psi.A(i+1),"Link"); 
    }

template<typename MPSType>
Index
leftLinkInd(MPSType const& psi, int i)
    { 
    return commonIndex(psi.A(i),psi.A(i-1),"Link"); 
    }

Real inline
averageM(MPS const& psi)
    {
    Real avgm = 0;
    for(int b = 1; b < psi.N(); ++b) 
        {
        avgm += linkInd(psi,b).m();
        }
    avgm /= (psi.N()-1);
    return avgm;
    }

int inline
maxM(MPS const& psi)
    {
    int maxM_ = 0;
    for(int b = 1; b < psi.N(); ++b) 
        {
        int mb = commonIndex(psi.A(b),psi.A(b+1)).m();
        maxM_ = std::max(mb,maxM_);
        }
    return maxM_;
    }

void inline
overlap(MPS const& psi, MPS const& phi, Real& re, Real& im)
    {
    auto z = overlapC(psi,phi);
    re = z.real();
    im = z.imag();
    }

Real inline
overlap(MPS const& psi, MPS const& phi) //Re[<psi|phi>]
    {
    Real re, im;
    overlap(psi,phi,re,im);
    if(std::fabs(im) > (1E-12 * std::fabs(re)) )
        printfln("Real overlap: WARNING, dropping non-zero imaginary part (=%.5E) of expectation value.",im);
    return re;
    }

Complex inline
psiphiC(MPS const& psi, 
        MPS const& phi)
    {
    return overlapC(psi,phi);
    }

void inline
psiphi(MPS const& psi, MPS const& phi, Real& re, Real& im)
    {
    overlap(psi,phi,re,im);
    }

Real inline
psiphi(MPS const& psi, MPS const& phi) //Re[<psi|phi>]
    {
    return overlap(psi,phi);
    }


} //namespace itensor

#endif
