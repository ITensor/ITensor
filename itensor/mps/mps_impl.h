//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MPS_IMPL_H_
#define __ITENSOR_MPS_IMPL_H_

namespace itensor {

template <typename BigMatrixT>
Spectrum MPS::
svdBond(int b, ITensor const& AA, Direction dir, 
        BigMatrixT const& PH, Args args)
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

    // Store the original tags for link b so that it can
    // be put back onto the newly introduced link index
    auto original_link_tags = tags(linkIndex(*this,b));

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

    // Put the old tags back onto the new index
    auto lb = commonIndex(A_[b],A_[b+1]);
    A_[b].setTags(original_link_tags,lb);
    A_[b+1].setTags(original_link_tags,lb);


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
    for(auto j : range1(length(psi)))
        {
        if(itensor::isComplex(psi(j))) return true;
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
    return itensor::norm(psi(orthoCenter(psi)));
    }

Real inline
normalize(MPS & psi)
    {
    Global::warnDeprecated("normalize(MPS) is deprecated in favor of .normalize()");
    auto nrm = psi.normalize();
    return nrm;
    }

template <typename MPSType>
IndexSet 
siteInds(MPSType const& W, int b)
    {
    return uniqueInds(W(b),{W(b-1),W(b+1)});
    }

Index inline
siteIndex(MPS const& psi, int j)
    { 
    return uniqueIndex(psi(j),{psi(j-1),psi(j+1)}); 
    }

template <typename MPSType>
Index
rightLinkIndex(MPSType const& psi, int i)
    {
    if( i == length(psi) ) return Index();
    return commonIndex(psi(i),psi(i+1));
    }

// Note: for ITensors with QNs, this is different
// from rightLinkIndex(psi,i-1) since indices 
// evaluate equal even in arrow directions are different
template <typename MPSType>
Index
leftLinkIndex(MPSType const& psi, int i)
    {
    if( i == 1 ) return Index();
    return commonIndex(psi(i),psi(i-1));
    }

// This is a shorthand for rightLinkIndex
template <typename MPSType>
Index 
linkIndex(MPSType const& psi, int i)
    {
    return rightLinkIndex(psi,i);
    }

template <typename MPSType>
IndexSet
linkInds(MPSType const& psi, int i)
    {
    if( i == 1 ) return IndexSet(rightLinkIndex(psi,i));
    else if( i == length(psi) ) return IndexSet(leftLinkIndex(psi,i));
    return IndexSet(leftLinkIndex(psi,i),rightLinkIndex(psi,i));
    }

Real inline
averageLinkDim(MPS const& psi)
    {
    Real avgdim = 0;
    for(int b = 1; b < length(psi); ++b)
        {
        avgdim += dim(linkIndex(psi,b));
        }
    avgdim /= (length(psi)-1);
    return avgdim;
    }

Real inline
averageM(MPS const& psi)
    {
    Global::warnDeprecated("averageM(MPS) is deprecated in favor of averageLinkDim(MPS)");
		return averageLinkDim(psi);
    }

int inline
maxLinkDim(MPS const& psi)
    {
    int maxdim_ = 0;
    for(int b = 1; b < length(psi); ++b)
        {
        int mb = dim(commonIndex(psi(b),psi(b+1)));
        maxdim_ = std::max(mb,maxdim_);
        }
    return maxdim_;
    }

int inline
maxM(MPS const& psi)
    {
    Global::warnDeprecated("maxM(MPS) is deprecated in favor of maxLinkDim(MPS)");
    return maxLinkDim(psi);
    }

template <typename MPSType>
void
overlap(MPSType const& psi, MPSType const& phi, Real& re, Real& im)
    {
    auto z = overlapC(psi,phi);
    re = real(z);
    im = imag(z);
    }

template <typename MPSType>
Real
overlap(MPSType const& psi, MPSType const& phi) //Re[<psi|phi>]
    {
    Real re, im;
    overlap(psi,phi,re,im);
    if(std::fabs(im) > (1E-12 * std::fabs(re)) )
        printfln("Real overlap: WARNING, dropping non-zero imaginary part (=%.5E) of expectation value.",im);
    return re;
    }

namespace detail
  {

  ITensor inline
  denseDelta(Index const& i, Index const& j)
      {
      auto del = ITensor(i,j);
      for( auto ii : range1(minDim(del)) )
        del.set(ii,ii,1);
      return del;
      }

  }

template <typename MPSType>
bool
checkOrtho(MPSType const& A,
           int i,
           Direction dir,
           Real threshold)
    {
    auto left = (dir==Fromleft);

    auto Adag = A;
    Adag.dag().primeLinks(1);

    auto lout_dag = (left ? leftLinkIndex(Adag,i) : rightLinkIndex(Adag,i));
    auto rho = A(i) * prime(Adag(i),-1,lout_dag);

    auto lin = (left ? rightLinkIndex(A,i) : leftLinkIndex(A,i));
    auto lin_dag = (left ? rightLinkIndex(Adag,i) : leftLinkIndex(Adag,i));
    auto Delta = detail::denseDelta(lin, lin_dag);

    auto Diff = rho - Delta;

    if(norm(Diff) < threshold)
        return true;

    //Print any helpful debugging info here:
    println("checkOrtho: on line ",__LINE__," of mps.h,");
    println("checkOrtho: Tensor at position ",i," failed to be ",left?"left":"right"," ortho.");
    printfln("checkOrtho: norm(Diff) = %E",norm(Diff));
    printfln("checkOrtho: Error threshold set to %E",threshold);
    //-----------------------------

    return false;
    }

template <typename MPSType>
bool
checkOrtho(MPSType const& A,
           Real threshold)
    {
    for(int i = 1; i <= A.leftLim(); ++i)
    if(!checkOrtho(A,i,Fromleft,threshold))
        {
        std::cout << "checkOrtho: A_[i] not left orthogonal at site i="
                  << i << std::endl;
        return false;
        }

    for(int i = length(A); i >= A.rightLim(); --i)
    if(!checkOrtho(A,i,Fromright,threshold))
        {
        std::cout << "checkOrtho: A_[i] not right orthogonal at site i="
                  << i << std::endl;
        return false;
        }
    return true;
    }

Complex inline
psiphiC(MPS const& psi, 
        MPS const& phi)
    {
    Global::warnDeprecated("psiphiC(MPS,MPS) is deprecated in favor of overlapC(MPS,MPS)");
    return overlapC(psi,phi);
    }

void inline
psiphi(MPS const& psi, MPS const& phi, Real& re, Real& im)
    {
    Global::warnDeprecated("psiphi(MPS,MPS,Real,Real) is deprecated in favor of overlap(MPS,MPS,Real,Real)");
    overlap(psi,phi,re,im);
    }

Real inline
psiphi(MPS const& psi, MPS const& phi) //Re[<psi|phi>]
    {
    Global::warnDeprecated("psiphi(MPS,MPS) is deprecated in favor of overlap(MPS,MPS)");
    return overlap(psi,phi);
    }


} //namespace itensor

#endif
