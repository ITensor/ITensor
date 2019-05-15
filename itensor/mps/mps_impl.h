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
    // Truncate blocks of degenerate singular values
    args.add("RespectDegenerate",args.getBool("RespectDegenerate",true));

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

} //namespace itensor

#endif
