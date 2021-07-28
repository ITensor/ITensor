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
#include <algorithm>
#include "itensor/util/stdx.h"
#include "itensor/tensor/algs.h"
#include "itensor/decomp.h"
#include "itensor/util/print_macro.h"
#include "itensor/itdata/qutil.h"

namespace itensor {

//const auto MAX_INT = std::numeric_limits<int>::max();

using std::swap;
using std::istream;
using std::ostream;
using std::vector;
using std::find;
using std::pair;
using std::make_pair;
using std::string;
using std::sqrt;
using std::move;
using std::tie;

template<typename T>
Spectrum
svdImpl(ITensor const& A,
        Index const& uI, 
        Index const& vI,
        ITensor & U, 
        ITensor & D, 
        ITensor & V,
        Args args)
    {
    if( args.defined("Minm") )
      {
      if( args.defined("MinDim") )
        {
        Global::warnDeprecated("Args Minm and MinDim are both defined. Minm is deprecated in favor of MinDim, MinDim will be used.");
        }
      else
        {
        Global::warnDeprecated("Arg Minm is deprecated in favor of MinDim.");
        args.add("MinDim",args.getInt("Minm"));
        }
      }

    if( args.defined("Maxm") )
      {
      if( args.defined("MaxDim") )
        {
        Global::warnDeprecated("Args Maxm and MaxDim are both defined. Maxm is deprecated in favor of MaxDim, MaxDim will be used.");
        }
      else
        {
        Global::warnDeprecated("Arg Maxm is deprecated in favor of MaxDim.");
        args.add("MaxDim",args.getInt("Maxm"));
        }
      }

    auto do_truncate = args.getBool("Truncate");
    auto cutoff = args.getReal("Cutoff",MIN_CUT);
    auto maxdim = args.getInt("MaxDim",MAX_DIM);
    auto mindim = args.getInt("MinDim",1);
    auto doRelCutoff = args.getBool("DoRelCutoff",true);
    auto absoluteCutoff = args.getBool("AbsoluteCutoff",false);
    auto show_eigs = args.getBool("ShowEigs",false);
    auto litagset = getTagSet(args,"LeftTags","Link,U");
    auto ritagset = getTagSet(args,"RightTags","Link,V");
    if(litagset == ritagset) 
        {
        Error("In SVD, must specify different tags for the new left and right indices (with Args 'LeftTags' and 'RightTags')");
        }

    if(not hasQNs(A))
        {
        auto M = toMatRefc<T>(A,uI,vI);

        Mat<T> UU,VV;
        Vector DD;

        SVD(M,UU,DD,VV,args);

        //conjugate VV so later we can just do
        //U*D*V to reconstruct ITensor A:
        conjugate(VV);

        //
        // Truncate
        //
        Vector probs;
        if(do_truncate || show_eigs)
            {
            probs = DD;
            for(auto j : range(probs)) probs(j) = sqr(probs(j));
            }

        Real truncerr = 0;
        Real docut_lower = -1;
        Real docut_upper = -1;
        int ndegen = 1;
        long m = DD.size();
        if(do_truncate)
            {
            tie(truncerr,docut_lower,docut_upper,ndegen) = truncate(probs,maxdim,mindim,cutoff,
                                                                    absoluteCutoff,doRelCutoff,args);
            m = probs.size();
            resize(DD,m);
            reduceCols(UU,m);
            reduceCols(VV,m);
            }


        if(show_eigs) 
            {
            auto showargs = args;
            showargs.add("Cutoff",cutoff);
            showargs.add("MaxDim",maxdim);
            showargs.add("MinDim",mindim);
            showargs.add("Truncate",do_truncate);
            showargs.add("DoRelCutoff",doRelCutoff);
            showargs.add("AbsoluteCutoff",absoluteCutoff);
            showEigs(probs,truncerr,A.scale(),showargs);
            }
        
        auto uL = Index(m,litagset);
        auto vL = Index(m,ritagset);

        //Fix sign to make sure D has positive elements
        Real signfix = (A.scale().sign() == -1) ? -1 : +1;
        D = ITensor({uL,vL},
                    Diag<Real>{DD.begin(),DD.end()},
                    A.scale()*signfix);
        U = ITensor({uI,uL},Dense<T>(move(UU.storage())),LogNum(signfix));
        V = ITensor({vI,vL},Dense<T>(move(VV.storage())));

        //Square all singular values
        //since convention is to report
        //density matrix eigs
        for(auto& el : DD) el = sqr(el);

#ifdef USESCALE
        if(A.scale().isFiniteReal()) 
            {
            DD *= sqr(A.scale().real0());
            }
        else                         
            {
            println("Warning: scale not finite real after svd");
            }
#endif

        return Spectrum(move(DD),{"Truncerr",truncerr});
        }
    else
        {
        auto compute_qn = args.getBool("ComputeQNs",false);

        auto blocks = doTask(GetBlocks<T>{A.inds(),uI,vI},A.store());

        auto Nblock = blocks.size();
        if(Nblock == 0) throw ResultIsZero("IQTensor has no blocks");

        //TODO: optimize allocation/lookup of Umats,Vmats
        //      etc. by allocating memory ahead of time (see algs.cc)
        //      and making Umats a vector of MatrixRef's to this memory
        auto Umats = vector<Mat<T>>(Nblock);
        auto Vmats = vector<Mat<T>>(Nblock);

        //TODO: allocate dvecs in a single allocation
        //      make dvecs a vector<VecRef>
        auto dvecs = vector<Vector>(Nblock);

        auto alleig = stdx::reserve_vector<Real>(std::min(dim(uI),dim(vI)));

        auto alleigqn = vector<EigQN>{};
        if(compute_qn)
            {
            alleigqn = stdx::reserve_vector<EigQN>(std::min(dim(uI),dim(vI)));
            }

        if(dim(uI) == 0) throw ResultIsZero("dim(uI) == 0");
        if(dim(vI) == 0) throw ResultIsZero("dim(vI) == 0");

        for(auto b : range(Nblock))
            {
            auto& M = blocks[b].M;
            auto& UU = Umats.at(b);
            auto& VV = Vmats.at(b);
            auto& d =  dvecs.at(b);

            SVD(M,UU,d,VV,args);

            //conjugate VV so later we can just do
            //U*D*V to reconstruct ITensor A:
            conjugate(VV);

            alleig.insert(alleig.end(),d.begin(),d.end());
            if(compute_qn)
                {
                auto bi = blocks[b].i1;
                auto q = qn(uI,1+bi);
                for(auto sval : d)
                    {
                    alleigqn.emplace_back(sqr(sval),q);
                    }
                }
            }

        //Square the singular values into probabilities
        //(density matrix eigenvalues)
        for(auto& sval : alleig) sval = sval*sval;

        //Sort all eigenvalues from largest to smallest
        //irrespective of quantum numbers
        stdx::sort(alleig,std::greater<Real>{});
        if(compute_qn) stdx::sort(alleigqn,std::greater<EigQN>{});

        auto probs = Vector(move(alleig),VecRange{alleig.size()});

        long m = probs.size();
        Real truncerr = 0;
        Real docut_lower = -1;
        Real docut_upper = -1;
        int ndegen = 1;
        if(do_truncate)
            {
            tie(truncerr,docut_lower,docut_upper,ndegen) = truncate(probs,maxdim,mindim,cutoff,
                                                                    absoluteCutoff,doRelCutoff,args);
            m = probs.size();
            alleigqn.resize(m);
            }

        if(show_eigs) 
            {
            auto showargs = args;
            showargs.add("Cutoff",cutoff);
            showargs.add("MaxDim",maxdim);
            showargs.add("MinDim",mindim);
            showargs.add("Truncate",do_truncate);
            showargs.add("DoRelCutoff",doRelCutoff);
            showargs.add("AbsoluteCutoff",absoluteCutoff);
            showEigs(probs,truncerr,A.scale(),showargs);
            }

        auto Liq = Index::qnstorage{};
        auto Riq = Index::qnstorage{};
        Liq.reserve(Nblock);
        Riq.reserve(Nblock);

        auto total_m = 0;
        for(auto b : range(Nblock))
            {
            auto& d = dvecs.at(b);
            auto& B = blocks[b];

            decltype(d.size()) this_m = 0;
            if(do_truncate)
                {
                //Keep all eigenvalues above docut_upper
                while(this_m < d.size() &&
                      total_m < m &&
                      sqr(d(this_m)) > docut_upper)
                    {
                    if(d(this_m) < 0) d(this_m) = 0;
                    ++this_m;
                    ++total_m;
                    }
                //Now check if there are any degenerate eigenvalues to keep
                //(ones above docut_lower)
                while(ndegen > 0 &&
                      this_m < d.size() &&
                      total_m < m &&
                      sqr(d(this_m)) > docut_lower)
                    {
                    if(d(this_m) < 0) d(this_m) = 0;
                    ++this_m;
                    ++total_m;
                    --ndegen;
                    }
                }
            else
                {
                this_m = d.size();
                total_m += this_m;
                }

            if(this_m == 0) 
                { 
                d.clear();
                B.M.clear();
                assert(not B.M);
                continue; 
                }

            resize(d,this_m);
            qn(uI,1+B.i1);
            Liq.emplace_back(qn(uI,1+B.i1),this_m);
            Riq.emplace_back(qn(vI,1+B.i2),this_m);
            }
        
#ifdef DEBUG
        if(Liq.empty()) throw std::runtime_error("New Index of U after SVD is empty");
        if(Riq.empty()) throw std::runtime_error("New Index of V after SVD is empty");
#endif

        auto L = Index(move(Liq),uI.dir(),litagset);
        auto R = Index(move(Riq),vI.dir(),ritagset);

        auto Uis = IndexSet(uI,dag(L));
        auto Dis = IndexSet(L,R);
        auto Vis = IndexSet(vI,dag(R));

        auto Ustore = QDense<T>(Uis,QN());
        auto Vstore = QDense<T>(Vis,QN());
        auto Dstore = QDiagReal(Dis);

        long n = 0;
        for(auto b : range(Nblock))
            {
            auto& B = blocks[b];
            auto& UU = Umats.at(b);
            auto& VV = Vmats.at(b);
            auto& d = dvecs.at(b);
            //Default-constructed B.M corresponds
            //to this_m==0 case above
            if(not B.M) continue;

            //println("block b = ",b);
            //printfln("{B.i1,n} = {%d,%d}",B.i1,n);
            //printfln("{n,n} = {%d,%d}",n,n);
            //printfln("{B.i2,n} = {%d,%d}",B.i2,n);

            auto uind = Block(2);
            uind[0] = B.i1;
            uind[1] = n;
            auto pU = getBlock(Ustore,Uis,uind);
            assert(pU.data() != nullptr);
            assert(uI.blocksize0(B.i1) == long(nrows(UU)));
            auto Uref = makeMatRef(pU,uI.blocksize0(B.i1),L.blocksize0(n));
            reduceCols(UU,L.blocksize0(n));
            Uref &= UU;

            auto dind = Block(2);
            dind[0] = n;
            dind[1] = n;
            auto pD = getBlock(Dstore,Dis,dind);
            assert(pD.data() != nullptr);
            auto Dref = makeVecRef(pD.data(),d.size());
            Dref &= d;

            auto vind = Block(2);
            vind[0] = B.i2;
            vind[1] = n;
            auto pV = getBlock(Vstore,Vis,vind);
            assert(pV.data() != nullptr);
            assert(vI.blocksize0(B.i2) == long(nrows(VV)));
            auto Vref = makeMatRef(pV.data(),pV.size(),vI.blocksize0(B.i2),R.blocksize0(n));
            reduceCols(VV,R.blocksize0(n));
            //println("Doing Vref &= VV");
            //Print(Vref.range());
            //Print(VV.range());
            Vref &= VV;

            /////////DEBUG
            //Matrix D(d.size(),d.size());
            //for(decltype(d.size()) n = 0; n < d.size(); ++n)
            //    {
            //    D(n,n) = d(n);
            //    }
            //D *= A.scale().real0();
            //auto AA = Uref * D * transpose(Vref);
            //Print(Uref);
            //Print(D);
            //Print(Vref);
            //printfln("Check %d = \n%s",b,AA);
            //printfln("Diff %d = %.10f",b,norm(AA-B.M));
            /////////DEBUG

            ++n;
            }

        //Fix sign to make sure D has positive elements
        Real signfix = (A.scale().sign() == -1) ? -1. : +1.;
        U = ITensor(Uis,move(Ustore));
        D = ITensor(Dis,move(Dstore),A.scale()*signfix);
        V = ITensor(Vis,move(Vstore),LogNum{signfix});
        
        //Originally eigs were found without including scale
        //so put the scale back in
        if(A.scale().isFiniteReal())
            {
            probs *= sqr(A.scale().real0());
            }
        else
            {
            println("Warning: scale not finite real after svd");
            }

        if(compute_qn)
            {
            auto qns = stdx::reserve_vector<QN>(alleigqn.size());
            for(auto& eq : alleigqn) qns.push_back(eq.qn);
            return Spectrum(move(probs),move(qns),{"Truncerr",truncerr});
            }

        return Spectrum(move(probs),{"Truncerr",truncerr});
        }
    return Spectrum{};
    }



Spectrum 
svdOrd2(ITensor const& A, 
        Index const& uI, 
        Index const& vI,
        ITensor & U, 
        ITensor & D, 
        ITensor & V,
        Args args)
    {
    if( args.defined("Maxm") )
      {
      if( args.defined("MaxDim") )
        {
        Global::warnDeprecated("Args Maxm and MaxDim are both defined. Maxm is deprecated in favor of MaxDim, MaxDim will be used.");
        }
      else
        {
        Global::warnDeprecated("Arg Maxm is deprecated in favor of MaxDim.");
        args.add("MaxDim",args.getInt("Maxm"));
        }
      }

    auto do_truncate = args.defined("Cutoff") || args.defined("MaxDim");
    if(not args.defined("Truncate")) 
        args.add("Truncate",do_truncate);

    if(A.order() != 2) 
        {
        Error("A must be matrix-like (order 2)");
        }
    if(isComplex(A))
        {
        return svdImpl<Cplx>(A,uI,vI,U,D,V,args);
        }
    return svdImpl<Real>(A,uI,vI,U,D,V,args);
    }

} //namespace itensor
