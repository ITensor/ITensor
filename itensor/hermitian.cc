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
//#include <tuple>
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
diagHImpl(ITensor H, 
          ITensor& U, 
          ITensor& D,
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

    auto origdim = dim(H.inds().front());
    auto cutoff = args.getReal("Cutoff",0.);
    auto maxdim = args.getInt("MaxDim",origdim);
    auto mindim = args.getInt("MinDim",1);
    auto def_do_trunc = args.defined("Cutoff") || args.defined("MaxDim");
    auto do_truncate = args.getBool("Truncate",def_do_trunc);
    auto doRelCutoff = args.getBool("DoRelCutoff",true);
    auto absoluteCutoff = args.getBool("AbsoluteCutoff",false);
    auto showeigs = args.getBool("ShowEigs",false);
    auto itagset = getTagSet(args,"Tags","Link");

    // If no truncation is occuring, reset MaxDim
    // to the full matrix dimension
    if(!do_truncate)
        {
        if(args.defined("MaxDim"))
            {
            args.add("MaxDim",origdim);
            maxdim = origdim;
            }
        }

    if(not hasQNs(H))
        {
        if(order(H) != 2)
            {
            Print(order(H));
            Print(H);
            Error("Tensor has more than 2 indices in diag_hermitian");
            }

        auto i1 = H.inds().front();
        auto i2 = H.inds().back();

        auto active = (primeLevel(i1) < primeLevel(i2)) ? i1 : i2;

        auto pdiff = std::abs(primeLevel(i1)-primeLevel(i2));


#ifdef USESCALE
        //Depending on the sign of the scale, calling .toMatrix11NoScale 
        //yields a matrix proportional to either H or -H.
        //If H (scale().sign() > 0) then want to temporarily reverse 
        //the sign of the matrix when calling the diagonalization routine
        //to ensure eigenvalues are ordered from largest to smallest.
        if(H.scale().sign() < 0) H.scaleTo(H.scale()*(-1));
#endif

        //Do the diagonalization
        Vector DD;
        Mat<T> UU,iUU;
        auto R = toMatRefc<T>(H,active,prime(active));
        diagHermitian(R,UU,DD);
        conjugate(UU);

        //Truncate
        Real truncerr = 0.0;
        long m = DD.size();
        Real docut_lower = -1;
        Real docut_upper = -1;
        int ndegen = 1;
        if(do_truncate)
            {
            //if(DD(1) < 0) DD *= -1; //DEBUG
            tie(truncerr,docut_lower,docut_upper,ndegen) = truncate(DD,maxdim,mindim,cutoff,absoluteCutoff,doRelCutoff,args);
            m = DD.size();
            reduceCols(UU,m);
            }

        if(m > maxdim)
            {
            printfln("m > maxdim; m = %d, maxdim = %d",m,maxdim);
            Error("m > maxdim");
            }
        if(m > 50000)
            {
            printfln("WARNING: very large m = %d in ITensor diag_hermitian");
            }

        if(showeigs)
            {
            auto showargs = args;
            showargs.add("Cutoff",cutoff);
            showargs.add("MaxDim",maxdim);
            showargs.add("MinDim",mindim);
            showargs.add("Truncate",do_truncate);
            showargs.add("DoRelCutoff",doRelCutoff);
            showargs.add("AbsoluteCutoff",absoluteCutoff);
            showEigs(DD,truncerr,sqrt(H.scale()),showargs);
            }

        auto newmid = Index(m,itagset);

        U = ITensor({active,newmid},Dense<T>{move(UU.storage())}); 
        D = ITensor({prime(newmid,pdiff),newmid},DiagReal{DD.begin(),DD.end()},H.scale());

        if(not H.scale().isTooBigForReal())
            {
            DD *= H.scale().real0();
            }
        else
            {
            println("diag_hermitian: scale too big for Real, omitting from returned spectrum.");
            }

        return Spectrum{move(DD),{"Truncerr",truncerr}};
        }
    else  // With QNs
        {
        auto compute_qns = args.getBool("ComputeQNs",false);

        if(H.order() != 2)
            {
            Print(H.inds());
            Error("diag_hermitian requires order 2 input tensor");
            }
        
        auto i1 = H.inds().front();
        auto i2 = H.inds().back();
        auto ai = (primeLevel(i1) < primeLevel(i2)) ? i1 : i2;
        auto pdiff = std::abs(primeLevel(i1)-primeLevel(i2));

#ifdef DEBUG
        auto Zero = QN();
        if(div(H) != Zero)
            { 
            Print(H); 
            Error("diagHermitian currently only defined for block-diagonal IQTensors");
            }
#endif

        if(H.scale().sign() < 0) H.scaleTo(H.scale()*(-1));

        auto blocks = doTask(GetBlocks<T>{H.inds(),ai,prime(ai,pdiff)},H.store());
        auto Nblock = blocks.size();

        size_t totaldsize = 0,
               totalUsize = 0;
        for(auto b : range(Nblock))
            {
            totaldsize += nrows(blocks[b].M);
            totalUsize += nrows(blocks[b].M)*ncols(blocks[b].M);
            }

        auto Udata = vector<T>(totalUsize);
        auto Umats = vector<MatRef<T>>(Nblock);

        auto ddata = vector<Real>(totaldsize);
        auto dvecs = vector<VectorRef>(Nblock);

        auto alleig = stdx::reserve_vector<Real>(dim(ai));
        auto alleigqn = vector<EigQN>{};
        if(compute_qns) alleigqn = stdx::reserve_vector<EigQN>(dim(ai));

        //1. Diagonalize each ITensor within H.
        //   Store results in mmatrix and mvector.
        totaldsize = 0;
        totalUsize = 0;
        for(auto b : range(Nblock))
            {
            auto& M = blocks[b].M;
            auto& UU = Umats.at(b);
            auto& d =  dvecs.at(b);
            auto rM = nrows(M),
                 cM = ncols(M);

            d = makeVecRef(ddata.data()+totaldsize,rM);
            UU = makeMatRef(Udata.data()+totalUsize,rM*cM,rM,cM);

            diagHermitian(M,UU,d);
            conjugate(UU);

            alleig.insert(alleig.end(),d.begin(),d.end());
            if(compute_qns)
                {
                auto bi = blocks[b].i1;
                auto q = qn(ai,1+bi);
                for(auto eig : d)
                    {
                    alleigqn.emplace_back(eig,q);
                    }
                }
            totaldsize += rM;
            totalUsize += rM*cM;
            }


        //2. Truncate eigenvalues

        stdx::sort(alleig,std::greater<Real>{});
        if(compute_qns) stdx::sort(alleigqn,std::greater<EigQN>{});

        auto probs = Vector{move(alleig),VecRange{alleig.size()}};

        //Determine number of states to keep m
        long m = probs.size();
        Real truncerr = 0;
        Real docut_lower = -1;
        Real docut_upper = -1;
        int ndegen = 0;
        if(do_truncate)
            {
            tie(truncerr,docut_lower,docut_upper,ndegen) = truncate(probs,maxdim,mindim,cutoff,
                                                                    absoluteCutoff,doRelCutoff,args);

            m = probs.size();
            alleigqn.resize(m);
            }

        if(showeigs)
            {
            auto showargs = args;
            showargs.add("Cutoff",cutoff);
            showargs.add("MaxDim",maxdim);
            showargs.add("MinDim",mindim);
            showargs.add("Truncate",do_truncate);
            showargs.add("DoRelCutoff",doRelCutoff);
            showargs.add("AbsoluteCutoff",absoluteCutoff);
            showEigs(probs,truncerr,H.scale(),showargs);
            }

        if(m > maxdim)
            {
            printfln("m > maxdim; m = %d, maxdim = %d",m,maxdim);
            Error("m > maxdim");
            }
        if(m > 20000)
            {
            printfln("WARNING: very large m = %d in diag_hermitian",m);
            }

        //3. Truncate eigenvalues and eigenvectors of H

        //Form new Link Index with appropriate m's for each block
        Index::qnstorage iq;
        iq.reserve(Nblock);

        auto total_m = 0;
        for(auto b : range(Nblock))
            {
            auto& UU = Umats.at(b);
            auto& d = dvecs.at(b);
            auto& B = blocks[b];

            decltype(d.size()) this_m = 0;
            if(do_truncate)
                {
                //Keep all eigenvalues above docut_upper
                while(this_m < d.size() && 
                      total_m < m && 
                      d(this_m) > docut_upper)
                    {
                    ++this_m;
                    ++total_m;
                    }
                //Now check if there are any degenerate eigenvalues to keep
                //(ones above docut_upper)
                while(ndegen > 0 && 
                      this_m < d.size() && 
                      total_m < m && 
                      d(this_m) > docut_lower)
                    {
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

            d = subVector(d,0,this_m);
            UU = columns(UU,0,this_m);

            iq.emplace_back(qn(ai,1+B.i1),this_m);
            }

        if(iq.empty())
            {
            if(blocks.empty()) Error("No blocks in IQTensor svd");
            auto& B = blocks.front();
            iq.emplace_back(qn(ai,1+B.i1),1l);
            }

        auto d = Index(move(iq),-ai.dir(),itagset);

        auto Uis = IndexSet(dag(ai),dag(d));
        auto Dis = IndexSet(prime(d,pdiff),dag(d));

        auto Ustore = QDense<T>(Uis,QN());
        auto Dstore = QDiagReal(Dis);

        long n = 0;
        for(auto b : range(Nblock))
            {
            auto& B = blocks[b];
            auto& UU = Umats.at(b);
            auto& dv = dvecs.at(b);
            auto mm = ncols(UU);
            //Default-constructed B.M corresponds
            //to this_m==0 case above
            if(not B.M) continue;

            auto uind = Block(2);
            uind[0] = B.i1;
            uind[1] = n;
            auto pU = getBlock(Ustore,Uis,uind);
            assert(pU.data() != nullptr);
            assert(ai.blocksize0(B.i1) == long(nrows(UU)));
            auto Uref = makeMatRef(pU,nrows(UU),mm);
            Uref &= UU;

            auto dind = Block(2);
            dind[0] = n;
            dind[1] = n;
            auto pD = getBlock(Dstore,Dis,dind);
            assert(pD.data() != nullptr);
            auto Dref = makeVecRef(pD.data(),mm);
            Dref &= dv;

            ++n;
            }

        U = ITensor(Uis,move(Ustore));
        D = ITensor(Dis,move(Dstore),H.scale());

        if(H.scale().isTooBigForReal())
            {
            println("scale too big, omitting from reported eigenvalues");
            }
        else
            {
            probs *= H.scale().real0();
            }

        if(compute_qns)
            {
            auto qns = stdx::reserve_vector<QN>(alleigqn.size());
            for(auto& eq : alleigqn) qns.push_back(eq.qn);
            return Spectrum(move(probs),move(qns),{"Truncerr",truncerr});
            }

        return Spectrum{move(probs),{"Truncerr",truncerr}};
        }
    }

Spectrum
diag_hermitian(ITensor    H, 
               ITensor  & U, 
               ITensor  & D,
               Args const& args)
    {
    if(isComplex(H))
        {
        return diagHImpl<Cplx>(H,U,D,args);
        }
    return diagHImpl<Real>(H,U,D,args);
    }

} //namespace itensor
