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
#include <tuple>
#include <limits>
#include "itensor/util/stdx.h"
#include "itensor/tensor/algs.h"
#include "itensor/decomp.h"
#include "itensor/util/print_macro.h"
#include "itensor/itdata/qutil.h"

namespace itensor {

//const auto MAX_INT = std::numeric_limits<int>::max();
const auto REAL_EPSILON = std::numeric_limits<Real>::epsilon();

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

///////////////

template<typename V>
struct ToMatRefc
    {
    using value_type = V;
    long nrows=0,
         ncols=0;
    bool transpose=false;
    ToMatRefc(long nr, long nc, bool trans=false) 
        : nrows(nr), ncols(nc), transpose(trans)
        { }
    };
template<typename V>
MatRefc<V>
doTask(ToMatRefc<V> const& T, 
       Dense<V> const& d)
    {
    auto res = makeMatRef(d.data(),d.size(),T.nrows,T.ncols);
    if(T.transpose) return transpose(res);
    return res;
    }

template<typename V>
MatRefc<V>
toMatRefc(ITensor const& T, 
          Index const& i1, 
          Index const& i2)
    {
    if(i1 == T.inds().front())
        {
        return doTask(ToMatRefc<V>{dim(i1),dim(i2)},T.store());
        }
    return doTask(ToMatRefc<V>{dim(i2),dim(i1),true},T.store());
    }
template MatRefc<Real>
toMatRefc(ITensor const& T, Index const& i1, Index const& i2);
template MatRefc<Cplx>
toMatRefc(ITensor const& T, Index const& i1, Index const& i2);

/////////////

template<typename T>
vector<Ord2Block<T>>
doTask(GetBlocks<T> const& G, 
       QDense<T> const& d)
    {
    if(G.is.order() != 2) Error("doTask(GetBlocks,QDenseReal) only supports 2-index tensors");
    auto res = vector<Ord2Block<T>>{d.offsets.size()};
    auto dblock = IntArray(2,0);
    size_t n = 0;
    for(auto& dio : d.offsets)
        {
        auto& R = res[n++];
        computeBlockInd(dio.block,G.is,dblock);
        auto nrow = G.is[0].blocksize0(dblock[0]);
        auto ncol = G.is[1].blocksize0(dblock[1]);
        R.i1 = dblock[0];
        R.i2 = dblock[1];
        R.M = makeMatRef(d.data()+dio.offset,d.size()-dio.offset,nrow,ncol);
        }
    if(G.transpose) 
        {
        for(auto& R : res) 
            {
            R.M = transpose(R.M);
            std::swap(R.i1,R.i2);
            }
        }
    return res;
    }
template vector<Ord2Block<Real>>
doTask(GetBlocks<Real> const& G, QDense<Real> const& d);
template vector<Ord2Block<Cplx>>
doTask(GetBlocks<Cplx> const& G, QDense<Cplx> const& d);

///////////////

Spectrum
svd(ITensor const& AA,
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

    if( args.defined("UseOrigM") )
      {
      if( args.defined("UseOrigDim") )
        {
        Global::warnDeprecated("Args UseOrigM and UseOrigDim are both defined. UseOrigM is deprecated in favor of UseOrigDim, UseOrigDim will be used.");
        }
      else
        {
        Global::warnDeprecated("Arg UseOrigDim is deprecated in favor of MaxDim.");
        args.add("UseOrigDim",args.getBool("UseOrigM"));
        }
      }

#ifdef DEBUG
    if(!U && !V)
        Error("U and V default-initialized in svd, must indicate at least one index on U or V");
#endif

    auto noise = args.getReal("Noise",0);
    auto useOrigDim = args.getBool("UseOrigDim",false);

    // Set the cutoff to 0 by default
    args.add("Cutoff",args.getReal("Cutoff",0.0));

    if(noise > 0)
        Error("Noise term not implemented for svd");

    //if(isZero(AA,Args("Fast"))) 
    //    throw ResultIsZero("svd: AA is zero");


    //Combiners which transform AA
    //into a order 2 tensor
    std::vector<Index> Uinds,
                       Vinds;
    Uinds.reserve(AA.order());
    Vinds.reserve(AA.order());
    //Divide up indices based on U
    //If U is null, use V instead
    auto &L = (U ? U : V);
    auto &Linds = (U ? Uinds : Vinds),
         &Rinds = (U ? Vinds : Uinds);
    for(const auto& I : AA.inds())
        {
        if(hasIndex(L,I)) Linds.push_back(I);
        else              Rinds.push_back(I);
        }
    ITensor Ucomb,
            Vcomb;
    Index ui,
          vi;

    auto AAcomb = AA;
    if(!Uinds.empty())
        {
        std::tie(Ucomb,ui) = combiner(std::move(Uinds));
        AAcomb *= Ucomb;
        }
    if(!Vinds.empty())
        {
        std::tie(Vcomb,vi) = combiner(std::move(Vinds));
        AAcomb *= Vcomb;
        }

    if(useOrigDim)
        {
        //Try to determine current m,
        //then set mindim_ and maxdim_ to this.
        args.add("Cutoff",-1);
        long mindim = 1,
             maxdim = MAX_DIM;
        if(D.order() == 0)
            {
            //auto mid = commonIndex(U,V,Link);
            //TODO: check this does the same thing
            auto mid = commonIndex(U,V);
            if(mid) mindim = maxdim = dim(mid);
            else    mindim = maxdim = 1;
            }
        else
            {
            mindim = maxdim = dim(D.inds().front());
            }
        args.add("MinDim",mindim);
        args.add("MaxDim",maxdim);
        }

    //auto ui = commonIndex(AAcomb,Ucomb);
    //auto vi = commonIndex(AAcomb,Vcomb);

    auto spec = svdOrd2(AAcomb,ui,vi,U,D,V,args);

    U = dag(Ucomb) * U;
    V = V * dag(Vcomb);

    return spec;
    } //svd

std::tuple<ITensor,ITensor,ITensor>
svd(ITensor const& AA, IndexSet const& Uis, IndexSet const& Vis,
    Args args)
    {
    if( !hasSameInds(inds(AA),IndexSet(Uis,Vis)) )
      Error("In svd, U indices and V indices must match the indices of the input ITensor");
    return svd(AA,Uis,args);
    }

std::tuple<ITensor,ITensor,ITensor>
svd(ITensor const& AA, IndexSet const& Uis,
    Args args)
    {
    ITensor U(Uis),S,V;
    svd(AA,U,S,V,args);
    auto u = commonIndex(U,S);
    auto v = commonIndex(S,V);
    return std::tuple<ITensor,ITensor,ITensor>(U,S,V);
    }

std::tuple<ITensor,ITensor>
polar(ITensor const& T,
      IndexSet const& Uis)
    {
    auto [U,S,V] = svd(T,Uis);
    auto u = commonIndex(S,U);
    auto v = commonIndex(S,V);
    auto Vis = uniqueInds(inds(V),{v});
    auto Vp = prime(V,Vis)*delta(v,u);
    auto Q = U*Vp;
    auto P = dag(Vp)*S*V;
    return std::tuple<ITensor,ITensor>(Q,P);
    }

// output: truncerr,docut_lower,docut_upper,ndegen_below
std::tuple<Real,Real,Real,int>
truncate(Vector & P,
         long maxdim,
         long mindim,
         Real cutoff,
         bool absoluteCutoff,
         bool doRelCutoff,
         Args const& args)
    {
    auto respectDegenerate = args.getBool("RespectDegenerate",false);

    long origm = P.size();
    long n = origm-1;
    Real docut_lower = 0;
    Real docut_upper = 0;

    //Special case if P's are zero
    if(P(0) == 0.0)
        {
        resize(P,1); 
        auto degen_cutoff = REAL_EPSILON;
        auto docut_lower = -degen_cutoff;
        auto docut_upper = degen_cutoff;
        auto ndegen_below = 1;
        return std::make_tuple(0.,docut_lower,docut_upper,ndegen_below);
        }
    
    if(origm == 1) 
        {
        auto degen_cutoff = 1E-3*P(0);
        auto docut_lower = P(0)-degen_cutoff;
        auto docut_upper = P(0)+degen_cutoff;
        auto ndegen_below = 1;
        return std::make_tuple(0.,docut_lower,docut_upper,ndegen_below);
        }

    // TODO: check this is correct
    // First normalize P by the largest value to scale things properly
    // We will undo this at the end
    auto P0 = P(0);
    P /= P0;
    // If we are not using a relative cutoff,
    // also normalize the cutoff value
    if(absoluteCutoff || !doRelCutoff)
      {
      cutoff /= P0;
      }

    //Zero out any negative weight
    for(auto zn = n; zn >= 0; --zn)
        {
        if(P(zn) >= 0) break;
        P(zn) = 0;
        }

    Real truncerr = 0;
    //Always truncate down to at least m==maxdim (m==n+1)
    while(n >= maxdim)
        {
        truncerr += P(n);
        --n;
        }

    if(absoluteCutoff) //absoluteCutoff is typically false
        {
        //Test if individual prob. weights fall below cutoff
        //rather than using *sum* of discarded weights
        for(; P(n) < cutoff && n >= mindim; --n) 
            {
            truncerr += P(n);
            }
        }
    else
        {
        Real scale = 1.0;
        //if doRelCutoff, use normalized P's when truncating
        if(doRelCutoff) 
            {
            scale = sumels(P);
            if(scale == 0.0) scale = 1.0;
            }

        //Continue truncating until *sum* of discarded probability 
        //weight reaches cutoff reached (or m==mindim)
        while(truncerr+P(n) < cutoff*scale && n >= mindim)
            {
            truncerr += P(n);
            --n;
            }
        truncerr = (scale == 0 ? 0 : truncerr/scale);
        }

    if(n < 0) n = 0;

    //P is 0-indexed, so add 1 to n to 
    //get correct state count m
    auto m = n+1;

    auto degen_cutoff = 1E-3*P(n);
    if(P(n) == 0.0)
        {
        degen_cutoff = REAL_EPSILON;
        }

    docut_lower = P(n)-degen_cutoff;
    docut_upper = P(n)+degen_cutoff;

    //Check for a degeneracy:
    auto ndegen_below = 1;
    if(n > 0)
        {
        while( (n-ndegen_below >= 0) && 
               (std::abs(P(n-ndegen_below)-P(n)) < degen_cutoff) )
            {
            ++ndegen_below;
            }
        }

    auto ndegen_above = 0;
    if(n+1 < origm) 
        {
        while( (n+1+ndegen_above < origm) && 
               (std::abs(P(n+1+ndegen_above)-P(n)) < degen_cutoff) )
            {
            ++ndegen_above;
            }
        }

    // If the new dimension falls within a degenerate subspace and we 
    // are respecting degeneracies (we are making sure to not truncate 
    // within a degenerate subspace)
    if( respectDegenerate && (ndegen_below + ndegen_above > 1) )
        {
        if( maxdim < m+ndegen_above ) // If the maximum dimension falls within the degenerate subspace
            {
            // Truncate the subspace
            m -= ndegen_below;
            ndegen_above += ndegen_below;
            ndegen_below = 0;
            }
        else // Otherwise, keep the subspace
            {
            m += ndegen_above;
            ndegen_below += ndegen_above;
            ndegen_above = 0;
            }
        }

    resize(P,m); 

    // Put back the normalization factor
    P *= P0;
    truncerr *= P0;
    docut_lower *= P0;
    docut_upper *= P0;

    return std::make_tuple(truncerr,docut_lower,docut_upper,ndegen_below);
    } // truncate

void
showEigs(Vector const& P,
         Real truncerr,
         LogNum const& scale,
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

    auto do_truncate = args.getBool("Truncate",true);
    auto cutoff = args.getReal("Cutoff",0.);
    auto maxdim = args.getInt("MaxDim",P.size());
    auto mindim = args.getInt("MinDim",1);
    auto doRelCutoff = args.getBool("DoRelCutoff",true);
    auto absoluteCutoff = args.getBool("AbsoluteCutoff",false);

    println();
    printfln("mindim = %d, maxdim = %d, cutoff = %.2E, truncate = %s",mindim,maxdim,cutoff,do_truncate);
    printfln("Kept m=%d states, trunc. err. = %.3E", P.size(),truncerr);
    printfln("doRelCutoff = %s, absoluteCutoff = %s",doRelCutoff,absoluteCutoff);
    IF_USESCALE(printfln("Scale is = %sexp(%.2f)",scale.sign() > 0 ? "" : "-",scale.logNum());)

    auto stop = std::min(size_t{10},P.size());
    auto Ps = Vector(subVector(P,0,stop));

#ifndef USESCALE
    print("Eigenvalues:");
#else
    if(scale.logNum() < 10 && scale.isFiniteReal())
        {
        Ps *= sqr(scale.real0());
        print("Eigenvalues:");
        }
    else
        {
        print("Eigenvalues [not including scale = ",scale.logNum(),"]:");
        }
#endif

    for(auto n : range(Ps))
        {
        auto eig = Ps(n);
        printf(( eig > 1E-3 && eig < 1000) ? (" %.4f") : (" %.3E") , eig); 
        }
    println();
    } // showEigs

Real
absSqrt(Real x) { return std::sqrt(std::fabs(x)); }

Spectrum
factor(ITensor const& T,
       ITensor      & A,
       ITensor      & B,
       Args args)
    {
    auto itagset = getTagSet(args,"Tags","Link");
    ITensor D;
    auto spec = svd(T,A,D,B,{args,"LeftTags=",itagset});
    auto dl = commonIndex(A,D);
    auto dr = commonIndex(B,D);
    D.apply(absSqrt);
    A *= D;
    B *= D;
    //Replace index dl with dr
    A *= delta(dl,dr);
    return spec;
    }

std::tuple<ITensor,ITensor>
factor(ITensor const& T,
       IndexSet const& Ais,
       IndexSet const& Bis,
       Args const& args)
    {
    if( !hasSameInds(inds(T),IndexSet(Ais,Bis)) )
      Error("In factor, A indices and B indices must match the indices of the input ITensor");
    return factor(T,Ais,args);
    }

std::tuple<ITensor,ITensor>
factor(ITensor const& T,
       IndexSet const& Ais,
       Args const& args)
    {
    ITensor A(Ais),B;
    factor(T,A,B,args);
    return std::tuple<ITensor,ITensor>(A,B);
    }

template<typename value_type>
void
eigDecompImpl(ITensor T, 
              ITensor & L, 
              ITensor & R, 
              ITensor & D,
              Args const& args)
    {
    auto itagset = getTagSet(args,"Tags","Link");
    // New index listing eigenvectors
    Index newmid;
    if(not hasQNs(T))
        {
        auto full = args.getBool("FullDecomp",false);

        if(order(T) != 2)
            {
            Print(order(T));
            Print(T);
            Error("eig_decomp requires 2-index tensor as input");
            }

        auto lind = noPrime(T.inds().front());

        //Do the diagonalization
        auto MM = toMatRefc<value_type>(T,prime(lind),lind);

        Vector Dr, Di;
        Matrix Rr, Ri;
        Matrix Lr, Li;
        if(!full) 
            {
            eigen(MM,Rr,Ri,Dr,Di);
            }
        else
            {
            eigDecomp(MM,Lr,Li,Dr,Di,Rr,Ri);
            }

        newmid = Index(dim(lind),itagset);

        //put right eigenvectors into an ITensor
        if(norm(Ri) > 1E-16*norm(Rr))
            {
            //complex eigenvectors
            auto store = DenseCplx(Rr.size());
            auto ri = Rr.begin();
            auto ii = Ri.begin();
            for(decltype(Rr.size()) n = 0; n < Rr.size(); ++ri, ++ii, ++n)
                {
#ifdef DEBUG
                if(ri == Rr.end() || ii == Ri.end()) Error("out of range iterator");
#endif
                store[n] = Cplx(*ri,*ii);
                }
            R = ITensor({lind,newmid},move(store));
            }
        else
            {
            //real eigenvectors
            R = ITensor({lind,newmid},DenseReal{move(Rr.storage())});
            }

        if(norm(Di) > 1E-16*norm(Dr))
            {
            //complex eigenvalues
            auto store = DiagCplx(Dr.size());
            for(auto n : range(Dr.size()))
                {
                store.store.at(n) = Cplx(Dr(n),Di(n));
                }
            D = ITensor({prime(newmid),newmid},move(store),T.scale());
            }
        else
            {
            //real eigenvectors
            D = ITensor({prime(newmid),newmid},DiagReal{move(Dr.storage())},T.scale());
            }

        if(full)
            {

            // If doing full decomp, prime R
            R.prime();


            //put left eigenvectors into an ITensor
            if(norm(Li) > 1E-16*norm(Lr))
                {
                //complex eigenvectors
                auto store = DenseCplx(Lr.size());
                auto ri = Lr.begin();
                auto ii = Li.begin();
                for(decltype(Lr.size()) n = 0; n < Lr.size(); ++ri, ++ii, ++n)
                    {
#ifdef DEBUG
                    if(ri == Lr.end() || ii == Li.end()) Error("out of range iterator");
#endif
                    store[n] = Cplx(*ri,*ii);
                    }
                L = ITensor({lind,newmid},move(store));
                }
            else
                {
                //real eigenvectors
                L = ITensor({lind,newmid},DenseReal{move(Lr.storage())});
                }
            }
        }
    else
        {
        Error("eigDecompImpl not implemented for QN ITensor");
        }
    }

Spectrum 
denmatDecomp(ITensor const& AA,
             ITensor & A,
             ITensor & B,
             Direction dir,
             Args const& args)
    {
    return denmatDecomp(AA,A,B,dir,NoOp(),args);
    }

std::tuple<ITensor,ITensor>
denmatDecomp(ITensor const& T,
             IndexSet const& Ais,
             IndexSet const& Bis,
             Direction dir,
             Args const& args)
    {
    if( !hasSameInds(inds(T),IndexSet(Ais,Bis)) )
      Error("In svd, A indices and B indices must match the indices of the input ITensor");
    return denmatDecomp(T,Ais,dir,args);
    }

std::tuple<ITensor,ITensor>
denmatDecomp(ITensor const& T,
             IndexSet const& Ais,
             Direction dir,
             Args const& args)
    {
    ITensor A(Ais),B;
    denmatDecomp(T,A,B,dir,args);
    return std::tuple<ITensor,ITensor>(A,B);
    }

Spectrum
diagPosSemiDef(ITensor const& M,
               ITensor      & U,
               ITensor      & D,
               Args args)
    {
    if(!args.defined("Tags")) args.add("Tags","Link");

    //
    // Pick an arbitrary index and do some analysis
    // on its prime level spacing
    //
    auto k = M.inds().front();
    auto kps = stdx::reserve_vector<int>(order(M));
    for(auto& i : M.inds()) if( noPrime(i)==noPrime(k) ) kps.push_back(i.primeLevel());
    if(kps.size() <= 1ul || kps.size()%2 != 0ul)
        {
        Error("Input tensor to diagHermitian should have pairs of indices with equally spaced prime levels");
        }
    auto nk = kps.size();
    std::sort(kps.begin(),kps.end());
    //idiff == "inner" difference between cluster of low-prime-level copies
    //         of k, if more than one
    auto idiff = kps.at(nk/2-1)-kps.front();
    //mdiff == max prime-level difference of copies of k
    auto mdiff = kps.back()-kps.front();
    //pdiff == spacing between lower and higher prime level index pairs
    auto pdiff = mdiff-idiff;

    auto inds = stdx::reserve_vector<Index>(order(M)/2);
    for(auto& i : M.inds())
    for(auto& j : M.inds())
        {
        if( noPrime(i)==noPrime(j) && i.primeLevel()+pdiff == j.primeLevel() )
            {
            inds.push_back(i);
            }
        }
    if(inds.empty() || order(M)/2 != (long)inds.size())
        {
        Error("Input tensor to diagHermitian should have pairs of indices with equally spaced prime levels");
        }

    auto [comb,cind] = combiner(std::move(inds),args);
    (void)cind;
    auto Mc = M*comb;

    auto combP = dag(prime(comb,pdiff));
    try {
        Mc = combP * Mc;
        }
    catch(ITError const& e)
        {
        println("Diagonalize expects opposite arrow directions for primed and unprimed indices.");
        throw e;
        }

    auto spec = diag_hermitian(Mc,U,D,args);

    U = comb * U;

    return spec;
    } //diagHermitian

std::tuple<ITensor,ITensor>
diagPosSemiDef(ITensor const& T,
               Args const& args)
    {
    ITensor U,D;
    diagPosSemiDef(T,U,D,args);
    return std::tuple<ITensor,ITensor>(U,D);
    }

Spectrum
diagHermitian(ITensor const& M,
              ITensor      & U,
              ITensor      & D,
              Args args)
    {
    args.add("Truncate",false);
    return diagPosSemiDef(M,U,D,args);
    }

std::tuple<ITensor,ITensor>
diagHermitian(ITensor const& T,
              Args const& args)
    {
    ITensor U,D;
    diagHermitian(T,U,D,args);
    return std::tuple<ITensor,ITensor>(U,D);
    }

void 
eigen(ITensor const& T, 
      ITensor & V, 
      ITensor & D,
      Args args)
    {
    if(!args.defined("Tags")) args.add("Tags","Link");
    auto colinds = std::vector<Index>{};
    for(auto& I : T.inds())
        { 
        if(I.primeLevel() == 0) colinds.push_back(I);
        }
    auto [comb,cind] = combiner(std::move(colinds),args);

    auto Tc = prime(dag(comb)) * T * comb; 

    // The version where Tc is ordered
    // {cind,prime(cind)} does not work,
    // here we will explicitly permute the
    // ITensor so the indices are ordered
    // {prime(cind),cind}.
    // This is not great since it is an
    // extra reordering of the data,
    // look into fixing eigen properly.
    if(inds(Tc)(1) != prime(cind))
      Tc = permute(Tc,{prime(cind),cind});

    ITensor L;
    Index eigvec_ind;
    if(isComplex(T))
        {
        eigDecompImpl<Cplx>(Tc,L,V,D,args);
        }
    else
        {
        eigDecompImpl<Real>(Tc,L,V,D,args);
        }

    V *= comb;
    }

std::tuple<ITensor,ITensor>
eigen(ITensor const& T,
      Args const& args)
    {
    ITensor V,D;
    eigen(T,V,D,args);
    return std::tuple<ITensor,ITensor>(V,D);
    }

void 
eigDecomp(ITensor const& T, 
          ITensor & R,
          ITensor & D,
          ITensor & Rinv,
          Args const& args)
    {
    auto colinds = std::vector<Index>{};
    for(auto& I : T.inds())
        { 
        if(I.primeLevel() == 0) colinds.push_back(I);
        }
    auto [comb,cind] = combiner(std::move(colinds));
    (void)cind;

    auto Tc = prime(comb) * T * comb; 

    if(isComplex(Tc))
        {
        eigDecompImpl<Cplx>(Tc,Rinv,R,D,{args,"FullDecomp",true});
        }
    else
        {
        eigDecompImpl<Real>(Tc,Rinv,R,D,{args,"FullDecomp",true});
        }

    R = R * prime(comb);
    Rinv = Rinv * comb;
    }


template<typename T>
struct Exp
    {
    T tt = 0.;
    Exp(T t_) : tt(t_) { }

    T
    operator()(Real x) const { return exp(tt*x); }
    };

ITensor
expHermitian(ITensor const& T, Cplx t)
    {
    ITensor U,d;
    diagHermitian(T,U,d);

    if(t.imag()==0.)
        {
        d.apply(Exp<Real>(t.real()));
        }
    else
        {
        d.apply(Exp<Cplx>(t));
        }

    return prime(U)*d*dag(U);
    }

} //namespace itensor
