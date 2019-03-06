//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include <algorithm>
#include <tuple>
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
    if(G.is.r() != 2) Error("doTask(GetBlocks,QDenseReal) only supports 2-index tensors");
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


std::tuple<Real,Real>
truncate(Vector & P,
         long maxm,
         long minm,
         Real cutoff,
         bool absoluteCutoff,
         bool doRelCutoff,
         Args const& args)
    {
    long origm = P.size();
    long n = origm-1;
    Real docut = 0;

    //Special case if P's are zero
    if(P(0) == 0.0)
        {
        resize(P,1); 
        return std::make_tuple(0.,0.);
        }
    
    if(origm == 1) 
        {
        docut = P(0)/2.;
        return std::make_tuple(0,0);
        }

    //Zero out any negative weight
    for(auto zn = n; zn >= 0; --zn)
        {
        if(P(zn) >= 0) break;
        P(zn) = 0;
        }

    Real truncerr = 0;
    //Always truncate down to at least m==maxm (m==n+1)
    while(n >= maxm)
        {
        truncerr += P(n);
        --n;
        }

    if(absoluteCutoff) //absoluteCutoff is typically false
        {
        //Test if individual prob. weights fall below cutoff
        //rather than using *sum* of discarded weights
        for(; P(n) < cutoff && n >= minm; --n) 
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
        //weight reaches cutoff reached (or m==minm)
        while(truncerr+P(n) < cutoff*scale && n >= minm)
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

    if(n+1 < origm) 
        {
        docut = (P(n+1) + P(n))/2.;
        //Check for a degeneracy:
        if(std::fabs(P(n+1)-P(n)) < 1E-3*P(n)) 
            {
            docut += 1E-3*P(n);
            }
        }

    resize(P,m); 

    return std::make_tuple(truncerr,docut);
    } // truncate

void
showEigs(Vector const& P,
         Real truncerr,
         LogNum const& scale,
         Args const& args)
    {
    auto do_truncate = args.getBool("Truncate",true);
    auto cutoff = args.getReal("Cutoff",0.);
    auto maxm = args.getInt("Maxm",P.size());
    auto minm = args.getInt("Minm",1);
    auto doRelCutoff = args.getBool("DoRelCutoff",true);
    auto absoluteCutoff = args.getBool("AbsoluteCutoff",false);

    println();
    printfln("minm = %d, maxm = %d, cutoff = %.2E, truncate = %s",minm,maxm,cutoff,do_truncate);
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



Spectrum
factor(ITensor const& T,
       ITensor      & A,
       ITensor      & B,
       Args const& args)
    {
    //TODO: make a standard TagSet for factor()
    //auto name = args.getString("IndexName","c");
    auto itagset = getTagSet(args,"Tags","Link,FAC");
    ITensor D;
    auto spec = svd(T,A,D,B,{args,"LeftTags=",toString(itagset)});
    auto dl = commonIndex(A,D);
    auto dr = commonIndex(B,D);
    D.apply([](Real x){ return std::sqrt(std::fabs(x)); });
    A *= D;
    B *= D;
    //Replace index dl with dr
    A *= delta(dl,dr);
    return spec;
    }

//TODO: create a tag convention
template<typename value_type>
void 
eigDecompImpl(ITensor T, 
              ITensor & L, 
              ITensor & R, 
              ITensor & D,
              Args const& args)
    {
    if(not hasQNs(T))
        {
        auto full = args.getBool("FullDecomp",false);

        if(ord(T) != 2)
            {
            Print(ord(T));
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

        auto newmid = Index(dim(lind));

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
//
//template<typename value_type>
//void 
//eigDecompImpl(IQTensor T, 
//              IQTensor & L, 
//              IQTensor & R, 
//              IQTensor & D,
//              Args const& args)
//    {
//    /*
//    const bool doRelCutoff = args.getBool("DoRelCutoff",false);
//    bool cplx = T.isComplex();
//
//#ifdef DEBUG
//    if(T.r() != 2)
//        {
//        Print(T.r());
//        Print(T);
//        Error("eig_decomp requires 2-index tensor as input");
//        }
//#endif
//
//    const int nblocks = T.blocks().size();
//
//    vector<Matrix> rmatrix(nblocks),
//                   imatrix(nblocks);
//    vector<Vec> reigs(nblocks),
//                   ieigs(nblocks);
//
//    if(T.empty())
//        throw ResultIsZero("T has no blocks");
//
//    LogNum refNorm(1);
//    if(doRelCutoff)
//        {
//        Real maxLogNum = -200;
//        T.scaleOutNorm();
//        for(const ITensor& t : T.blocks())
//            {
//            maxLogNum = std::max(maxLogNum,t.scale().logNum());
//            }
//        refNorm = LogNumber(maxLogNum,1);
//        }
//    T.scaleTo(refNorm);
//
//    //1. Diagonalize each ITensor within rho.
//    //   Store results in mmatrix and mvector.
//    int itenind = 0;
//    for(const ITensor& t : T.blocks())
//        {
//        Index li = t.indices().front(),
//              ri = t.indices().back();
//
//        if(!hasIndex(L,li))
//            swap(li,ri);
//
//        Matrix &Ur = rmatrix.at(itenind),
//               &Ui = imatrix.at(itenind);
//        Vec &dr = reigs.at(itenind),
//               &di = ieigs.at(itenind);
//
//        //Diag ITensors within rho
//        if(!cplx)
//            {
//            Matrix M;
//            t.toMatrix11NoScale(li,ri,M);
//            GenEigenValues(M,dr,di,Ur,Ui);
//            }
//        else
//            {
//            ITensor ret = realPart(t),
//                    imt = imagPart(t);
//            ret.scaleTo(refNorm);
//            imt.scaleTo(refNorm);
//            Matrix Mr,Mi;
//            ret.toMatrix11NoScale(li,ri,Mr);
//            imt.toMatrix11NoScale(li,ri,Mi);
//            ComplexEigenvalues(Mr,Mi,dr,di,Ur,Ui);
//            }
//
//        ++itenind;
//        }
//
//
//    //Build blocks for unitary diagonalizing rho
//    vector<ITensor> Vblocks,
//                    Dblocks;
//
//    //Also form new Link IQIndex with appropriate m's for each block
//    IQIndex::Storage iq;
//    iq.reserve(T.blocks().size());
//
//    itenind = 0;
//    for(const ITensor& t : T.blocks())
//        {
//        Vec &dr = reigs.at(itenind),
//               &di = ieigs.at(itenind);
//        Matrix &Ur = rmatrix.at(itenind),
//               &Ui = imatrix.at(itenind);
//
//        Index nm("d",dr.Length());
//
//        Index act = t.indices().front();
//        if(!hasIndex(R,act))
//            act = t.indices().back();
//
//        iq.push_back(IndexQN(nm,qn(R,act)));
//
//        ITensor blk(act,nm,Ur);
//        if(Norm(Ui.TreatAsVector()) > 1E-12)
//            {
//            blk += Complex_i*ITensor(act,nm,Ui);
//            }
//        Vblocks.push_back(blk);
//
//        ITensor Dblk(prime(nm),nm,dr);
//        if(Norm(di) > 1E-12)
//            {
//            Dblk += Complex_i*ITensor(prime(nm),nm,di);
//            }
//        Dblocks.push_back(Dblk);
//
//        ++itenind;
//        }
//
//    if(iq.size() == 0)
//        {
//        throw ResultIsZero("iq.size() == 0");
//        }
//
//    IQIndex newmid("L",iq,-R.dir());
//
//    V = IQTensor(dag(R),dag(newmid));
//    for(const ITensor& t : Vblocks)
//        {
//        V += t;
//        }
//
//    D = IQTensor(prime(newmid),dag(newmid));
//    for(const ITensor& t : Dblocks)
//        {
//        D += t;
//        }
//
//    D *= refNorm;
//
//    */
//    }

void 
eigen(ITensor const& T, 
      ITensor & V, 
      ITensor & D,
      Args const& args)
    {
    auto colinds = std::vector<Index>{};
    for(auto& I : T.inds())
        { 
        if(I.primeLevel() == 0) colinds.push_back(I);
        }
    auto comb = combiner(std::move(colinds));

    auto Tc = prime(comb) * T * comb; 

    ITensor L;
    if(isComplex(T))
        {
        eigDecompImpl<Cplx>(Tc,L,V,D,args);
        }
    else
        {
        eigDecompImpl<Real>(Tc,L,V,D,args);
        }

    V = V * comb;
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
    auto comb = combiner(std::move(colinds));

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
