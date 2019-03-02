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

const auto MAX_INT = std::numeric_limits<int>::max();

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
        Index const& ui, 
        Index const& vi,
        ITensor & U, 
        ITensor & D, 
        ITensor & V,
        Args const& args)
    {
    SCOPED_TIMER(7);
    auto do_truncate = args.getBool("Truncate");
    auto thresh = args.getReal("SVDThreshold",1E-3);
    auto cutoff = args.getReal("Cutoff",MIN_CUT);
    auto maxm = args.getInt("Maxm",MAX_M);
    auto minm = args.getInt("Minm",1);
    auto doRelCutoff = args.getBool("DoRelCutoff",true);
    auto absoluteCutoff = args.getBool("AbsoluteCutoff",false);
    auto ignore_degeneracy = args.getBool("IgnoreDegeneracy",true);
    auto lname = args.getString("LeftIndexName","ul");
    auto rname = args.getString("RightIndexName","vl");
    auto itype = getIndexType(args,"IndexType",Link);
    auto litype = getIndexType(args,"LeftIndexType",itype);
    auto ritype = getIndexType(args,"RightIndexType",itype);
    auto show_eigs = args.getBool("ShowEigs",false);

    auto M = toMatRefc<T>(A,ui,vi);

    Mat<T> UU,VV;
    Vector DD;

    TIMER_START(6)
    SVD(M,UU,DD,VV,thresh);
    TIMER_STOP(6)

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
    Real docut = -1;
    long m = DD.size();
    if(do_truncate)
        {
        tie(truncerr,docut) = truncate(probs,maxm,minm,cutoff,
                                       absoluteCutoff,doRelCutoff,args);
        if(ignore_degeneracy)
            {
            m = probs.size();
            }
        else
            {
            long total_m = 0;
            for(decltype(probs.size()) n = 0; n < probs.size() && probs(n) > docut; ++n)
                {
                total_m += 1;
                }
            m = total_m;
            }

#ifdef DEBUG
        if(m==0) throw std::runtime_error("Index of S after SVD is empty. Consider raising Maxm or Cutoff, or making IgnoreDegeneracy true");
#endif
        resize(DD,m);
        reduceCols(UU,m);
        reduceCols(VV,m);
        }


    if(show_eigs) 
        {
        auto showargs = args;
        showargs.add("Cutoff",cutoff);
        showargs.add("Maxm",maxm);
        showargs.add("Minm",minm);
        showargs.add("Truncate",do_truncate);
        showargs.add("DoRelCutoff",doRelCutoff);
        showargs.add("AbsoluteCutoff",absoluteCutoff);
        showEigs(probs,truncerr,A.scale(),showargs);
        }
    
    Index uL(lname,m,litype),
          vL(rname,m,ritype);

    //Fix sign to make sure D has positive elements
    Real signfix = (A.scale().sign() == -1) ? -1 : +1;
    D = ITensor({uL,vL},
                Diag<Real>{DD.begin(),DD.end()},
                A.scale()*signfix);
    U = ITensor({ui,uL},Dense<T>(move(UU.storage())),LogNum(signfix));
    V = ITensor({vi,vL},Dense<T>(move(VV.storage())));

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



template<typename T>
Spectrum
svdImpl(IQTensor A, 
        IQIndex const& uI, 
        IQIndex const& vI,
        IQTensor & U, 
        IQTensor & D, 
        IQTensor & V,
        Args       args)
    {
    auto do_truncate = args.getBool("Truncate");
    auto thresh = args.getReal("SVDThreshold",1E-3);
    auto cutoff = args.getReal("Cutoff",0);
    auto maxm = args.getInt("Maxm",MAX_INT);
    auto minm = args.getInt("Minm",1);
    auto doRelCutoff = args.getBool("DoRelCutoff",true);
    auto absoluteCutoff = args.getBool("AbsoluteCutoff",false);
    auto ignore_degeneracy = args.getBool("IgnoreDegeneracy",true);
    auto show_eigs = args.getBool("ShowEigs",false);
    auto lname = args.getString("LeftIndexName","ul");
    auto rname = args.getString("RightIndexName","vl");
    auto itype = getIndexType(args,"IndexType",Link);
    auto litype = getIndexType(args,"LeftIndexType",itype);
    auto ritype = getIndexType(args,"RightIndexType",itype);
    auto compute_qn = args.getBool("ComputeQNs",false);

    args.add("IgnoreDegeneracy",ignore_degeneracy);

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

    auto alleig = stdx::reserve_vector<Real>(std::min(uI.m(),vI.m()));

    auto alleigqn = vector<EigQN>{};
    if(compute_qn)
        {
        alleigqn = stdx::reserve_vector<EigQN>(std::min(uI.m(),vI.m()));
        }

    if(uI.m() == 0) throw ResultIsZero("uI.m() == 0");
    if(vI.m() == 0) throw ResultIsZero("vI.m() == 0");

    for(auto b : range(Nblock))
        {
        auto& M = blocks[b].M;
        auto& UU = Umats.at(b);
        auto& VV = Vmats.at(b);
        auto& d =  dvecs.at(b);

        SVD(M,UU,d,VV,thresh);

        //conjugate VV so later we can just do
        //U*D*V to reconstruct ITensor A:
        conjugate(VV);

        alleig.insert(alleig.end(),d.begin(),d.end());
        if(compute_qn)
            {
            auto bi = blocks[b].i1;
            auto q = uI.qn(1+bi);
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
    Real docut = -1;
    if(do_truncate)
        {
        tie(truncerr,docut) = truncate(probs,maxm,minm,cutoff,
                                       absoluteCutoff,doRelCutoff,args);
        m = probs.size();
        alleigqn.resize(m);
        }

    if(show_eigs) 
        {
        auto showargs = args;
        showargs.add("Cutoff",cutoff);
        showargs.add("Maxm",maxm);
        showargs.add("Minm",minm);
        showargs.add("Truncate",do_truncate);
        showargs.add("DoRelCutoff",doRelCutoff);
        showargs.add("AbsoluteCutoff",absoluteCutoff);
        showEigs(probs,truncerr,A.scale(),showargs);
        }

    auto Liq = IQIndex::storage{};
    auto Riq = IQIndex::storage{};
    Liq.reserve(Nblock);
    Riq.reserve(Nblock);

    long total_m = 0;
    for(auto b : range(Nblock))
        {
        auto& d = dvecs.at(b);
        auto& B = blocks[b];

        //Count number of eigenvalues in the sector above docut
        long this_m = 0;
        for(decltype(d.size()) n = 0; n < d.size() && sqr(d(n)) > docut; ++n)
            {
            //We need to check that the number of states doesn't
            //go above m, which can happen if there are degeneracies
            if(m > total_m)
                {
                total_m += 1;
                this_m += 1;
                }
            if(d(n) < 0) d(n) = 0;
            }

        if(m == 0 && d.size() >= 1) // zero mps, just keep one arb state
            { 
            this_m = 1; 
            m = 1; 
            docut = 1; 
            }

        if(this_m == 0) 
            { 
            d.clear();
            B.M.clear();
            assert(not B.M);
            continue; 
            }

        resize(d,this_m);

        Liq.emplace_back(Index("l",this_m,litype),uI.qn(1+B.i1));
        Riq.emplace_back(Index("r",this_m,ritype),vI.qn(1+B.i2));
        }

#ifdef DEBUG
    if(Liq.empty() || Riq.empty()) throw std::runtime_error("IQIndex of S after SVD is empty. Consider raising Maxm or Cutoff, or making IgnoreDegeneracy true");
#endif

    auto L = IQIndex(lname,move(Liq),uI.dir());
    auto R = IQIndex(rname,move(Riq),vI.dir());

    auto Uis = IQIndexSet(uI,dag(L));
    auto Dis = IQIndexSet(L,R);
    auto Vis = IQIndexSet(vI,dag(R));

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
        //Print(uI[B.i1].m());
        //Print(L[n].m());

        auto uind = stdx::make_array(B.i1,n);
        auto pU = getBlock(Ustore,Uis,uind);
        assert(pU.data() != nullptr);
        assert(uI[B.i1].m() == long(nrows(UU)));
        auto Uref = makeMatRef(pU,uI[B.i1].m(),L[n].m());
        reduceCols(UU,L[n].m());
        Uref &= UU;

        auto dind = stdx::make_array(n,n);
        auto pD = getBlock(Dstore,Dis,dind);
        assert(pD.data() != nullptr);
        auto Dref = makeVecRef(pD.data(),d.size());
        Dref &= d;

        auto vind = stdx::make_array(B.i2,n);
        auto pV = getBlock(Vstore,Vis,vind);
        assert(pV.data() != nullptr);
        assert(vI[B.i2].m() == long(nrows(VV)));
        auto Vref = makeMatRef(pV.data(),pV.size(),vI[B.i2].m(),R[n].m());
        reduceCols(VV,R[n].m());
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
    U = IQTensor(Uis,move(Ustore));
    D = IQTensor(Dis,move(Dstore),A.scale()*signfix);
    V = IQTensor(Vis,move(Vstore),LogNum{signfix});
    
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

    } // svdImpl IQTensor

template<typename IndexT>
Spectrum 
svdRank2(ITensorT<IndexT> const& A, 
         IndexT const& ui, 
         IndexT const& vi,
         ITensorT<IndexT> & U, 
         ITensorT<IndexT> & D, 
         ITensorT<IndexT> & V,
         Args args)
    {
    auto do_truncate = args.defined("Cutoff") 
                    || args.defined("Maxm");
    if(not args.defined("Truncate")) 
        {
        args.add("Truncate",do_truncate);
        }

    if(A.r() != 2) 
        {
        Print(A);
        Error("A must be matrix-like (rank 2)");
        }
    if(isComplex(A))
        {
        return svdImpl<Cplx>(A,ui,vi,U,D,V,args);
        }
    return svdImpl<Real>(A,ui,vi,U,D,V,args);
    }
template Spectrum 
svdRank2(ITensor const&,Index const&,Index const&,
         ITensor &,ITensor &,ITensor &,Args );
template Spectrum 
svdRank2(IQTensor const&,IQIndex const&,IQIndex const&,
         IQTensor &,IQTensor &,IQTensor &,Args );

} //namespace itensor
