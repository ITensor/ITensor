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
diagHImpl(ITensor H, 
          ITensor& U, 
          ITensor& D,
          Args const& args)
    {
    auto cutoff = args.getReal("Cutoff",0.);
    auto maxm = args.getInt("Maxm",H.inds().front().m());
    auto minm = args.getInt("Minm",1);
    auto def_do_trunc = args.defined("Cutoff") || args.defined("Maxm");
    auto do_truncate = args.getBool("Truncate",def_do_trunc);
    auto doRelCutoff = args.getBool("DoRelCutoff",true);
    auto absoluteCutoff = args.getBool("AbsoluteCutoff",false);
    auto ignore_degeneracy = args.getBool("IgnoreDegeneracy",true);
    auto showeigs = args.getBool("ShowEigs",false);
    auto iname = args.getString("IndexName","d");

    if(H.r() != 2)
        {
        Print(H.r());
        Print(H);
        Error("Rank greater than 2 in diag_hermitian");
        }

    auto i1 = H.inds().front();
    auto i2 = H.inds().back();

    auto active = (i1.primeLevel() < i2.primeLevel()) ? i1 : i2;

    auto pdiff = std::abs(i1.primeLevel()-i2.primeLevel());


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
    Real docut = -1;
    if(do_truncate)
        {
        //if(DD(1) < 0) DD *= -1; //DEBUG
        tie(truncerr,docut) = truncate(DD,maxm,minm,cutoff,absoluteCutoff,doRelCutoff,args);
        if(ignore_degeneracy)
            {
            m = DD.size();
            }
        else
            {
            long total_m = 0;
            for(decltype(DD.size()) n = 0; n < DD.size() && DD(n) > docut; ++n)
                {
                total_m += 1;
                }
            m = total_m;
            }
#ifdef DEBUG
        if(m==0) throw std::runtime_error("Index of D after diagHermitian is empty. Consider raising Maxm or Cutoff, or making IgnoreDegeneracy true");
#endif
        reduceCols(UU,m);
        }

    if(m > maxm)
        {
        printfln("m > maxm; m = %d, maxm = %d",m,maxm);
        Error("m > maxm");
        }
    if(m > 50000)
        {
        printfln("WARNING: very large m = %d in ITensor diag_hermitian");
        }

    if(showeigs)
        {
        auto showargs = args;
        showargs.add("Cutoff",cutoff);
        showargs.add("Maxm",maxm);
        showargs.add("Minm",minm);
        showargs.add("Truncate",do_truncate);
        showargs.add("DoRelCutoff",doRelCutoff);
        showargs.add("AbsoluteCutoff",absoluteCutoff);
        showEigs(DD,truncerr,sqrt(H.scale()),showargs);
        }

    auto newmid = Index(iname,m,active.type());

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

template<typename T>
Spectrum
diagHImpl(IQTensor    H, 
          IQTensor  & U, 
          IQTensor  & D,
          Args        args)
    {
    SCOPED_TIMER(7)
    auto cutoff = args.getReal("Cutoff",0.);
    auto maxm = args.getInt("Maxm",MAX_INT);
    auto minm = args.getInt("Minm",1);
    auto def_do_trunc = args.defined("Cutoff") || args.defined("Maxm");
    auto do_truncate = args.getBool("Truncate",def_do_trunc);
    auto doRelCutoff = args.getBool("DoRelCutoff",true);
    auto absoluteCutoff = args.getBool("AbsoluteCutoff",false);
    auto ignore_degeneracy = args.getBool("IgnoreDegeneracy",true);
    auto showeigs = args.getBool("ShowEigs",false);
    auto compute_qns = args.getBool("ComputeQNs",false);
    auto iname = args.getString("IndexName","d");
    args.add("IgnoreDegeneracy",ignore_degeneracy);
    auto itype = getIndexType(args,"IndexType",Link);
    if(H.r() != 2)
        {
        Print(H.inds());
        Error("diag_hermitian requires rank 2 input tensor");
        }
    
    auto i1 = H.inds().front();
    auto i2 = H.inds().back();
    auto ai = (i1.primeLevel() < i2.primeLevel()) ? i1 : i2;
    auto pdiff = std::abs(i1.primeLevel()-i2.primeLevel());

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

    auto alleig = stdx::reserve_vector<Real>(ai.m());
    auto alleigqn = vector<EigQN>{};
    if(compute_qns) alleigqn = stdx::reserve_vector<EigQN>(ai.m());

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
            auto q = ai.qn(1+bi);
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
    Real docut = -1;
    if(do_truncate)
        {
        tie(truncerr,docut) = truncate(probs,maxm,minm,cutoff,
                                       absoluteCutoff,doRelCutoff,args);
        m = probs.size();
        alleigqn.resize(m);
        }

    if(showeigs)
        {
        auto showargs = args;
        showargs.add("Cutoff",cutoff);
        showargs.add("Maxm",maxm);
        showargs.add("Minm",minm);
        showargs.add("Truncate",do_truncate);
        showargs.add("DoRelCutoff",doRelCutoff);
        showargs.add("AbsoluteCutoff",absoluteCutoff);
        showEigs(probs,truncerr,H.scale(),showargs);
        }

    if(m > maxm)
        {
        printfln("m > maxm; m = %d, maxm = %d",m,maxm);
        Error("m > maxm");
        }
    if(m > 20000)
        {
        printfln("WARNING: very large m = %d in diag_hermitian",m);
        }

    //3. Truncate eigenvalues and eigenvectors of H

    //Form new Link IQIndex with appropriate m's for each block
    IQIndex::storage iq;
    iq.reserve(Nblock);

    long total_m = 0;
    for(auto b : range(Nblock))
        {
        auto& UU = Umats.at(b);
        auto& d = dvecs.at(b);
        auto& B = blocks[b];

        long this_m = d.size();
        if(do_truncate)
            {
            //Truncate all elems of d falling below docut
            while(this_m > 0 && d(this_m-1) <= docut) --this_m;
            }
        //We need to check that the number of states doesn't
        //go above m, which can happen if there are degeneracies
        total_m += this_m;
        if(total_m > m)
          {
          this_m += (m-total_m);
          total_m = m;
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

        iq.emplace_back(Index(iname+nameint("_",b),this_m, itype),ai.qn(1+B.i1));
        }

    if(iq.empty())
        {
        if(blocks.empty()) Error("No blocks in IQTensor svd");
        auto& B = blocks.front();
        iq.emplace_back(Index(iname+nameint("_",0),1, itype),ai.qn(1+B.i1));
        }

    auto d = IQIndex(iname,move(iq),-ai.dir());

    auto Uis = IQIndexSet(dag(ai),dag(d));
    auto Dis = IQIndexSet(prime(d,pdiff),dag(d));

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

        auto uind = stdx::make_array(B.i1,n);
        auto pU = getBlock(Ustore,Uis,uind);
        assert(pU.data() != nullptr);
        assert(ai[B.i1].m() == long(nrows(UU)));
        auto Uref = makeMatRef(pU,nrows(UU),mm);
        Uref &= UU;

        auto dind = stdx::make_array(n,n);
        auto pD = getBlock(Dstore,Dis,dind);
        assert(pD.data() != nullptr);
        auto Dref = makeVecRef(pD.data(),mm);
        Dref &= dv;

        ++n;
        }

    U = IQTensor(Uis,move(Ustore));
    D = IQTensor(Dis,move(Dstore),H.scale());

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

template<typename I>
Spectrum
diag_hermitian(ITensorT<I>    H, 
               ITensorT<I>  & U, 
               ITensorT<I>  & D,
               Args const& args)
    {
    if(isComplex(H))
        {
        return diagHImpl<Cplx>(H,U,D,args);
        }
    return diagHImpl<Real>(H,U,D,args);
    }
template
Spectrum
diag_hermitian(ITensor    H, 
               ITensor  & U, 
               ITensor  & D,
               Args const& args);
template
Spectrum
diag_hermitian(IQTensor    H, 
               IQTensor  & U, 
               IQTensor  & D,
               Args const& args);

} //namespace itensor
