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

namespace itensor {

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
        return doTask(ToMatRefc<V>{i1.m(),i2.m()},T.store());
        }
    return doTask(ToMatRefc<V>{i2.m(),i1.m(),true},T.store());
    }

/////////////

template<typename T>
struct GetBlocks
    {
    using value_type = T;
    IQIndexSet const& is;
    bool transpose = false;

    GetBlocks(IQIndexSet const& is_, 
              IQIndex const& i1_, 
              IQIndex const& i2_)
      : is(is_)
        { 
        if(is.r() != 2) Error("GetBlocks only supports rank 2 currently");
        transpose = (i2_ == is.front());
        }
    };

template<typename T>
struct Rank2Block
    {
    MatRefc<T> M;
    long i1 = 0,
         i2 = 0;
    };

template<typename T>
vector<Rank2Block<T>>
doTask(GetBlocks<T> const& G, 
       QDense<T> const& d)
    {
    if(G.is.r() != 2) Error("doTask(GetBlocks,QDenseReal) only supports rank 2");
    auto res = vector<Rank2Block<T>>{d.offsets.size()};
    auto dblock = IntArray(2,0);
    size_t n = 0;
    for(auto& dio : d.offsets)
        {
        auto& R = res[n++];
        computeBlockInd(dio.block,G.is,dblock);
        auto nrow = G.is[0][dblock[0]].m();
        auto ncol = G.is[1][dblock[1]].m();
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

///////////////


std::tuple<Real,Real>
truncate(Vector & P,
         long maxm,
         long minm,
         Real cutoff,
         bool absoluteCutoff,
         bool doRelCutoff)
    {
    long origm = P.size();
    long n = origm-1;
    Real docut = 0;
    
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
    //Always truncate down to at least m==maxm
    for(; n >= maxm; --n) truncerr += P(n);

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
        if(doRelCutoff) scale = sumels(P);

        //Continue truncating until *sum* of discarded probability 
        //weight reaches cutoff reached (or m==minm)
        for(;truncerr+P(n) < cutoff*scale && n >= minm; --n)
            {
            truncerr += P(n);
            }
        truncerr = (scale == 0 ? 0 : truncerr/scale);
        }

    if(n < 0) n = 0;

    //P is 0-indexed, so add 1 to n to 
    //get correct state count m
    auto m = n+1;

    if(m < origm) docut = (P(m) + P(m-1))/2. - 1E-5*P(m);

    resize(P,m); 

    return std::make_tuple(truncerr,docut);
    } // truncate

void
showEigs(Vector const& P,
         Real truncerr,
         LogNum const& scale,
         Args const& args)
    {
    auto cutoff = args.getReal("Cutoff",MIN_CUT);
    auto maxm = args.getInt("Maxm",MAX_M);
    auto minm = args.getInt("Minm",1);
    auto do_truncate = args.getBool("Truncate",true);
    auto doRelCutoff = args.getBool("DoRelCutoff",true);
    auto absoluteCutoff = args.getBool("AbsoluteCutoff",false);

    println();
    printfln("minm = %d, maxm = %d, cutoff = %.2E, truncate = %s",minm,maxm,cutoff,do_truncate);
    printfln("Kept m=%d states, trunc. err. = %.3E", P.size(),truncerr);
    printfln("doRelCutoff = %s, absoluteCutoff = %s",doRelCutoff,absoluteCutoff);
    printfln("Scale is = %sexp(%.2f)",scale.sign() > 0 ? "" : "-",scale.logNum());

    auto stop = std::min(size_t{10},P.size());
    auto Ps = Vector(subVector(P,0,stop));

    //Real orderMag = log(std::fabs(P(0))) + scale.logNum();
    if(scale.logNum() < 10 && scale.isFiniteReal())
        {
        Ps *= sqr(scale.real0());
        print("Density matrix evals:");
        }
    else
        {
        print("Density matrix evals [not including scale = ",scale.logNum(),"]:");
        }

    for(auto n : range(Ps))
        {
        auto eig = Ps(n);
        printf(( eig > 1E-3 && eig < 1000) ? (" %.4f") : (" %.3E") , eig); 
        }
    println();
    } // showEigs


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
    auto thresh = args.getReal("SVDThreshold",1E-3);
    auto cutoff = args.getReal("Cutoff",MIN_CUT);
    auto maxm = args.getInt("Maxm",MAX_M);
    auto minm = args.getInt("Minm",1);
    auto do_truncate = args.getBool("Truncate",true);
    auto doRelCutoff = args.getBool("DoRelCutoff",true);
    auto absoluteCutoff = args.getBool("AbsoluteCutoff",false);
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
                                       absoluteCutoff,doRelCutoff);
        m = probs.size();
        resize(DD,m);
        reduceCols(UU,m);
        reduceCols(VV,m);
        }

    Spectrum spec;
    spec.truncerr(truncerr);

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

    if(A.scale().isFiniteReal()) 
        {
        DD *= sqr(A.scale().real0());
        }
    else                         
        {
        println("Warning: scale not finite real after svd");
        }

    spec.eigsKept(move(DD));

    return spec;
    }


template<typename T>
Spectrum
svdImpl(IQTensor A, 
        IQIndex const& uI, 
        IQIndex const& vI,
        IQTensor & U, 
        IQTensor & D, 
        IQTensor & V,
        Args const& args)
    {
    auto thresh = args.getReal("SVDThreshold",1E-4);
    auto cutoff = args.getReal("Cutoff",MIN_CUT);
    auto maxm = args.getInt("Maxm",MAX_M);
    auto minm = args.getInt("Minm",1);
    auto do_truncate = args.getBool("Truncate",true);
    auto doRelCutoff = args.getBool("DoRelCutoff",true);
    auto absoluteCutoff = args.getBool("AbsoluteCutoff",false);
    auto show_eigs = args.getBool("ShowEigs",false);

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
        }

    //Square the singular values into probabilities
    //(density matrix eigenvalues)
    for(auto& sval : alleig) sval = sval*sval;
    //Sort all eigenvalues from largest to smallest
    //irrespective of quantum numbers
    stdx::sort(alleig,std::greater<Real>{});

    auto probs = Vector(move(alleig),VecRange{alleig.size()});

    long m = probs.size();
    Real truncerr = 0;
    Real docut = -1;
    if(do_truncate)
        {
        tie(truncerr,docut) = truncate(probs,maxm,minm,cutoff,
                                       absoluteCutoff,doRelCutoff);
        m = probs.size();
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

    for(auto b : range(Nblock))
        {
        auto& d = dvecs.at(b);
        auto& B = blocks[b];

        //Count number of eigenvalues in the sector above docut
        long this_m = 0;
        for(decltype(d.size()) n = 0; n < d.size() && sqr(d(n)) > docut; ++n)
            {
            this_m += 1;
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

        Liq.emplace_back(Index("l",this_m),uI.qn(1+B.i1));
        Riq.emplace_back(Index("r",this_m),vI.qn(1+B.i2));
        }
    
    auto L = IQIndex("L",move(Liq),uI.dir());
    auto R = IQIndex("R",move(Riq),vI.dir());

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

    return Spectrum(move(probs),Args("Truncerr",truncerr));

    } // svdImpl IQTensor

template<typename IndexT>
Spectrum 
svdRank2(ITensorT<IndexT> const& A, 
         IndexT const& ui, 
         IndexT const& vi,
         ITensorT<IndexT> & U, 
         ITensorT<IndexT> & D, 
         ITensorT<IndexT> & V,
         Args const& args)
    {
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
         ITensor &,ITensor &,ITensor &,Args const&);
template Spectrum 
svdRank2(IQTensor const&,IQIndex const&,IQIndex const&,
         IQTensor &,IQTensor &,IQTensor &,Args const&);


template<typename T>
Spectrum
diagHImpl(ITensor rho, 
          ITensor& U, 
          ITensor& D,
          Args const& args)
    {
    auto cutoff = args.getReal("Cutoff",MIN_CUT);
    auto maxm = args.getInt("Maxm",MAX_M);
    auto minm = args.getInt("Minm",1);
    auto do_truncate = args.getBool("Truncate",false);
    auto doRelCutoff = args.getBool("DoRelCutoff",true);
    auto absoluteCutoff = args.getBool("AbsoluteCutoff",false);
    auto showeigs = args.getBool("ShowEigs",false);

    Index active;
    for(auto& I : rho.inds())
        if(I.primeLevel() == 0)
            {
            active = I;
            break;
            }

    if(!active)
        {
        Print(rho.inds());
        Error("Tensor must have one unprimed index");
        }

    if(rho.r() != 2)
        {
        Print(rho.r());
        Print(rho);
        Error("Rank greater than 2 in diag_hermitian");
        }

    //Depending on the sign of the scale, calling .toMatrix11NoScale 
    //yields a matrix proportional to either rho or -rho.
    //If rho (scale().sign() > 0) then want to temporarily reverse 
    //the sign of the matrix when calling the diagonalization routine
    //to ensure eigenvalues are ordered from largest to smallest.
    if(rho.scale().sign() < 0) rho.scaleTo(rho.scale()*(-1));

    //Do the diagonalization
    Vector DD;
    Mat<T> UU,iUU;
    auto R = toMatRefc<T>(rho,active,prime(active));
    diagHermitian(R,UU,DD);
    conjugate(UU);

    //Truncate
    Real truncerr = 0.0;
    long m = DD.size();
    Real docut = -1;
    if(do_truncate)
        {
        if(DD(1) < 0) DD *= -1; //DEBUG
        tie(truncerr,docut) = truncate(DD,maxm,minm,cutoff,absoluteCutoff,doRelCutoff);
        m = DD.size();
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
        showEigs(DD,truncerr,rho.scale(),showargs);
        }

    auto newmid = Index(active.rawname(),m,active.type());

    U = ITensor({active,newmid},Dense<T>{move(UU.storage())}); 
    D = ITensor({prime(newmid),newmid},DiagReal{DD.begin(),DD.end()},rho.scale());

    if(not rho.scale().isTooBigForReal())
        {
        DD *= rho.scale().real0();
        }
    else
        {
        println("diag_hermitian: scale too big for Real, omitting from returned spectrum.");
        }

    return Spectrum{move(DD),{"Truncerr",truncerr}};
    }

template<typename T>
Spectrum
diagHImpl(IQTensor    rho, 
          IQTensor  & U, 
          IQTensor  & D,
          Args const& args)
    {
    SCOPED_TIMER(7)
    auto cutoff = args.getReal("Cutoff",MIN_CUT);
    auto maxm = args.getInt("Maxm",MAX_M);
    auto minm = args.getInt("Minm",1);
    auto do_truncate = args.getBool("Truncate",false);
    auto doRelCutoff = args.getBool("DoRelCutoff",true);
    auto absoluteCutoff = args.getBool("AbsoluteCutoff",false);
    auto showeigs = args.getBool("ShowEigs",false);

    if(rho.r() != 2)
        {
        Print(rho.inds());
        Error("diag_hermitian requires rank 2 input tensor");
        }
    
    IQIndex ai;
    for(auto& I : rho.inds())
        {
        if(I.primeLevel()==0)
            {
            ai = I;
            break;
            }
        }

    if(not ai) Error("in diag_hermitian rho should have one primed and one unprimed IQIndex");

#ifdef DEBUG
    auto Zero = QN();
    if(div(rho) != Zero)
        { 
        Print(rho); 
        Error("Non-zero divergence of rho, QNs not conserved by Hamiltonian?");
        }
#endif

    if(rho.scale().sign() < 0) rho.scaleTo(rho.scale()*(-1));

    auto blocks = doTask(GetBlocks<T>{rho.inds(),ai,prime(ai)},rho.store());
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

    auto alleig = stdx::reserve_vector<Real>(rho.inds().front().m());

    //1. Diagonalize each ITensor within rho.
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
        totaldsize += rM;
        totalUsize += rM*cM;
        }


    //2. Truncate eigenvalues

    stdx::sort(alleig,std::greater<Real>{});

    auto probs = Vector{move(alleig),VecRange{alleig.size()}};

    //Determine number of states to keep m
    long m = probs.size();
    Real truncerr = 0;
    Real docut = -1;
    if(do_truncate)
        {
        tie(truncerr,docut) = truncate(probs,maxm,minm,cutoff,
                                       absoluteCutoff,doRelCutoff);
        m = probs.size();
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
        showEigs(probs,truncerr,rho.scale(),showargs);
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

    //3. Truncate eigenvalues and eigenvectors of rho

    //Form new Link IQIndex with appropriate m's for each block
    IQIndex::storage iq;
    iq.reserve(Nblock);

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

        if(this_m == 0) 
            { 
            d.clear();
            B.M.clear();
            assert(not B.M);
            continue; 
            }

        d = subVector(d,0,this_m);
        UU = columns(UU,0,this_m);

        iq.emplace_back(Index(nameint("d",b),this_m),ai.qn(1+B.i1));
        }

    if(iq.empty())
        {
        if(blocks.empty()) Error("No blocks in IQTensor svd");
        auto& B = blocks.front();
        iq.emplace_back(Index(nameint("d",0),1),ai.qn(1+B.i1));
        }

    auto d = IQIndex("d",move(iq),-ai.dir());

    auto Uis = IQIndexSet(dag(ai),dag(d));
    auto Dis = IQIndexSet(prime(d),dag(d));

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
    D = IQTensor(Dis,move(Dstore),rho.scale());

    if(rho.scale().isTooBigForReal())
        {
        println("scale too big, omitting from reported eigenvalues");
        }
    else
        {
        probs *= rho.scale().real0();
        }

    return Spectrum{move(probs),{"Truncerr",truncerr}};
    }

template<typename I>
Spectrum
diag_hermitian(ITensorT<I>    rho, 
               ITensorT<I>  & U, 
               ITensorT<I>  & D,
               Args const& args)
    {
    if(isComplex(rho))
        {
        return diagHImpl<Cplx>(rho,U,D,args);
        }
    return diagHImpl<Real>(rho,U,D,args);
    }
template
Spectrum
diag_hermitian(ITensor    rho, 
               ITensor  & U, 
               ITensor  & D,
               Args const& args);
template
Spectrum
diag_hermitian(IQTensor    rho, 
               IQTensor  & U, 
               IQTensor  & D,
               Args const& args);


template<typename Tensor>
void
factor(Tensor const& T,
       Tensor      & A,
       Tensor      & B,
       Args const& args)
    {
    auto name = args.getString("IndexName","c");
    Tensor D;
    svd(T,A,D,B,{args,"LeftIndexName=",name});
    auto dl = commonIndex(A,D);
    auto dr = commonIndex(B,D);
    D.apply([](Real x){ return std::sqrt(std::fabs(x)); });
    A *= D;
    B *= D;
    //Replace index dl with dr
    A *= delta(dl,dr);
    }
template void
factor(ITensor const& T,ITensor& A,ITensor & B,Args const& args);
template void
factor(IQTensor const& T,IQTensor& A,IQTensor & B,Args const& args);

template<typename value_type>
void 
eigDecompImpl(ITensor T, 
              ITensor & L, 
              ITensor & R, 
              ITensor & D,
              Args const& args)
    {
    auto full = args.getBool("FullDecomp",false);

    if(rank(T) != 2)
        {
        Print(rank(T));
        Print(T);
        Error("eig_decomp requires rank 2 tensor as input");
        }

    auto lind = noprime(T.inds().front());

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

    auto newmid = Index("C",lind.m(),lind.type());

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

template<typename value_type>
void 
eigDecompImpl(IQTensor T, 
              IQTensor & L, 
              IQTensor & R, 
              IQTensor & D,
              Args const& args)
    {
    /*
    const bool doRelCutoff = args.getBool("DoRelCutoff",false);
    bool cplx = T.isComplex();

#ifdef DEBUG
    if(T.r() != 2)
        {
        Print(T.r());
        Print(T);
        Error("eig_decomp requires rank 2 tensor as input");
        }
#endif

    const int nblocks = T.blocks().size();

    vector<Matrix> rmatrix(nblocks),
                   imatrix(nblocks);
    vector<Vec> reigs(nblocks),
                   ieigs(nblocks);

    if(T.empty())
        throw ResultIsZero("T has no blocks");

    LogNum refNorm(1);
    if(doRelCutoff)
        {
        Real maxLogNum = -200;
        T.scaleOutNorm();
        for(const ITensor& t : T.blocks())
            {
            maxLogNum = std::max(maxLogNum,t.scale().logNum());
            }
        refNorm = LogNumber(maxLogNum,1);
        }
    T.scaleTo(refNorm);

    //1. Diagonalize each ITensor within rho.
    //   Store results in mmatrix and mvector.
    int itenind = 0;
    for(const ITensor& t : T.blocks())
        {
        Index li = t.indices().front(),
              ri = t.indices().back();

        if(!hasindex(L,li))
            swap(li,ri);

        Matrix &Ur = rmatrix.at(itenind),
               &Ui = imatrix.at(itenind);
        Vec &dr = reigs.at(itenind),
               &di = ieigs.at(itenind);

        //Diag ITensors within rho
        if(!cplx)
            {
            Matrix M;
            t.toMatrix11NoScale(li,ri,M);
            GenEigenValues(M,dr,di,Ur,Ui);
            }
        else
            {
            ITensor ret = realPart(t),
                    imt = imagPart(t);
            ret.scaleTo(refNorm);
            imt.scaleTo(refNorm);
            Matrix Mr,Mi;
            ret.toMatrix11NoScale(li,ri,Mr);
            imt.toMatrix11NoScale(li,ri,Mi);
            ComplexEigenvalues(Mr,Mi,dr,di,Ur,Ui);
            }

        ++itenind;
        }


    //Build blocks for unitary diagonalizing rho
    vector<ITensor> Vblocks,
                    Dblocks;

    //Also form new Link IQIndex with appropriate m's for each block
    IQIndex::Storage iq;
    iq.reserve(T.blocks().size());

    itenind = 0;
    for(const ITensor& t : T.blocks())
        {
        Vec &dr = reigs.at(itenind),
               &di = ieigs.at(itenind);
        Matrix &Ur = rmatrix.at(itenind),
               &Ui = imatrix.at(itenind);

        Index nm("d",dr.Length());

        Index act = t.indices().front();
        if(!hasindex(R,act))
            act = t.indices().back();

        iq.push_back(IndexQN(nm,qn(R,act)));

        ITensor blk(act,nm,Ur);
        if(Norm(Ui.TreatAsVector()) > 1E-12)
            {
            blk += Complex_i*ITensor(act,nm,Ui);
            }
        Vblocks.push_back(blk);

        ITensor Dblk(prime(nm),nm,dr);
        if(Norm(di) > 1E-12)
            {
            Dblk += Complex_i*ITensor(prime(nm),nm,di);
            }
        Dblocks.push_back(Dblk);

        ++itenind;
        }

    if(iq.size() == 0)
        {
        throw ResultIsZero("iq.size() == 0");
        }

    IQIndex newmid("L",iq,-R.dir());

    V = IQTensor(dag(R),dag(newmid));
    for(const ITensor& t : Vblocks)
        {
        V += t;
        }

    D = IQTensor(prime(newmid),dag(newmid));
    for(const ITensor& t : Dblocks)
        {
        D += t;
        }

    D *= refNorm;

    */
    }

template<typename index_type>
void 
eigen(ITensorT<index_type> const& T, 
      ITensorT<index_type> & V, 
      ITensorT<index_type> & D,
      Args const& args)
    {
    auto colinds = std::vector<index_type>{};
    for(auto& I : T.inds())
        { 
        if(I.primeLevel() == 0) colinds.push_back(I);
        }
    auto comb = combiner(std::move(colinds));

    auto Tc = prime(comb) * T * comb; 

    ITensorT<index_type> L;
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
template void 
eigen(ITensor const&, ITensor&, ITensor&, Args const&);
template void 
eigen(IQTensor const&, IQTensor&,IQTensor&, Args const&);

template<typename index_type>
void 
eigDecomp(ITensorT<index_type> const& T, 
          ITensorT<index_type> & R,
          ITensorT<index_type> & D,
          ITensorT<index_type> & Rinv,
          Args const& args)
    {
    auto colinds = std::vector<index_type>{};
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
template void 
eigDecomp(ITensor const&, ITensor &, ITensor & , ITensor & , Args const& );
template void 
eigDecomp(IQTensor const&, IQTensor &,IQTensor & , IQTensor & , Args const& );


template<typename I>
ITensorT<I>
expHermitian(ITensorT<I> const& T)
    {
    ITensorT<I> U;
    ITensorT<I> d;
    diagHermitian(T,U,d);

    struct Exp
        {
        Real
        operator()(Real x) const { return exp(x); }
        Cplx
        operator()(Cplx z) const { return exp(z); }
        };
    d.apply(Exp());

    return prime(U)*d*dag(U);
    }
template ITensor expHermitian(ITensor const& T);
template IQTensor expHermitian(IQTensor const& T);

} //namespace itensor
