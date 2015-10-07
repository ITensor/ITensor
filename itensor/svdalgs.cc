//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include <algorithm>
#include <tuple>
#include "itensor/util/stdx.h"
#include "itensor/tensor/algs.h"
#include "itensor/svdalgs.h"

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

struct ToMatrixRefc
    {
    long nrows=0,
         ncols=0;
    bool transpose=false;
    ToMatrixRefc(long nr, long nc, bool trans=false) 
        : nrows(nr), ncols(nc), transpose(trans)
        { }
    };

MatrixRefc
doTask(ToMatrixRefc const& T, 
       ITReal const& d)
    {
    auto res = makeMatRef(d.data(),d.size(),T.nrows,T.ncols);
    if(T.transpose) return transpose(res);
    return res;
    }


MatrixRefc
toMatrixRefc(ITensor const& T, 
          Index const& i1, 
          Index const& i2)
    {
    if(i1 == T.inds().front())
        {
        return doTask(ToMatrixRefc{i1.m(),i2.m()},T.store());
        }
    return doTask(ToMatrixRefc{i2.m(),i1.m(),true},T.store());
    }

/////////////

struct GetBlocks
    {
    IQIndexSet const& is;
    bool transpose = false;

    GetBlocks(IQIndexSet const& is_, 
              IQIndex const& i1_, 
              IQIndex const& i2_)
      : is(is_)
        { 
#ifdef DEBUG
        if(is.r() != 2) Error("GetBlocks only supports rank 2 currently");
#endif
        transpose = (i2_ == is.front());
        }
    };

struct Rank2Block
    {
    MatrixRefc M;
    long i1 = 0,
         i2 = 0;
    };

vector<Rank2Block>
doTask(GetBlocks const& G, 
       IQTReal const& d)
    {
#ifdef DEBUG
    if(G.is.r() != 2) Error("doTask(GetBlocks,IQTReal) only supports rank 2");
#endif
    vector<Rank2Block> res{d.offsets.size()};
    Label dblock(2,0);
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
        docut = P(1)/2.;
        return 0;
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
        //Truncate all probability weights below cutoff (or m==minm)
        for(; P(n) < cutoff && n >= minm; --n) truncerr += P(n);
        }
    else
        {
        Real scale = doRelCutoff ? P(1) : 1.0;
        //Continue truncating until *sum* of discarded probability 
        //weight reaches cutoff reached (or m==minm)
        for(;truncerr+P(n) < cutoff*scale && n >= minm; --n)
            {
            truncerr += P(n);
            }
        truncerr = (P(1) == 0 ? 0 : truncerr/scale);
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
    auto doRelCutoff = args.getBool("DoRelCutoff",false);
    auto absoluteCutoff = args.getBool("AbsoluteCutoff",false);

    println();
    printfln("minm = %d, maxm = %d, cutoff = %.2E, truncate = %s",minm,maxm,cutoff,do_truncate);
    printfln("Kept m=%d states, trunc. err. = %.3E", P.size(),truncerr);
    printfln("doRelCutoff = %s, absoluteCutoff = %s",doRelCutoff,absoluteCutoff);
    printfln("Scale is = %sexp(%.2f)",scale.sign() > 0 ? "" : "-",scale.logNum());

    auto stop = std::min(10ul,P.size());
    auto Ps = Vector(subVector(P,0,stop));

    Real orderMag = log(std::fabs(P(1))) + scale.logNum();
    if(std::fabs(orderMag) < 5 && scale.isFiniteReal())
        {
        Ps *= sqr(scale.real0());
        print("Denmat evals: ");
        }
    else
        {
        printf("Denmat evals (not including log(scale) = %.2f): ",scale.logNum());
        }

    for(decltype(stop) j = 0; j < stop; ++j)
        {
        auto eig = Ps(j);
        printf(( eig > 1E-3 && eig < 1000) ? ("%.4f") : ("%.3E") , eig); 
        print((j != stop) ? ", " : "\n");
        }
    println();
    } // showEigs

Spectrum 
svdRank2(ITensor A, 
         Index const& ui, 
         Index const& vi,
         ITensor & U, 
         ITensor & D, 
         ITensor & V,
         Args const& args)
    {
    auto thresh = args.getReal("SVDThreshold",1E-3);
    auto cutoff = args.getReal("Cutoff",MIN_CUT);
    auto maxm = args.getInt("Maxm",MAX_M);
    auto minm = args.getInt("Minm",1);
    auto do_truncate = args.getBool("Truncate",true);
    auto doRelCutoff = args.getBool("DoRelCutoff",true);
    auto absoluteCutoff = args.getBool("AbsoluteCutoff",false);
    auto cplx = isComplex(A);
    auto lname = args.getString("LeftIndexName","ul");
    auto rname = args.getString("RightIndexName","vl");
    auto itype = getIndexType(args,"IndexType",Link);
    auto litype = getIndexType(args,"LeftIndexType",itype);
    auto ritype = getIndexType(args,"RightIndexType",itype);
    auto show_eigs = args.getBool("ShowEigs",false);

#ifdef DEBUG
    if(A.r() != 2) Error("A must be matrix-like (rank 2)");
#endif

    Matrix UU,VV,
           iUU,iVV;
    Vector DD;

    if(!cplx)
        {
        auto M = toMatrixRefc(A,ui,vi);
        SCOPED_TIMER(6)
        SVD(M,UU,DD,VV,thresh);
        }
    else
        {
        Error("Complex ITensor SVD not yet implemented");
        ITensor Are = realPart(A),
                Aim = imagPart(A);
        Are.scaleTo(A.scale());
        Aim.scaleTo(A.scale());
        //Matrix Mre,Mim;
        //Are.toMatrix11NoScale(ui,vi,Mre);
        //Aim.toMatrix11NoScale(ui,vi,Mim);

        //SVD(Mre,Mim,UU,iUU,DD,VV,iVV,thresh);
        }

    //
    // Truncate
    //


    Vector probs;
    if(do_truncate || show_eigs)
        {
        probs = DD;
        for(auto j : index(probs)) probs(j) = sqr(probs(j));
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

    if(cplx)
        {
        Error("Complex ITensor SVD not yet implemented (2)");
        //ITensor iU(ui,uL,iUU.Columns(1,m)),
        //        iV(vL,vi,iVV.Rows(1,m));
        //if(iU.norm() > 1E-14)
        //    U = U + iU*Complex_i;
        //if(iV.norm() > 1E-14)
        //    V = V + iV*Complex_i;
        }
    else
        {
        D = ITensor({uL,vL},
                    ITDiag<Real>{DD.begin(),DD.end()},
                    A.scale()*signfix);
        U = ITensor({ui,uL},ITReal(move(UU.storage())),LogNum(signfix));
        V = ITensor({vi,vL},ITReal(move(VV.storage())));
        }


    //Square all singular values
    //since convention is to report
    //density matrix eigs

    for(auto& el : DD) el = sqr(el);

    if(A.scale().isFiniteReal()) DD *= sqr(A.scale().real0());
    else                         println("Warning: scale not finite real");

    spec.eigsKept(move(DD));

    return spec;

    } // svdRank2 ITensor

Spectrum
svdRank2(IQTensor A, 
         IQIndex const& uI, 
         IQIndex const& vI,
         IQTensor & U, 
         IQTensor & D, 
         IQTensor & V,
         Args const& args)
    {
    auto cplx = isComplex(A);
    auto thresh = args.getReal("SVDThreshold",1E-4);
    auto cutoff = args.getReal("Cutoff",MIN_CUT);
    auto maxm = args.getInt("Maxm",MAX_M);
    auto minm = args.getInt("Minm",1);
    auto do_truncate = args.getBool("Truncate",true);
    auto doRelCutoff = args.getBool("DoRelCutoff",false);
    auto absoluteCutoff = args.getBool("AbsoluteCutoff",false);
    auto show_eigs = args.getBool("ShowEigs",false);

    if(A.r() != 2) Error("A must be matrix-like");

    auto blocks = doTask(GetBlocks{A.inds(),uI,vI},A.store());

    auto Nblock = blocks.size();
    if(Nblock == 0) throw ResultIsZero("IQTensor has no blocks");

    //TODO: optimize allocation/lookup of Umats,Vmats
    //      etc. by allocating memory ahead of time (see algs.cc)
    //      and making Umats a vector of MatrixRef's to this memory
    vector<Matrix> Umats(Nblock),
                   Vmats(Nblock),
                   iUmats,
                   iVmats;
    if(cplx)
        {
        iUmats.resize(Nblock);
        iVmats.resize(Nblock);
        }

    //TODO: allocate dvecs in a single allocation
    //      make dvecs a vector<VecRef>
    auto dvecs = vector<Vector>(Nblock);

    auto alleig = stdx::reserve_vector<Real>(std::min(uI.m(),vI.m()));

    if(uI.m() == 0) throw ResultIsZero("uI.m() == 0");
    if(vI.m() == 0) throw ResultIsZero("vI.m() == 0");

    for(decltype(Nblock) b = 0; b < Nblock; ++b)
        {
        auto& M = blocks[b].M;
        auto& UU = Umats.at(b);
        auto& VV = Vmats.at(b);
        auto& d =  dvecs.at(b);

        //printfln("Block %d = \n%s",b,M);

        if(!cplx)
            {
            SVD(M,UU,d,VV,thresh);
            }
        else
            {
            Error("Complex IQTensor SVD not yet implemented");

            //ITensor ret = realPart(t),
            //        imt = imagPart(t);
            //Matrix Mre(ui->m(),vi->m()),
            //       Mim(ui->m(),vi->m());
            //ret.toMatrix11NoScale(*ui,*vi,Mre);
            //imt.toMatrix11NoScale(*ui,*vi,Mim);

            ////SVDComplex(Mre,Mim,
            ////           UU,iUmatrix.at(itenind),
            ////           d,
            ////           VV,iVmatrix.at(itenind));
            //SVD(Mre,Mim,
            //    UU,iUmatrix.at(itenind),
            //    d,
            //    VV,iVmatrix.at(itenind),
            //    thresh);
            }

        alleig.insert(alleig.end(),d.begin(),d.end());
        }

    //Square the singular values into probabilities
    //(density matrix eigenvalues)
    for(auto& sval : alleig) sval = sval*sval;
    //Sort all eigenvalues from largest to smallest
    //irrespective of quantum numbers
    stdx::sort(alleig,std::greater<Real>{});

    //print("alleig = "); for(auto& el : alleig) print(" ",el); println();

    auto pstore = alleig;
    auto probs = Vector(move(pstore),VecRange{pstore.size()});

    //print("probs = "); for(auto& el : probs) print(" ",el); println();

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

    IQIndex::storage Liq,
                     Riq;
    Liq.reserve(Nblock);
    Riq.reserve(Nblock);

    for(decltype(Nblock) b = 0; b < Nblock; ++b)
        {
        auto& d = dvecs.at(b);
        auto& B = blocks[b];

        //print("d^2 = "); for(auto&& el : d) print(" ",sqr(el)); println();

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
    
    IQIndex L("L",move(Liq),uI.dir()), 
            R("R",move(Riq),vI.dir());

    IQIndexSet Uis(uI,dag(L)),
               Dis(L,R),
               Vis(vI,dag(R));

    //Print(Uis);
    //Print(Dis);
    //Print(Vis);

    IQTReal Ustore(Uis,QN()),
            Vstore(Vis,QN());

    IQTDiag Dstore(Dis,div(A));

    long n = 0;
    for(auto b: count(Nblock))
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
        //println("Doing Uref &= UU");
        //Print(Uref.range());
        //Print(UU.range());
        Uref &= UU;

        auto dind = stdx::make_array(n,n);
        auto pD = getBlock(Dstore,Dis,dind);
        assert(pD.data() != nullptr);
        auto Dref = makeVecRef(pD.data(),d.size());
        //println("Doing Dref &= d");
        //Print(Dref.range());
        //Print(d.range());
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
    D = IQTensor({L,R},move(Dstore),A.scale()*signfix);
    V = IQTensor(Vis,move(Vstore));

    //Originally eigs were found without including scale
    //so put the scale back in
    probs *= sqr(A.scale().real0());

    return Spectrum(move(probs),Args("Truncerr",truncerr));

    } // svdRank2 IQTensor


Spectrum
diag_hermitian(ITensor rho, 
               ITensor& U, 
               ITensor& D,
               const Args& args)
    {
    auto cutoff = args.getReal("Cutoff",MIN_CUT);
    auto maxm = args.getInt("Maxm",MAX_M);
    auto minm = args.getInt("Minm",1);
    auto do_truncate = args.getBool("Truncate",false);
    auto doRelCutoff = args.getBool("DoRelCutoff",false);
    auto absoluteCutoff = args.getBool("AbsoluteCutoff",false);
    auto cplx = isComplex(rho);
    auto showeigs = args.getBool("ShowEigs",false);

#ifdef DEBUG
    if(rho.r() != 2)
        {
        Print(rho.r());
        Print(rho);
        Error("Rank greater than 2 in diag_hermitian");
        }
#endif

    Index active;
    for(const Index& I : rho.inds())
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

    //Depending on the sign of the scale, calling .toMatrix11NoScale 
    //yields a matrix proportional to either rho or -rho.
    //If rho (scale().sign() > 0) then want to temporarily reverse 
    //the sign of the matrix when calling the diagonalization routine
    //to ensure eigenvalues are ordered from largest to smallest.
    if(rho.scale().sign() < 0) rho.scaleTo(rho.scale()*(-1));

    //Do the diagonalization
    Vector DD;
    Matrix UU,iUU;
    if(!cplx)
        {
        auto R = toMatrixRefc(rho,active,prime(active));
        diagSymmetric(R,UU,DD);
        }
    else
        {
        Error("Complex diag_hermitian not yet implemented");
        //Matrix Mr,Mi;
        //ITensor rrho = realPart(rho),
        //        irho = imagPart(rho);
        //rrho.scaleTo(rho.scale());
        //irho.scaleTo(rho.scale());
        //rrho.toMatrix11NoScale(prime(active),active,Mr);
        //irho.toMatrix11NoScale(prime(active),active,Mi);
        //if(flipSign)
        //    {
        //    Mr *= -1.0; 
        //    Mi *= -1.0; 
        //    }
        //HermitianEigenvalues(Mr,Mi,DD,UU,iUU); 
        //if(flipSign) DD *= -1.0;
        }


    //Include rho's scale to get the actual eigenvalues kept
    //as long as the leading eigenvalue is within a few orders
    //of magnitude of 1.0. Otherwise just report the scaled eigs.
    //Real orderMag = log(std::fabs(DD(1))) + rho.scale().logNum();
    //if(std::fabs(orderMag) < 5 && rho.scale().isFiniteReal())
    //    {
    //    DD *= rho.scale().real();
    //    }

    if(showeigs)
        {
        println("Before truncating, m = ",DD.size());
        println("DD = ",DD);
        printfln("maxm=%d,minm=%d,cutoff=%.2E",maxm,minm,cutoff);
        }

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
        if(showeigs)
            {
            printfln("Truncated to m=%d, trunc. err. = %.2E",m,truncerr);
            }
        }
    Spectrum spec;
    spec.truncerr(truncerr);

#ifdef DEBUG
    if(m > maxm)
        {
        printfln("m > maxm; m = %d, maxm = %d",m,maxm);
        Error("m > maxm");
        }
    if(m > 50000)
        {
        printfln("WARNING: very large m = %d in ITensor diag_hermitian");
        }
#endif

    if(args.getBool("ShowEigs",false))
        {
        printfln("\nminm = %d, maxm = %d, cutoff = %.3E",minm,maxm,cutoff);
        printfln("Kept %d states in diag_denmat",m);
        printfln("Truncation error = %.3E",truncerr);
        auto stop = DD.size();
        print("Eigs: ");
        for(decltype(stop) j = 1; j <= stop; ++j)
            {
            printf(DD(j) > 1E-3 ? ("%.3f") : ("%.3E"),DD(j));
            print((j != stop) ? ", " : "\n");
            }
        println();
        }

    Index newmid(active.rawname(),m,active.type());

    if(cplx)
        {
        //ITensor iU(active,newmid,iUU.Columns(1,m));
        //U = U + iU*Complex_i;
        }
    else
        {
        U = ITensor({active,newmid},ITReal{move(UU.storage())}); 
        D = ITensor({prime(newmid),newmid},ITDiag<Real>{DD.begin(),DD.end()},rho.scale());
        }

    if(rho.scale().isFiniteReal())
        DD *= rho.scale().real();
    else
        println("Scale not a finite Real, omitting from returned spectrum.");

    spec.eigsKept(move(DD));

    return spec;
    }

Spectrum
diag_hermitian(IQTensor rho, IQTensor& U, IQTensor& D,
               const Args& args)
    {
    /*
    const Real cutoff = args.getReal("Cutoff",MIN_CUT);
    const int maxm = args.getInt("Maxm",MAX_M);
    const int minm = args.getInt("Minm",1);
    const bool do_truncate = args.getBool("Truncate",false);
    const bool doRelCutoff = args.getBool("DoRelCutoff",false);
    const bool absoluteCutoff = args.getBool("AbsoluteCutoff",false);
    const bool cplx = rho.isComplex();

    if(rho.r() != 2)
        {
        Print(rho.indices());
        Error("Density matrix doesn't have rank 2");
        }

#ifdef DEBUG
    const QN Zero;
    if(div(rho) != Zero)
        { 
        Print(rho); 
        Error("Non-zero divergence of rho, QNs not conserved by Hamiltonian?");
        }
#endif

    vector<Matrix> mmatrix(rho.blocks().size()),
                   imatrix;
    vector<Vector> mvector(rho.blocks().size());
    vector<Real> alleig;
    alleig.reserve(rho.indices().front().m());

    if(cplx)
        imatrix.resize(rho.blocks().size());

    if(rho.indices().front().m() == 0)
        throw ResultIsZero("rho.index(1).m()");
    if(rho.empty())
        throw ResultIsZero("rho.empty()");

    LogNum refNorm(1);
    if(doRelCutoff)
        {
        //DO_IF_DEBUG(cout << "Doing relative cutoff\n";)
        Real maxLogNum = -200;
        rho.scaleOutNorm();
        for(const ITensor& t : rho.blocks())
            {
            maxLogNum = std::max(maxLogNum,t.scale().logNum());
            }
        refNorm = LogNum(maxLogNum,1);
        }
    rho.scaleTo(refNorm);


    //1. Diagonalize each ITensor within rho.
    //   Store results in mmatrix and mvector.
    int itenind = 0;
    for(const ITensor& t : rho.blocks())
        {
        Index a;
        for(const Index& I : t.indices())
            {
            if(I.primeLevel() == 0)
                {
                a = I;
                break;
                }
            }

        Matrix &UU = mmatrix.at(itenind);
        Vector &d =  mvector.at(itenind);

        //Depending on the sign of the scale, calling .toMatrix11NoScale 
        //yields a matrix proportional to either t or -t.
        //If t (scale().sign() > 0) then want to temporarily reverse 
        //the sign of the matrix when calling the diagonalization routine
        //to ensure eigenvalues are ordered from largest to smallest.
        bool flipSign = t.scale().sign() > 0;

        //Diag ITensors within rho
        const int n = a.m();
        if(!cplx)
            {
            Matrix M;
            t.toMatrix11NoScale(a,prime(a),M);
            if(flipSign) M *= -1;
            EigenValues(M,d,UU);
            if(flipSign) d *= -1;
            }
        else
            {
            ITensor ret = realPart(t),
                    imt = imagPart(t);
            ret.scaleTo(refNorm);
            imt.scaleTo(refNorm);
            Matrix Mr,Mi;
            Matrix &iUU = imatrix.at(itenind);
            ret.toMatrix11NoScale(prime(a),a,Mr);
            imt.toMatrix11NoScale(prime(a),a,Mi);
            if(flipSign)
                {
                Mr *= -1;
                Mi *= -1;
                }
            HermitianEigenvalues(Mr,Mi,d,UU,iUU);
            if(flipSign) d *= -1;
            }

        for(int j = 1; j <= n; ++j) 
            alleig.push_back(d(j));

        ++itenind;

#ifdef STRONG_DEBUG
	Real maxM = 1.0;
        for(int r = 1; r <= n; ++r)
	    for(int c = r+1; c <= n; ++c)
		maxM = std::max(maxM,std::fabs(M(r,c)));
	Real maxcheck = 1e-13 * maxM;
        for(int r = 1; r <= n; ++r)
	    for(int c = r+1; c <= n; ++c)
            {
            if(std::fabs(M(r,c)-M(c,r)) > maxcheck)
                {
                Print(M);
                Error("M not symmetric in diag_denmat");
                }
            }

        Matrix Id(UU.Nrows(),UU.Nrows()); Id = 1;
        Matrix Diff = Id-(UU.t()*UU);
        if(Norm(Diff.TreatAsVector()) > 1E-12)
            {
            printfln("\ndiff=%.2E",Norm(Diff.TreatAsVector()));
            Print(UU.t()*UU);
            Error("UU not unitary in diag_denmat");
            }
        
#endif //STRONG_DEBUG
        }

    //2. Truncate eigenvalues

    //Determine number of states to keep m
    Real svdtruncerr = 0;
    Real docut = -1;
    int m = (int)alleig.size();

    sort(alleig.begin(),alleig.end());

    if(do_truncate)
        {
        //Sort all eigenvalues from smallest to largest
        //irrespective of quantum numbers

        svdtruncerr = truncate(alleig,m,docut,maxm,minm,cutoff,absoluteCutoff,doRelCutoff);
        }
    Spectrum spec;
    spec.truncerr(svdtruncerr);

    if(args.getBool("ShowEigs",false))
        {
        cout << endl;
        printfln("Kept %d states in diag_denmat line 721", m);
        printfln("svdtruncerr = %.2E",svdtruncerr);
        printfln("docut = %.2E",docut);
        printfln("cutoff=%.2E, minm=%d, maxm=%d",cutoff,minm,maxm);
        cout << "doRelCutoff is " << (doRelCutoff ? "true" : "false") << endl;
        cout << "absoluteCutoff is " << (absoluteCutoff ? "true" : "false") << endl;
        cout << "refNorm is " << refNorm << endl;
        int s = alleig.size();
        const int max_show = 20;
        int stop = s-std::min(s,max_show);
        cout << "Eigs: ";
        for(int j = s-1; j >= stop; --j)
            {
            printf(alleig.at(j) > 1E-3 ? ("%.3f") : ("%.3E"), alleig.at(j));
            cout << ((j != stop) ? ", " : "\n");
            }
        cout << endl;
        }

#ifdef DEBUG
    if(m > maxm)
        {
        printfln("m > maxm; m = %d, maxm = %d",m,maxm);
        Error("m > maxm");
        }
    if(m > 20000)
        {
        cout << "WARNING: very large m = " << m << " in diag_hermitian" << endl;
        }
#endif

    IQIndex active;
    for(const IQIndex& I : rho.indices())
        {
        if(I.primeLevel() == 0)
            {
            active = I;
            break;
            }
        }

    //
    //Truncate eigenvalues and eigenvectors of rho
    //

    //Build blocks for unitary diagonalizing rho
    vector<ITensor> blocks,
                    iblocks;
    vector<ITensor> Dblocks;
    blocks.reserve(rho.blocks().size());
    Dblocks.reserve(rho.blocks().size());
    if(cplx) iblocks.reserve(rho.blocks().size());

    //Also form new Link IQIndex with appropriate m's for each block
    IQIndex::Storage iq;
    iq.reserve(rho.blocks().size());

    itenind = 0;
    for(const ITensor& t : rho.blocks())
        {
        Vec& thisD = mvector.at(itenind);
        Matrix& thisU = mmatrix.at(itenind);

        int this_m = 1;
        if(do_truncate)
            {
            while(this_m <= thisD.Length() && thisD(this_m) > docut) 
                {
                if(thisD(this_m) < 0) thisD(this_m) = 0;
                ++this_m;
                }
            --this_m; //since the loop overshoots by 1

            if(m == 0 && thisD.Length() >= 1) // zero mps, just keep one arb state
                { this_m = 1; m = 1; docut = 1; }

            thisD.ReduceDimension(this_m);
            }
        else
            {
            this_m = thisD.Length();
            }

        if(this_m == 0) { ++itenind; continue; }

        Index nm("qlink",this_m);

        Index act;
        for(const Index& I : t.indices())
            {
            if(I.primeLevel() == 0)
                {
                act = I;
                break;
                }
            }

        iq.push_back(IndexQN(nm,qn(active,act)));

        MatrixRef Utrunc = thisU.Columns(1,this_m);

        ITensor block(act,nm);
        block.fromMatrix11(act,nm,Utrunc);
        blocks.push_back(block);

        if(cplx)
            {
            iblocks.push_back(ITensor(act,nm,imatrix.at(itenind).Columns(1,this_m)));
            }

        Dblocks.push_back(ITensor(prime(nm),nm,thisD.SubVector(1,this_m)));

        ++itenind;
        }

    if(iq.size() == 0)
        {
        Print(m);
        Print(docut);
        throw ResultIsZero("iq.size() == 0");
        }

    IQIndex newmid("qlink",iq, -active.dir());

    U = IQTensor(dag(active),dag(newmid));
    D = IQTensor(prime(newmid),dag(newmid));
    for(long j = 0; j < blocks.size(); ++j)
        {
        D += Dblocks.at(j);
        U += blocks.at(j);
        }

    if(cplx)
        {
        IQTensor iU(dag(active),dag(newmid));
        for(size_t j = 0; j < iblocks.size(); ++j)
            {
            if(iblocks.at(j).norm() > 1E-14)
                iU += iblocks.at(j);
            }
        if(!iU.blocks().empty())
            U = U + iU*Complex_i;
        }

    D *= refNorm;

    Vec DD(newmid.m());
    const size_t aesize = alleig.size();
    for(int i = 1; i <= newmid.m(); ++i) 
        DD(i) = alleig.at(aesize-i);

    spec.eigsKept(DD);

    return spec;
    */

    //TODO: remove this, just here to make it compile
    return Spectrum();

    } //void diag_hermitian

void 
eig_decomp(ITensor T, 
           const Index& L, const Index& R,
           ITensor& V, ITensor& D,
           const Args& args)
    {
    /*
    //const bool doRelCutoff = args.getBool("DoRelCutoff",false);
    bool cplx = T.isComplex();

#ifdef DEBUG
    if(T.r() != 2)
        {
        Print(T.r());
        Print(T);
        Error("eig_decomp requires rank 2 tensor as input");
        }
#endif

    //Do the diagonalization
    Vec Dr,Di;
    Matrix Ur,Ui;
    if(!cplx)
        {
        Matrix M;
        T.toMatrix11NoScale(L,R,M);
        GenEigenValues(M,Dr,Di,Ur,Ui); 
        }
    else
        {
        Matrix Mr,Mi;
        ITensor rT = realPart(T),
                iT = imagPart(T);
        rT.scaleTo(T.scale());
        iT.scaleTo(T.scale());
        rT.toMatrix11NoScale(L,R,Mr);
        iT.toMatrix11NoScale(L,R,Mi);
        ComplexEigenvalues(Mr,Mi,Dr,Di,Ur,Ui); 
        }


    Index newmid("d",R.m(),R.type());
    V = ITensor(R,newmid,Ur);
    D = ITensor(prime(newmid),newmid,Dr);

    if(Norm(Ui.TreatAsVector()) > 1E-12)
        {
        V += ITensor(R,newmid,Ui)*Complex_i;
        }

    if(Norm(Di) > 1E-12)
        {
        D += ITensor(prime(newmid),newmid,Di)*Complex_i;
        }

    D *= T.scale();
    */

    }

void 
eig_decomp(IQTensor T, 
           const IQIndex& L, const IQIndex& R,
           IQTensor& V, IQTensor& D,
           const Args& args)
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

} //namespace itensor
