//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#include "svdalgs.h"
#include <algorithm>

namespace itensor {

using std::swap;
using std::istream;
using std::ostream;
using std::cout;
using std::endl;
using std::vector;
using std::find;
using std::pair;
using std::make_pair;
using std::string;
using std::sqrt;

//Vector static
//sqrt(Vector V)
//    {
//    for(int j = 1; j <= V.Length(); ++j)
//        V(j) = sqrt(fabs(V(j)));
//    return V;
//    }


Real static
truncate(Vector& D,
         int maxm,
         int minm,
         Real cutoff,
         bool absoluteCutoff,
         bool doRelCutoff)
    {
    int m = D.Length();
    if(m == 1) return 0;

    Real truncerr = 0;

    //Zero out any negative weight
    for(int zerom = m; zerom > 0; --zerom)
        {
        if(D(zerom) >= 0) break;
        D(zerom) = 0;
        }

    if(absoluteCutoff)
        {
        for(;m > maxm || (D(m) < cutoff && m > minm); --m)
            {
            truncerr += D(m);
            }
        }
    else
        {
        const Real scale = doRelCutoff ? D(1) : 1.0;
        for(;m > maxm || (truncerr+D(m) < cutoff*scale && m > minm); --m)
            {
            truncerr += D(m);
            }
        truncerr = (D(1) == 0 ? 0 : truncerr/scale);
        }

    D.ReduceDimension(m); 

    return truncerr;
    }

Real static
truncate(vector<Real>& alleig, 
         int& m, 
         Real& docut, 
         int maxm,
         int minm,
         Real cutoff,
         bool absoluteCutoff,
         bool doRelCutoff)
    {
    m = (int)alleig.size();
    if(m == 1)
        {
        docut = alleig.front()/2.;
        return 0;
        }
    int mdisc = 0;

    Real truncerr = 0;

    if(absoluteCutoff)
        {
        while(m > maxm || ( (alleig.at(mdisc) < cutoff && m > minm)
            && mdisc < (int)alleig.size() ) )
            {
            if(alleig.at(mdisc) > 0)
                truncerr += alleig.at(mdisc);
            else
                alleig.at(mdisc) = 0;

            ++mdisc;
            --m;
            }
        docut = (mdisc > 0 
                ? (alleig.at(mdisc-1) + alleig.at(mdisc))*0.5 - 1E-5*alleig.at(mdisc-1)
                : -1);
        }
    else
        {
        Real scale = doRelCutoff ? alleig.back() : 1.0;
        while(   m > maxm 
             || ( (mdisc < (int)alleig.size()) && (truncerr+alleig.at(mdisc) < cutoff*scale && m > minm))
             )
            {
            if(alleig.at(mdisc) > 0)
                truncerr += alleig.at(mdisc);
            else
                alleig.at(mdisc) = 0;

            ++mdisc;
            --m;
            }
        if(mdisc >= int(alleig.size())) mdisc = alleig.size() - 1;
        docut = (mdisc > 0 
                ? (alleig.at(mdisc-1) + alleig.at(mdisc))*0.5 - 1E-5*alleig.at(mdisc-1)
                : -1);
        truncerr = (alleig.back() == 0 ? 0 : truncerr/scale);
        }


    return truncerr;
    }



Spectrum 
svdRank2(ITensor A, const Index& ui, const Index& vi,
         ITensor& U, ITensor& D, ITensor& V,
         const OptSet& opts)
    {
    const Real thresh = opts.getReal("SVDThreshold",1E-4);
    const Real cutoff = opts.getReal("Cutoff",MIN_CUT);
    const int maxm = opts.getInt("Maxm",MAX_M);
    const int minm = opts.getInt("Minm",1);
    const bool do_truncate = opts.getBool("Truncate",true);
    const bool doRelCutoff = opts.getBool("DoRelCutoff",false);
    const bool absoluteCutoff = opts.getBool("AbsoluteCutoff",false);
    const bool cplx = A.isComplex();

    if(A.r() != 2)
        {
        Error("A must be matrix-like");
        }

    Matrix UU,VV,
           iUU,iVV;
    Vector DD;

    if(!cplx)
        {
        Matrix M;
        A.toMatrix11NoScale(ui,vi,M);

        SVD(M,UU,DD,VV,thresh);
        }
    else
        {
        ITensor Are = realPart(A),
                Aim = imagPart(A);
        Are.scaleTo(A.scale());
        Aim.scaleTo(A.scale());
        Matrix Mre,Mim;
        Are.toMatrix11NoScale(ui,vi,Mre);
        Aim.toMatrix11NoScale(ui,vi,Mim);

        //SVDComplex(Mre,Mim,UU,iUU,DD,VV,iVV);
        SVD(Mre,Mim,UU,iUU,DD,VV,iVV,thresh);
        }

    //Truncate

    Spectrum spec;

    int m = DD.Length();
    Real terr = 0;
    if(do_truncate)
        {
        Vector sqrD(DD);
        for(int j = 1; j <= sqrD.Length(); ++j)
            sqrD(j) = sqr(DD(j));
        terr = truncate(sqrD,maxm,minm,cutoff,absoluteCutoff,doRelCutoff);
        m = sqrD.Length();
        DD.ReduceDimension(m);
        }


    if(opts.getBool("ShowEigs",false))
        {
        cout << endl;
        printfln("minm = %d, maxm = %d, cutoff = %.3E",minm,maxm,cutoff);
        printfln("truncate = %s",(do_truncate?"true":"false"));
        printfln("doRelCutoff = %s",(doRelCutoff?"true":"false"));
        printfln("absoluteCutoff = %s",(absoluteCutoff?"true":"false"));
        printfln("Kept m=%d states in svdRank2 line 169", m);
        printfln("svdtruncerr = %.3E",spec.truncerr());

        int stop = min(10,DD.Length());
        Vector Ds = DD.SubVector(1,stop);

        Real orderMag = log(fabs(DD(1))) + A.scale().logNum();
        if(fabs(orderMag) < 5 && A.scale().isFiniteReal())
            {
            Ds *= fabs(A.scale().real0());
            cout << "Singular values: ";
            }
        else
            {
            cout << "Singular values (not including scale = " << A.scale() << "):";
            }

        for(int j = 1; j <= stop; ++j)
            {
            const Real sval = Ds(j);
            printf(( sval > 1E-3 && sval < 1000) ? ("%.3f") : ("%.3E") , sval); 
            print((j != stop) ? ", " : "\n");
            }
        cout << endl;
        }
    
    Index uL("ul",m,Link),vL("vl",m,Link);

    D = ITensor(uL,vL,DD);
    D *= A.scale();
    U = ITensor(ui,uL,UU.Columns(1,m));
    V = ITensor(vL,vi,VV.Rows(1,m));

    if(cplx)
        {
        ITensor iU(ui,uL,iUU.Columns(1,m)),
                iV(vL,vi,iVV.Rows(1,m));
        if(iU.norm() > 1E-14)
            U = U + iU*Complex_i;
        if(iV.norm() > 1E-14)
            V = V + iV*Complex_i;
        }

    //Square all singular values
    //since convention is to report
    //density matrix eigs
    for(int j = 1; j <= m; ++j)
        {
        DD(j) *= DD(j);
        }

    if(A.scale().isFiniteReal())
        {
        DD *= sqr(A.scale().real0());
        }
    else
        {
        cout << "Warning: scale not finite real" << endl;
        }

    //Global::lastd() = DD;

    //Include A's scale to get the actual eigenvalues kept
    //as long as the leading eigenvalue is within a few orders
    //of magnitude of 1.0. Otherwise just report the scaled eigs.
    //Real orderMag = log(fabs(DD(1))) + A.scale().logNum();
    //if(fabs(orderMag) < 5 && A.scale().isFiniteReal())
    //    {
    //    Global::lastd() *= A.scale().real();
    //    }

    return Spectrum(DD,Opt("Truncerr",terr));

    } // void svdRank2

Spectrum
svdRank2(IQTensor A, const IQIndex& uI, const IQIndex& vI,
         IQTensor& U, IQTensor& D, IQTensor& V,
         const OptSet& opts)
    {
    const bool cplx = A.isComplex();
    const Real thresh = opts.getReal("SVDThreshold",1E-4);
    const Real cutoff = opts.getReal("Cutoff",MIN_CUT);
    const int maxm = opts.getInt("Maxm",MAX_M);
    const int minm = opts.getInt("Minm",1);
    const bool do_truncate = opts.getBool("Truncate",true);
    const bool doRelCutoff = opts.getBool("DoRelCutoff",false);
    const bool absoluteCutoff = opts.getBool("AbsoluteCutoff",false);
    const Real logrefNorm = opts.getReal("LogRefNorm",0.);

    if(A.r() != 2)
        {
        Error("A must be matrix-like");
        }

    const int Nblock = A.blocks().size();
    if(Nblock == 0)
        throw ResultIsZero("A has no blocks");

    vector<Matrix> Umatrix(Nblock),
                   Vmatrix(Nblock),
                   iUmatrix,
                   iVmatrix;
    if(cplx)
        {
        iUmatrix.resize(Nblock);
        iVmatrix.resize(Nblock);
        }

    vector<Vector> dvector(Nblock);

    vector<Real> alleig;
    alleig.reserve(min(uI.m(),vI.m()));

    if(uI.m() == 0)
        throw ResultIsZero("uI.m() == 0");
    if(vI.m() == 0)
        throw ResultIsZero("vI.m() == 0");

    LogNumber refNorm(logrefNorm,1);
    if(doRelCutoff)
        {
        Real maxLogNum = -200;
        A.scaleOutNorm();
        Foreach(const ITensor& t, A.blocks())
            {
            maxLogNum = max(maxLogNum,t.scale().logNum());
            }
        refNorm = LogNumber(maxLogNum,1);
        }
    A.scaleTo(refNorm);

    //1. SVD each ITensor within A.
    //   Store results in mmatrix and mvector.
    int itenind = 0;
    Foreach(const ITensor& t, A.blocks())
        {
        Matrix &UU = Umatrix.at(itenind);
        Matrix &VV = Vmatrix.at(itenind);
        Vector &d =  dvector.at(itenind);

        const Index *ui=0,*vi=0;
        bool gotui = false;
        Foreach(const Index& I, t.indices())
            {
            if(!gotui) 
                {
                ui = &I;
                gotui = true;
                }
            else       
                {
                vi = &I;
                break;
                }
            }

        if(!hasindex(uI,*ui))
            swap(ui,vi);

        if(!cplx)
            {
            Matrix M(ui->m(),vi->m());
            t.toMatrix11NoScale(*ui,*vi,M);

            SVD(M,UU,d,VV,thresh);
            }
        else
            {
            ITensor ret = realPart(t),
                    imt = imagPart(t);
            ret.scaleTo(refNorm);
            imt.scaleTo(refNorm);
            Matrix Mre(ui->m(),vi->m()),
                   Mim(ui->m(),vi->m());
            ret.toMatrix11NoScale(*ui,*vi,Mre);
            imt.toMatrix11NoScale(*ui,*vi,Mim);

            //SVDComplex(Mre,Mim,
            //           UU,iUmatrix.at(itenind),
            //           d,
            //           VV,iVmatrix.at(itenind));
            SVD(Mre,Mim,
                UU,iUmatrix.at(itenind),
                d,
                VV,iVmatrix.at(itenind),
                thresh);
            }

        //Store the squared singular values
        //(denmat eigenvalues) in alleig
        for(int j = 1; j <= d.Length(); ++j) 
            alleig.push_back(sqr(d(j)));

        ++itenind;
        }

    //2. Truncate eigenvalues

    //Determine number of states to keep m
    int m = (int)alleig.size();
    Real svdtruncerr = 0;
    Real docut = -1;

    if(do_truncate)
        {
        //Sort all eigenvalues from smallest to largest
        //irrespective of quantum numbers
        sort(alleig.begin(),alleig.end());

        svdtruncerr = truncate(alleig,m,docut,maxm,minm,cutoff,
                               absoluteCutoff,doRelCutoff);
        }

    if(opts.getBool("ShowEigs",false))
        {
        cout << endl;
        println("svdRank2 (IQTensor):");
        printfln("    minm = %d, maxm = %d, cutoff = %.3E",minm,maxm,cutoff);
        printfln("    Kept m = %d states in svdRank2",m);
        printfln("    svdtruncerr = %.2E",svdtruncerr);
        printfln("    docut = %.2E",docut);
        println("    doRelCutoff is ",(doRelCutoff ? "true" : "false"));
        println("    absoluteCutoff is ",(absoluteCutoff ? "true" : "false"));
        println("    refNorm is ",refNorm);

        const int s = alleig.size();
        const int max_show = 20;
        int stop = s-min(s,max_show);

        //Include refNorm in printed eigs as long as
        //the leading eig is within a few orders of magnitude
        //of 1.0. Otherwise just print the scaled eigs.
        Real orderMag = log(fabs(alleig.at(s-1))) + refNorm.logNum();
        Real real_fac = 1;
        if(fabs(orderMag) < 5 && refNorm.isFiniteReal())
            {
            real_fac = refNorm.real();
            cout << "    Singular values: ";
            }
        else
            {
            cout << "    Singular values [omitting scale factor " << refNorm << "]: \n";
            if(alleig.at(s-1) > 1.e10)
                {
                Error("bad alleig");
                }
            cout << "    ";
            }

        for(int j = s-1; j >= stop; --j)
            {
            const Real sval = sqrt(alleig.at(j))*real_fac;
            printf( (sval >= 1E-3 && sval < 1E3) ? ("%.3f") : ("%.3E"), sval);
            print((j != stop) ? ", " : "\n");
            }
        cout << endl;
        } //end if(showeigs_)

    //Truncate denmat eigenvalue vectors
    //Also form new Link index with appropriate m's for each block
    IQIndex::Storage Liq, Riq;
    Liq.reserve(Nblock);
    Riq.reserve(Nblock);

    vector<ITensor> Ublock,
                    Vblock,
                    iUblock,
                    iVblock;
    Ublock.reserve(Nblock);
    Vblock.reserve(Nblock);
    if(cplx)
        {
        iUblock.reserve(Nblock);
        iVblock.reserve(Nblock);
        }

    vector<ITensor> Dblock;
    Dblock.reserve(Nblock);

    itenind = 0;
    int total_m = 0;
    Foreach(const ITensor& t, A.blocks())
        {
        const Matrix& UU = Umatrix.at(itenind);
        const Matrix& VV = Vmatrix.at(itenind);
        Vector& thisD = dvector.at(itenind);

        int this_m = 1;
        while(this_m <= thisD.Length() && sqr(thisD(this_m)) > docut) 
            {
            ++total_m;
            if(thisD(this_m) < 0) thisD(this_m) = 0;
            ++this_m;
            }
        --this_m; //since the loop overshoots by 1

        if(m == 0 && thisD.Length() >= 1) // zero mps, just keep one arb state
            { this_m = 1; m = 1; docut = 1; }

        if(this_m == 0) { ++itenind; continue; }

        const Index *ui=0,*vi=0;
        bool gotui = false;
        Foreach(const Index& I, t.indices())
            {
            if(!gotui) 
                {
                ui = &I;
                gotui = true;
                }
            else       
                {
                vi = &I;
                break;
                }
            }

        if(!hasindex(uI,*ui))
            swap(ui,vi);

        Index l("l",this_m);
        Liq.push_back(IndexQN(l,qn(uI,*ui)));

        Index r("r",this_m);
        Riq.push_back(IndexQN(r,qn(vI,*vi)));

        Dblock.push_back(ITensor(l,r,thisD.SubVector(1,this_m)));

        Ublock.push_back(ITensor(*ui,l,UU.Columns(1,this_m)));
        Vblock.push_back(ITensor(r,*vi,VV.Rows(1,this_m)));

        if(cplx)
            {
            iUblock.push_back(ITensor(*ui,l,iUmatrix.at(itenind).Columns(1,this_m)));
            iVblock.push_back(ITensor(r,*vi,iVmatrix.at(itenind).Rows(1,this_m)));
            }

        ++itenind;
        }

    if(Liq.size() == 0)
        throw ResultIsZero("Liq.size() == 0");

    IQIndex L("L",Liq,uI.dir()), R("R",Riq,vI.dir());

    D = IQTensor(L,R);
    U = IQTensor(uI,dag(L));
    V = IQTensor(dag(R),vI);

    //Load blocks into D,U, and V
    for(size_t j = 0; j < Dblock.size(); ++j)
        {
        D += Dblock.at(j);
        U += Ublock.at(j);
        V += Vblock.at(j);
        }

    if(cplx)
        {
        IQTensor iU(uI,dag(L));
        IQTensor iV(dag(R),vI);
        for(size_t j = 0; j < Dblock.size(); ++j)
            {
            if(iUblock.at(j).norm() > 1E-14)
                {
                iU += iUblock.at(j);
                }
            if(iVblock.at(j).norm() > 1E-14)
                {
                iV += iVblock.at(j);
                }
            }
        if(!iU.blocks().empty())
            {
            U = U + iU*Complex_i;
            }
        if(!iV.blocks().empty())
            {
            V = V + iV*Complex_i;
            }
        }

    //Originally eigs were found by calling
    //toMatrix11NoScale, so put the scale back in
    D *= refNorm;

//    Vector DD(L.m());
//    for(int i = 1; i <= L.m(); ++i)
//        {
//        DD(i) = alleig.at(alleig.size()-i);
//        }

    return Spectrum(D,Opt("Truncerr",svdtruncerr));

    } //void svdRank2


Spectrum
diag_hermitian(ITensor rho, ITensor& U, ITensor& D,
               const OptSet& opts)
    {
    const Real cutoff = opts.getReal("Cutoff",MIN_CUT);
    const int maxm = opts.getInt("Maxm",MAX_M);
    const int minm = opts.getInt("Minm",1);
    const bool do_truncate = opts.getBool("Truncate",true);
    const bool doRelCutoff = opts.getBool("DoRelCutoff",false);
    const bool absoluteCutoff = opts.getBool("AbsoluteCutoff",false);
    const bool cplx = rho.isComplex();

#ifdef DEBUG
    if(rho.r() != 2)
        {
        Print(rho.r());
        Print(rho);
        Error("Too many indices for density matrix");
        }
#endif

    Index active;
    Foreach(const Index& I, rho.indices())
        {
        if(I.primeLevel() == 0)
            {
            active = I;
            break;
            }
        }

    if(!active)
        {
        Print(rho.indices());
        Error("Tensor must have one unprimed index");
        }

    //Do the diagonalization
    Vector DD;
    Matrix UU,iUU;
    if(!cplx)
        {
        Matrix R;
        rho.toMatrix11NoScale(active,prime(active),R);
        R *= -1.0; 
        EigenValues(R,DD,UU); 
        DD *= -1.0;
        }
    else
        {
        Matrix Mr,Mi;
        ITensor rrho = realPart(rho),
                irho = imagPart(rho);
        rrho.scaleTo(rho.scale());
        irho.scaleTo(rho.scale());
        rrho.toMatrix11NoScale(prime(active),active,Mr);
        irho.toMatrix11NoScale(prime(active),active,Mi);
        Mr *= -1.0; 
        Mi *= -1.0; 
        HermitianEigenvalues(Mr,Mi,DD,UU,iUU); 
        DD *= -1.0;
        }


    //Include rho's scale to get the actual eigenvalues kept
    //as long as the leading eigenvalue is within a few orders
    //of magnitude of 1.0. Otherwise just report the scaled eigs.
    //Real orderMag = log(fabs(DD(1))) + rho.scale().logNum();
    //if(fabs(orderMag) < 5 && rho.scale().isFiniteReal())
    //    {
    //    DD *= rho.scale().real();
    //    }

    //Truncate
    Real svdtruncerr = 0.0;
    if(opts.getBool("ShowEigs",false))
        cout << "Before truncating, m = " << DD.Length() << endl;
    if(do_truncate)
        {
        svdtruncerr = truncate(DD,maxm,minm,cutoff,absoluteCutoff,doRelCutoff);
        }
    Spectrum spec;
    spec.truncerr(svdtruncerr);
    const int m = DD.Length();

#ifdef DEBUG
    if(m > maxm)
        {
        printfln("m > maxm; m = %d, maxm = %d",m,maxm);
        Error("m > maxm");
        }
    if(m > 20000)
        {
        cout << "WARNING: very large m = " << m << " in ITensor diag_hermitian" << endl;
        }
#endif

    if(opts.getBool("ShowEigs",false))
        {
        cout << endl;
        printfln("minm = %d, maxm = %d, cutoff = %.3E",minm,maxm,cutoff);
        printfln("Kept %d states in diag_denmat",m);
        printfln("svdtruncerr = %.3E",svdtruncerr);
        //cout << "doRelCutoff is " << doRelCutoff << endl;
        //int stop = min(D.Length(),10);
        int stop = DD.Length();
        cout << "Eigs: ";
        for(int j = 1; j <= stop; ++j)
            {
            printf(DD(j) > 1E-3 ? ("%.3f") : ("%.3E"),DD(j));
            print((j != stop) ? ", " : "\n");
            }
        cout << endl;
        }

    Index newmid(active.rawname(),m,active.type());
    U = ITensor(active,newmid,UU.Columns(1,m));
    D = ITensor(prime(newmid),newmid,DD);
    D *= rho.scale();

    if(cplx)
        {
        ITensor iU(active,newmid,iUU.Columns(1,m));
        U = U + iU*Complex_i;
        }

    spec.eigsKept(DD);

    return spec;
    }

Spectrum
diag_hermitian(IQTensor rho, IQTensor& U, IQTensor& D,
               const OptSet& opts)
    {
    const Real cutoff = opts.getReal("Cutoff",MIN_CUT);
    const int maxm = opts.getInt("Maxm",MAX_M);
    const int minm = opts.getInt("Minm",1);
    const bool do_truncate = opts.getBool("Truncate",true);
    const bool doRelCutoff = opts.getBool("DoRelCutoff",false);
    const bool absoluteCutoff = opts.getBool("AbsoluteCutoff",false);
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

    LogNumber refNorm(1);
    if(doRelCutoff)
        {
        //DO_IF_DEBUG(cout << "Doing relative cutoff\n";)
        Real maxLogNum = -200;
        rho.scaleOutNorm();
        Foreach(const ITensor& t, rho.blocks())
            {
            maxLogNum = max(maxLogNum,t.scale().logNum());
            }
        refNorm = LogNumber(maxLogNum,1);
        }
    rho.scaleTo(refNorm);

    //1. Diagonalize each ITensor within rho.
    //   Store results in mmatrix and mvector.
    int itenind = 0;
    Foreach(const ITensor& t, rho.blocks())
        {
        Index a;
        Foreach(const Index& I, t.indices())
            {
            if(I.primeLevel() == 0)
                {
                a = I;
                break;
                }
            }

        Matrix &UU = mmatrix.at(itenind);
        Vector &d =  mvector.at(itenind);

        //Diag ITensors within rho
        const int n = a.m();
        if(!cplx)
            {
            Matrix M;
            t.toMatrix11NoScale(a,prime(a),M);
            M *= -1;
            EigenValues(M,d,UU);
            d *= -1;
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
            Mr *= -1;
            Mi *= -1;
            HermitianEigenvalues(Mr,Mi,d,UU,iUU);
            d *= -1;
            }

        for(int j = 1; j <= n; ++j) 
            alleig.push_back(d(j));

        ++itenind;

#ifdef STRONG_DEBUG
	Real maxM = 1.0;
        for(int r = 1; r <= n; ++r)
	    for(int c = r+1; c <= n; ++c)
		maxM = max(maxM,fabs(M(r,c)));
	Real maxcheck = 1e-13 * maxM;
        for(int r = 1; r <= n; ++r)
	    for(int c = r+1; c <= n; ++c)
            {
            if(fabs(M(r,c)-M(c,r)) > maxcheck)
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

    if(do_truncate)
        {
        //Sort all eigenvalues from smallest to largest
        //irrespective of quantum numbers
        sort(alleig.begin(),alleig.end());

        svdtruncerr = truncate(alleig,m,docut,maxm,minm,cutoff,absoluteCutoff,doRelCutoff);
        }
    Spectrum spec;
    spec.truncerr(svdtruncerr);

    if(opts.getBool("ShowEigs",false))
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
        int stop = s-min(s,max_show);
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
    Foreach(const IQIndex& I, rho.indices())
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
    Foreach(const ITensor& t, rho.blocks())
        {
        Vector& thisD = mvector.at(itenind);
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
        Foreach(const Index& I, t.indices())
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
    for(size_t j = 0; j < blocks.size(); ++j)
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

    Vector DD(newmid.m());
    const size_t aesize = alleig.size();
    for(int i = 1; i <= newmid.m(); ++i) 
        DD(i) = alleig.at(aesize-i);

    spec.eigsKept(DD);

    //Include spec.refNorm() to get the actual eigenvalues kept
    //as long as the leading eigenvalue is within a few orders
    //of magnitude of 1.0. Otherwise just report the scaled eigs.
    /*
    Real orderMag = log(fabs(D(1))) + spec.refNorm().logNum();
    if(fabs(orderMag) < 5 && spec.refNorm().isFiniteReal())
        {
        Global::lastd() *= spec.refNorm().real();
        }
        */

    return spec;

    } //void diag_hermitian

void 
eig_decomp(ITensor T, 
           const Index& L, const Index& R,
           ITensor& V, ITensor& D,
           const OptSet& opts)
    {
    //const bool doRelCutoff = opts.getBool("DoRelCutoff",false);
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
    Vector Dr,Di;
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

    }

void 
eig_decomp(IQTensor T, 
           const IQIndex& L, const IQIndex& R,
           IQTensor& V, IQTensor& D,
           const OptSet& opts)
    {
    const bool doRelCutoff = opts.getBool("DoRelCutoff",false);
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
    vector<Vector> reigs(nblocks),
                   ieigs(nblocks);

    if(T.empty())
        throw ResultIsZero("T has no blocks");

    LogNumber refNorm(1);
    if(doRelCutoff)
        {
        Real maxLogNum = -200;
        T.scaleOutNorm();
        Foreach(const ITensor& t, T.blocks())
            {
            maxLogNum = max(maxLogNum,t.scale().logNum());
            }
        refNorm = LogNumber(maxLogNum,1);
        }
    T.scaleTo(refNorm);

    //1. Diagonalize each ITensor within rho.
    //   Store results in mmatrix and mvector.
    int itenind = 0;
    Foreach(const ITensor& t, T.blocks())
        {
        Index li = t.indices().front(),
              ri = t.indices().back();

        if(!hasindex(L,li))
            swap(li,ri);

        Matrix &Ur = rmatrix.at(itenind),
               &Ui = imatrix.at(itenind);
        Vector &dr = reigs.at(itenind),
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
    Foreach(const ITensor& t, T.blocks())
        {
        Vector &dr = reigs.at(itenind),
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
    Foreach(const ITensor& t, Vblocks)
        {
        V += t;
        }

    D = IQTensor(prime(newmid),dag(newmid));
    Foreach(const ITensor& t, Dblocks)
        {
        D += t;
        }

    D *= refNorm;

    }

}; //namespace itensor
