//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#include "svdworker.h"
#include "localop.h"

using namespace std;
using boost::format;

inline Vector
sqrt(Vector V)
    {
    for(int j = 1; j <= V.Length(); ++j)
        V(j) = sqrt(fabs(V(j)));
    return V;
    }

SVDWorker::
SVDWorker(const OptSet& opts) 
    : 
    refNorm_(1)
    { 
    initOpts(opts);

    N_ = opts.getInt("N",1);
    
    truncerr_ = vector<Real>(N_+1,NAN);
    eigsKept_ = vector<Vector>(N_+1);
    }

SVDWorker::
SVDWorker(int N, const OptSet& opts)
    : 
    N_(N), 
    truncerr_(N_+1,NAN), 
    refNorm_(1), 
    eigsKept_(N_+1)
    { 
    initOpts(opts);
    }

void SVDWorker::
initOpts(const OptSet& opts)
    {
    absoluteCutoff_ = opts.getBool("AbsoluteCutoff",false);
    cutoff_ = opts.getReal("Cutoff",MIN_CUT);
    doRelCutoff_ = opts.getBool("DoRelCutoff",false);
    maxm_ = opts.getInt("Maxm",MAX_M);
    minm_ = opts.getInt("Minm",1);
    noise_ = opts.getReal("Noise",0.);
    showeigs_ = opts.getBool("ShowEigs",false);
    use_orig_m_ = opts.getBool("UseOrigM",false);
    }

SVDWorker::
SVDWorker(int N, Real cutoff, int minm, int maxm, 
          bool doRelCutoff, const LogNumber& refNorm)
    : N_(N), 
      truncerr_(N_+1), 
      cutoff_(cutoff), 
      minm_(minm), 
      maxm_(maxm),
      use_orig_m_(false), 
      showeigs_(false), 
      doRelCutoff_(doRelCutoff),
      absoluteCutoff_(false), 
      refNorm_(refNorm), 
      eigsKept_(N_+1),
      noise_(0)
    { }


template <class Tensor, class SparseT> 
void SVDWorker::
csvd(int b, const Tensor& AA, Tensor& L, SparseT& V, Tensor& R)
    { 
    csvd<Tensor>(b,AA,L,V,R,LocalOp<Tensor>::Null()); 
    }
template
void SVDWorker::
csvd(int b, const ITensor& AA, ITensor& L, ITSparse& V, ITensor& R);
template
void SVDWorker::
csvd(int b, const IQTensor& AA, IQTensor& L, IQTSparse& V, IQTensor& R);


template <class Tensor> 
void SVDWorker::
denmatDecomp(int b, const Tensor& AA, Tensor& A, Tensor& B, Direction dir)
    { 
    denmatDecomp<Tensor>(b,AA,A,B,dir,LocalOp<Tensor>::Null()); 
    }
template
void SVDWorker::
denmatDecomp(int b, const ITensor& AA, ITensor& A, ITensor& B, Direction dir);
template
void SVDWorker::
denmatDecomp(int b, const IQTensor& AA, IQTensor& A, IQTensor& B, Direction dir);



template <class Tensor, class SparseT>
void SVDWorker::
svd(int b, const Tensor& AA, Tensor& U, SparseT& D, Tensor& V)
    { 
    svd<Tensor>(b,AA,U,D,V,LocalOp<Tensor>::Null()); 
    }
template
void SVDWorker::
svd(int b, const ITensor& AA, ITensor& U, ITSparse& D, ITensor& V);
template
void SVDWorker::
svd(int b, const IQTensor& AA, IQTensor& U, IQTSparse& D, IQTensor& V);



void SVDWorker::
svdRank2(ITensor A, const Index& ui, const Index& vi,
         ITensor& U, ITSparse& D, ITensor& V, int b)
    {
    if(A.r() != 2)
        {
        Error("A.r() must be 2");
        }

    //const Index &ui = (U.hasindex(A.index(1)) ? A.index(1) : A.index(2)),
    //            &vi = (V.hasindex(A.index(2)) ? A.index(2) : A.index(1));

    Matrix M;
    A.toMatrix11NoScale(ui,vi,M);

    Matrix UU,VV;
    Vector& DD = eigsKept_.at(b);
    SVD(M,UU,DD,VV);

    //Truncate
    int m = DD.Length();
    Real& svdtruncerr = truncerr_.at(b);
    svdtruncerr = 0;

    //Zero out any negative weight
    for(int zerom = m; zerom > 0; --zerom)
        {
        if(DD(zerom) >= 0) break;
        else              DD(zerom) = 0;
        }

    if(absoluteCutoff_)
        {
        while(m > maxm_ || (sqr(DD(m)) < cutoff_ && m > minm_))
            {
            svdtruncerr += sqr(DD(m--));
            }
        }
    else
        {
        Real scale = doRelCutoff_ ? sqr(DD(1)) : 1.0;
        while(m > maxm_ || (svdtruncerr+sqr(DD(m)) < cutoff_*scale && m > minm_))
            {
            svdtruncerr += sqr(DD(m--));
            }
        svdtruncerr = (DD(1) == 0 ? 0 : svdtruncerr/scale);
        }

    DD.ReduceDimension(m); 

    if(showeigs_)
        {
        cout << endl;
        cout << format("minm_ = %d, maxm_ = %d, cutoff_ = %.3E")
                       %minm_%maxm_%cutoff_ << endl;
        cout << format("use_orig_m_ = %s")%(use_orig_m_?"true":"false")<<endl;
        cout << format("Kept m=%d states in svdRank2 line 169") % m << endl;
        cout << format("svdtruncerr = %.3E")%svdtruncerr << endl;


        int stop = min(10,DD.Length());
        Vector Ds = DD.SubVector(1,stop);

        Real orderMag = log(fabs(DD(1))) + A.scale().logNum();
        if(fabs(orderMag) < 5 && A.scale().isFiniteReal())
            {
            Ds *= A.scale().real();
            cout << "Eigs: ";
            }
        else
            {
            cout << "Eigs (not including scale = " << A.scale() << "):";
            }

        for(int j = 1; j <= stop; ++j)
            {
            const Real eig = sqr(Ds(j));
            cout << format( (eig > 1E-3 && eig < 1000) ? ("%.3f") : ("%.3E")) 
                    % eig;
            cout << ((j != stop) ? ", " : "\n");
            }
        cout << endl;
        }
    
    Index uL("ul",m,Link),vL("vl",m,Link);

    D = ITSparse(uL,vL,DD);
    D *= A.scale();
    U = ITensor(ui,uL,UU.Columns(1,m));
    V = ITensor(vL,vi,VV.Rows(1,m));

    //Square all singular values
    //since convention is to report
    //density matrix eigs
    for(int j = 1; j <= m; ++j)
        {
        DD(j) *= DD(j);
        }

    Global::lastd() = DD;

    //Include A's scale to get the actual eigenvalues kept
    //as long as the leading eigenvalue is within a few orders
    //of magnitude of 1.0. Otherwise just report the scaled eigs.
    Real orderMag = log(fabs(DD(1))) + A.scale().logNum();
    if(fabs(orderMag) < 5 && A.scale().isFiniteReal())
        {
        Global::lastd() *= A.scale().real();
        }

    } // void SVDWorker::svdRank2

void SVDWorker::
svdRank2(IQTensor A, const IQIndex& uI, const IQIndex& vI,
         IQTensor& U, IQTSparse& D, IQTensor& V, int b)
    {
    if(A.r() != 2)
        {
        Error("A.r() must be 2");
        }

    const int Nblock = A.iten_size();
    if(Nblock == 0)
        throw ResultIsZero("A.iten_size() == 0");

    vector<Matrix> Umatrix(Nblock),
                   Vmatrix(Nblock);
    vector<Vector> dvector(Nblock);

    vector<Real> alleig;
    alleig.reserve(min(uI.m(),vI.m()));

    if(uI.m() == 0)
        throw ResultIsZero("uI.m() == 0");
    if(vI.m() == 0)
        throw ResultIsZero("vI.m() == 0");

    if(doRelCutoff_)
        {
        Real maxLogNum = -200;
        A.scaleOutNorm();
        Foreach(const ITensor& t, A.blocks())
            {
            maxLogNum = max(maxLogNum,t.scale().logNum());
            }
        refNorm_ = LogNumber(maxLogNum,1);
        }

    A.scaleTo(refNorm_);

    //1. SVD each ITensor within A.
    //   Store results in mmatrix and mvector.
    int itenind = 0;
    Foreach(const ITensor& t, A.blocks())
        {

        Matrix &UU = Umatrix.at(itenind);
        Matrix &VV = Vmatrix.at(itenind);
        Vector &d =  dvector.at(itenind);

        const 
        Index &ui = t.index(uI.hasindex(t.index(1)) ? 1 : 2),
              &vi = t.index(uI.hasindex(t.index(1)) ? 2 : 1);

        Matrix M(ui.m(),vi.m());
        t.toMatrix11NoScale(ui,vi,M);

        SVD(M,UU,d,VV);

        //Store the squared singular values
        //(denmat eigenvalues) in alleig
        for(int j = 1; j <= d.Length(); ++j) 
            alleig.push_back(sqr(d(j)));

        ++itenind;
        }

    //2. Truncate eigenvalues

    //Sort all eigenvalues from smallest to largest
    //irrespective of quantum numbers
    sort(alleig.begin(),alleig.end());

    //Determine number of states to keep m
    int m = (int)alleig.size();
    Real svdtruncerr = 0;
    Real docut = -1;
    int mdisc = 0; 

    if(absoluteCutoff_)
        {
        while(m > maxm_ || ( (alleig[mdisc] < cutoff_ && m > minm_)
            && mdisc < (int)alleig.size() ) )
            {
            if(alleig[mdisc] > 0)
                svdtruncerr += alleig[mdisc];
            else
                alleig[mdisc] = 0;

            ++mdisc;
            --m;
            }
        docut = (mdisc > 0 
                ? (alleig[mdisc-1] + alleig[mdisc])*0.5 - 1E-5*alleig[mdisc-1]
                : -1);
        }
    else
	    {
	    Real scale = doRelCutoff_ ? alleig.back() : 1.0;
        while(m > maxm_ 
            || ( (svdtruncerr+alleig[mdisc] < cutoff_*scale && m > minm_)
            && mdisc < (int)alleig.size() ) )
            {
            if(alleig[mdisc] > 0)
                svdtruncerr += alleig[mdisc];
            else
                alleig[mdisc] = 0;

            ++mdisc;
            --m;
            }
        docut = (mdisc > 0 
                ? (alleig[mdisc-1] + alleig[mdisc])*0.5 - 1E-5*alleig[mdisc-1]
                : -1);
        svdtruncerr = (alleig.back() == 0 ? 0 : svdtruncerr/scale);
	    }

    if(showeigs_)
        {
        cout << endl;
        cout << "svdRank2 (IQTensor):" << endl;
        cout << format("    minm_ = %d, maxm_ = %d, cutoff_ = %.3E")
                       %minm_%maxm_%cutoff_ << endl;
        cout << format("    use_orig_m_ = %s")
                %(use_orig_m_?"true":"false")<<endl;
        cout << format("    Kept m = %d, discarded %d states in svdRank2")
                                % m % mdisc << endl;
        cout << format("    svdtruncerr = %.2E")%svdtruncerr << endl;
        cout << format("    docut = %.2E")%docut << endl;
        cout << format("    cutoff=%.2E, minm=%d, maxm=%d")
                %cutoff_%minm_%maxm_ << endl;
        cout << "    doRelCutoff is " << (doRelCutoff_ ? "true" : "false") << endl;
        cout << "    absoluteCutoff is " << (absoluteCutoff_ ? "true" : "false") << endl;
        cout << "    refNorm is " << refNorm_ << endl;

        const int s = alleig.size();
        const int max_show = 20;
        int stop = s-min(s,max_show);

        //Include refNorm_ in printed eigs as long as
        //the leading eig is within a few orders of magnitude
        //of 1.0. Otherwise just print the scaled eigs.
        Real orderMag = log(fabs(alleig.at(s-1))) + refNorm_.logNum();
        Real real_fac = 1;
        if(fabs(orderMag) < 5 && refNorm_.isFiniteReal())
            {
            real_fac = refNorm_.real();
            cout << "    Eigs: ";
            }
        else
            {
            cout << "    Eigs [omitting scale factor " << refNorm_ << "]: \n";
            if(alleig.at(s-1) > 1.e10)
                {
                Error("bad alleig");
                }
            }

        cout << "    ";
        for(int j = s-1; j >= stop; --j)
            {
            cout << format( (alleig.at(j) > 1E-6 && alleig.at(j) < 1E3) ? ("%.3f") : ("%.3E")) 
                                % (alleig[j] * real_fac);
            cout << ((j != stop) ? ", " : "\n");
            }
        cout << endl;
        } //end if(showeigs_)

    //Truncate denmat eigenvalue vectors
    //Also form new Link index with appropriate m's for each block
    vector<inqn> Liq,Riq;
    Liq.reserve(Nblock);
    Riq.reserve(Nblock);

    vector<ITensor> Ublock,
                    Vblock;
    Ublock.reserve(Nblock);
    Vblock.reserve(Nblock);

    vector<ITSparse> Dblock;
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
            //if(Global::debug1())
            //    {
            //    cout << format("    %d Keeping eig %.3E, %.3E > %.3E") % total_m % thisD(this_m) % sqr(thisD(this_m)) % docut << endl;
            //    }
            if(thisD(this_m) < 0) thisD(this_m) = 0;
            ++this_m;
            }
        --this_m; //since the loop overshoots by 1

        if(m == 0 && thisD.Length() >= 1) // zero mps, just keep one arb state
            { this_m = 1; m = 1; docut = 1; }

        if(this_m == 0) { ++itenind; continue; }

        const 
        Index &ui = t.index(uI.hasindex(t.index(1)) ? 1 : 2),
              &vi = t.index(uI.hasindex(t.index(1)) ? 2 : 1);

        Index l("l",this_m);
        Liq.push_back(inqn(l,uI.qn(ui)));

        Index r("r",this_m);
        Riq.push_back(inqn(r,vI.qn(vi)));

        Dblock.push_back(ITSparse(l,r,thisD.SubVector(1,this_m)));

        Ublock.push_back(ITensor(ui,l,UU.Columns(1,this_m)));
        Vblock.push_back(ITensor(r,vi,VV.Rows(1,this_m)));

        ++itenind;
        }

    if(Liq.size() == 0)
        throw ResultIsZero("Liq.size() == 0");

    IQIndex L("L",Liq,uI.dir()), R("R",Riq,vI.dir());

    D = IQTSparse(L,R);
    U = IQTensor(uI,conj(L));
    V = IQTensor(conj(R),vI);

    //Load blocks into D,U, and V
    for(size_t j = 0; j < Dblock.size(); ++j)
        {
        D += Dblock[j];
        U += Ublock[j];
        V += Vblock[j];
        }

    //Originally eigs were found by calling
    //toMatrix11NoScale, so put the scale back in
    D *= refNorm_;

    //Update truncerr_ and eigsKept_
    truncerr_.at(b) = svdtruncerr;

    Vector& DD = eigsKept_.at(b);
    DD.ReDimension(L.m());
    for(int i = 1; i <= L.m(); ++i) 
        DD(i) = alleig.at(alleig.size()-i);

    Global::lastd() = DD;

    /*
    {
    IQTensor Ach = U * D * V;
    Ach -= A;
    Real nor = A.norm();
    cout << "relative error in SVD is " << Ach.norm()/nor SP cutoff_ << endl;
    }
    */

    } //void SVDWorker::svdRank2


Real SVDWorker::
diag_denmat(ITensor rho, Vector& D, Index& newmid, ITensor& U)
    {
    if(isComplex(rho))
        {
        rho = realPart(rho);
        }
#ifdef DEBUG
    if(rho.r() != 2)
        {
        Print(rho.r());
        Print(rho);
        Error("Too many indices for density matrix");
        }
#endif

    Index active = rho.index(1); 
    active.noprime();

    if(!doRelCutoff_) rho.scaleTo(refNorm_);

    //Do the diagonalization
    Index ri = rho.index(1); 
    ri.noprime();
    Matrix R,UU; 
    rho.toMatrix11NoScale(ri,primed(ri),R);
    R *= -1.0; 
    EigenValues(R,D,UU); 
    D *= -1.0;

    //Include rho's scale to get the actual eigenvalues kept
    //as long as the leading eigenvalue is within a few orders
    //of magnitude of 1.0. Otherwise just report the scaled eigs.
    Real orderMag = log(fabs(D(1))) + rho.scale().logNum();
    if(fabs(orderMag) < 5 && rho.scale().isFiniteReal())
        {
        D *= rho.scale().real();
        }

    //Truncate
    int m = D.Length();
    Real svdtruncerr = 0.0;

    //Zero out any negative weight
    for(int zerom = m; zerom > 0; --zerom)
        {
        if(D(zerom) >= 0) break;
        else              D(zerom) = 0;
        }

    if(absoluteCutoff_)
        {
        while(m > maxm_ || (D(m) < cutoff_ && m > minm_))
            {
            svdtruncerr += D(m--);
            }
        }
    else
        {
        Real scale = doRelCutoff_ ? D(1) : 1.0;
        while(m > maxm_ || (svdtruncerr+D(m) < cutoff_*scale && m > minm_))
            {
            svdtruncerr += D(m--);
            }
        svdtruncerr = (D(1) == 0 ? 0 : svdtruncerr/scale);
        }

    D.ReduceDimension(m); 

    if(showeigs_)
        {
        cout << endl;
        cout << format("minm_ = %d, maxm_ = %d, cutoff_ = %.3E")
                       %minm_%maxm_%cutoff_ << endl;
        cout << format("use_orig_m_ = %s")%(use_orig_m_?"true":"false")<<endl;
        cout << format("Kept %d states in diag_denmat")% m << endl;
        cout << format("svdtruncerr = %.3E")%svdtruncerr << endl;
        //cout << "doRelCutoff is " << doRelCutoff_ << endl;
        //cout << "refNorm is " << refNorm_ << endl;
        //int stop = min(D.Length(),10);
        int stop = D.Length();
        cout << "Eigs: ";
        for(int j = 1; j <= stop; ++j)
            {
            cout << format(D(j) > 1E-3 ? ("%.3f") : ("%.3E")) % D(j);
            cout << ((j != stop) ? ", " : "\n");
            }
        cout << endl;
        }

    newmid = Index(active.rawname(),m,active.type());
    U = ITensor(active,newmid,UU.Columns(1,m));
    Global::lastd() = D;
    return svdtruncerr;
    }

Real SVDWorker::
diag_denmat(IQTensor rho, Vector& D, IQIndex& newmid, IQTensor& U)
    {
    if(isComplex(rho))
        {
        rho = realPart(rho);
        }

    if(rho.r() != 2)
        {
        rho.printIndices("rho");
        Error("Density matrix doesn't have rank 2");
        }

    vector<Matrix> mmatrix(rho.iten_size());
    vector<Vector> mvector(rho.iten_size());
    vector<Real> alleig;
    alleig.reserve(rho.index(1).m());

    if(rho.index(1).m() == 0)
        throw ResultIsZero("rho.index(1).m()");
    if(rho.iten_empty())
        throw ResultIsZero("rho.iten_empty()");

    if(doRelCutoff_)
        {
        //DO_IF_DEBUG(cout << "Doing relative cutoff\n";)
        Real maxLogNum = -200;
        rho.scaleOutNorm();
        Foreach(const ITensor& t, rho.blocks())
            {
            maxLogNum = max(maxLogNum,t.scale().logNum());
            }
        refNorm_ = LogNumber(maxLogNum,1);
        }
    //DO_IF_DEBUG(cout << "refNorm = " << refNorm_ << endl; )
    //cout << "WARNING - SETTING REFNORM TO 10\n";
    //refNorm_ = LogNumber(log(10),-1);
    //else DO_IF_DEBUG(cout << "Not doing relative cutoff\n";);

    //cerr << boost::format("refNorm = %.1E (lognum = %f, sign = %d)\n\n")
    //%Real(refNorm)%refNorm.logNum()%refNorm.sign();

    rho.scaleTo(refNorm_);

    //1. Diagonalize each ITensor within rho.
    //   Store results in mmatrix and mvector.
    int itenind = 0;
    Foreach(const ITensor& t, rho.blocks())
        {
        if(!t.index(1).noprimeEquals(t.index(2)))
            { 
            Print(rho); 
            Print(t); 
            Error("Non-symmetric ITensor in density matrix, perhaps QNs not conserved?");
            }

        Matrix &UU = mmatrix.at(itenind);
        Vector &d =  mvector.at(itenind);

        //Diag ITensors within rho
        int n = t.index(1).m();
        Matrix M(n,n);
        t.toMatrix11NoScale(t.index(1),t.index(2),M);

        M *= -1;
        EigenValues(M,d,UU);
        d *= -1;

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
            cerr << boost::format("\ndiff=%.2E\n")%Norm(Diff.TreatAsVector());
            Print(UU.t()*UU);
            Error("UU not unitary in diag_denmat");
            }
        
#endif //STRONG_DEBUG
        }

    //2. Truncate eigenvalues

    //Sort all eigenvalues from smallest to largest
    //irrespective of quantum numbers
    sort(alleig.begin(),alleig.end());

    //Determine number of states to keep m
    int m = (int)alleig.size();
    Real svdtruncerr = 0;
    Real docut = -1;
    int mdisc = 0; 

    if(absoluteCutoff_)
        {
        while(m > maxm_ || ( (alleig[mdisc] < cutoff_ && m > minm_)
            && mdisc < (int)alleig.size() ) )
            {
            if(alleig[mdisc] > 0)
                svdtruncerr += alleig[mdisc];
            else
                alleig[mdisc] = 0;

            ++mdisc;
            --m;
            }
        docut = (mdisc > 0 
                ? (alleig[mdisc-1] + alleig[mdisc])*0.5 - 1E-5*alleig[mdisc-1]
                : -1);
        }
    else
	    {
	    Real scale = doRelCutoff_ ? alleig.back() : 1.0;
        while(m > maxm_ 
            || ( (svdtruncerr+alleig[mdisc] < cutoff_*scale && m > minm_)
            && mdisc < (int)alleig.size() ) )
            {
            if(alleig[mdisc] > 0)
                svdtruncerr += alleig[mdisc];
            else
                alleig[mdisc] = 0;

            ++mdisc;
            --m;
            }
        docut = (mdisc > 0 
                ? (alleig[mdisc-1] + alleig[mdisc])*0.5 - 1E-5*alleig[mdisc-1]
                : -1);
        svdtruncerr = (alleig.back() == 0 ? 0 : svdtruncerr/scale);
	    }

    if(showeigs_)
        {
        cout << endl;
        cout << boost::format("use_orig_m_ = %s")
                %(use_orig_m_?"true":"false")<<endl;
        cout << boost::format("Kept %d, discarded %d states in diag_denmat line 721")
                                % m % mdisc << endl;
        cout << boost::format("svdtruncerr = %.2E")%svdtruncerr << endl;
        cout << boost::format("docut = %.2E")%docut << endl;
        cout << boost::format("cutoff=%.2E, minm=%d, maxm=%d")
                %cutoff_%minm_%maxm_ << endl;
        cout << "doRelCutoff is " << (doRelCutoff_ ? "true" : "false") << endl;
        cout << "absoluteCutoff is " << (absoluteCutoff_ ? "true" : "false") << endl;
        cout << "refNorm is " << refNorm_ << endl;
        int s = alleig.size();
        const int max_show = 20;
        int stop = s-min(s,max_show);
        cout << "Eigs: ";
        for(int j = s-1; j >= stop; --j)
            {
            cout << format(alleig[j] > 1E-3 ? ("%.3f") : ("%.3E")) 
                           % alleig[j];
            cout << ((j != stop) ? ", " : "\n");
            }
        cout << endl;
        }

    assert(m <= maxm_); 
    assert(m < 20000);

    IQIndex active = (rho.index(1).primeLevel() == 0 ? rho.index(1)
                                                     : rho.index(2));

    //
    //Truncate eigenvalues and eigenvectors of rho
    //

    //Build blocks for unitary diagonalizing rho
    vector<ITensor> blocks; 
    blocks.reserve(rho.iten_size());

    //Also form new Link IQIndex with appropriate m's for each block
    vector<inqn> iq; 
    iq.reserve(rho.iten_size());

    itenind = 0;
    Foreach(const ITensor& t, rho.blocks())
        {
        Vector& thisD = mvector.at(itenind);
        Matrix& thisU = mmatrix.at(itenind);

        int this_m = 1;
        while(this_m <= thisD.Length() && thisD(this_m) > docut) 
            {
            if(thisD(this_m) < 0) thisD(this_m) = 0;
            ++this_m;
            }
        --this_m; //since the loop overshoots by 1

        if(m == 0 && thisD.Length() >= 1) // zero mps, just keep one arb state
            { this_m = 1; m = 1; docut = 1; }

        thisD.ReduceDimension(this_m);

        if(this_m == 0) { ++itenind; continue; }

        Index nm("qlink",this_m);
        Index act = deprimed(t.index(1));
        iq.push_back(inqn(nm,active.qn(act)));

        MatrixRef Utrunc = thisU.Columns(1,this_m);

        ITensor block(act,nm);
        block.fromMatrix11(act,nm,Utrunc);
        blocks.push_back(block);

        ++itenind;
        }

    if(iq.size() == 0)
        {
        Print(m);
        Print(docut);
        throw ResultIsZero("iq.size() == 0");
        }

    newmid = IQIndex("qlink",iq,active.dir()*Switch);

    U = IQTensor(active,newmid);
    Foreach(const ITensor& block, blocks)
        U += block;

    D.ReDimension(newmid.m());
    for(int i = 1; i <= newmid.m(); ++i) 
        D(i) = alleig.at(alleig.size()-i);

    Global::lastd() = D;
    //Include refNorm_ to get the actual eigenvalues kept
    //as long as the leading eigenvalue is within a few orders
    //of magnitude of 1.0. Otherwise just report the scaled eigs.
    Real orderMag = log(fabs(D(1))) + refNorm_.logNum();
    if(fabs(orderMag) < 5 && refNorm_.isFiniteReal())
        {
        Global::lastd() *= refNorm_.real();
        }

    return svdtruncerr;

    } //void SVDWorker::diag_denmat


int SVDWorker::
maxEigsKept() const
    {
    int res = -1;
    Foreach(const Vector& eigs,eigsKept_)
        res = max(res,eigs.Length());
    return res;
    }

Real SVDWorker::
maxTruncerr() const
    {
    Real res = -1;
    Foreach(const Real& te,truncerr_)
        res = max(res,te);
    return res;
    }


void SVDWorker::
read(std::istream& s)
    {
    s.read((char*) &N_,sizeof(N_));
    truncerr_.resize(N_+1);
    for(int j = 1; j <= N_; ++j)
        s.read((char*)&truncerr_[j],sizeof(truncerr_[j]));
    s.read((char*)&cutoff_,sizeof(cutoff_));
    s.read((char*)&minm_,sizeof(minm_));
    s.read((char*)&maxm_,sizeof(maxm_));
    s.read((char*)&use_orig_m_,sizeof(use_orig_m_));
    s.read((char*)&showeigs_,sizeof(showeigs_));
    s.read((char*)&doRelCutoff_,sizeof(doRelCutoff_));
    s.read((char*)&absoluteCutoff_,sizeof(absoluteCutoff_));
    s.read((char*)&refNorm_,sizeof(refNorm_));
    eigsKept_.resize(N_+1);
    for(int j = 1; j <= N_; ++j)
        readVec(s,eigsKept_.at(j));
    }

void SVDWorker::
write(std::ostream& s) const
    {
    s.write((char*) &N_,sizeof(N_));
    for(int j = 1; j <= N_; ++j)
        s.write((char*)&truncerr_[j],sizeof(truncerr_[j]));
    s.write((char*)&cutoff_,sizeof(cutoff_));
    s.write((char*)&minm_,sizeof(minm_));
    s.write((char*)&maxm_,sizeof(maxm_));
    s.write((char*)&use_orig_m_,sizeof(use_orig_m_));
    s.write((char*)&showeigs_,sizeof(showeigs_));
    s.write((char*)&doRelCutoff_,sizeof(doRelCutoff_));
    s.write((char*)&absoluteCutoff_,sizeof(absoluteCutoff_));
    s.write((char*)&refNorm_,sizeof(refNorm_));
    for(int j = 1; j <= N_; ++j)
        writeVec(s,eigsKept_[j]);
    }



ITensor SVDWorker::
pseudoInverse(const ITensor& C, Real cutoff)
    {
    if(C.r() != 1)
        {
        Print(C);
        Error("pseudoInverse only defined for rank 1 ITensors");
        }
    ITensor res(C);
    res.pseudoInvert(cutoff);
    return res;
    }

IQTensor SVDWorker::
pseudoInverse(const IQTensor& C, Real cutoff)
    {
    if(C.r() != 1)
        {
        C.printIndices("C");
        Error("pseudoInverse only defined for rank 1 ITensors");
        }
    IQTensor res(C.index(1));
    Foreach(const ITensor& t, C.blocks())
        { res += pseudoInverse(t,cutoff); }
    return res;
    }


/*
Real SVDWorker::
diag_denmat(const ITensor& rho, Vector& D, Index& newmid, ITensor& C, ITensor& U)
    {
    Real svdtruncerr = diag_denmat(rho,D,newmid,U);
    C = ITensor(newmid,sqrt(D));
    return svdtruncerr;
    }
*/

/*
Real SVDWorker::
diag_denmat_complex(const ITensor& rho, Vector& D, Index& newmid, ITensor& U)
    {
    // Need to fix this to put in Hermitian case!
    ITensor rhore(realPart(rho)),
            rhoim(imagPart(rho));
    return diag_denmat(rhore,D,newmid,U);
    }
    */


//
// Helper method for SVDWorker::diag_denmat(const IQTensor& rho,...)
// Diagonalizes and truncates the density matrix but doesn't create 
// any IQTensors such as U
//
/*
void SVDWorker::
diag_and_truncate(const IQTensor& rho, vector<Matrix>& mmatrix, 
                  vector<Vector>& mvector,
                  vector<Real>& alleig, Real& svdtruncerr, IQIndex& newmid)
    {
    if(rho.r() != 2)
        {
        rho.printIndices("rho");
        Error("Density matrix doesn't have rank 2");
        }

    mmatrix = vector<Matrix>(rho.iten_size());
    mvector = vector<Vector>(rho.iten_size());
    alleig.clear();
    alleig.reserve(rho.index(1).m());
    if(rho.index(1).m() == 0)
	throw ResultIsZero("rho.index(1).m()");
    if(rho.iten_size() == 0)
	throw ResultIsZero("rho.iten_size() == 0");

    if(doRelCutoff_)
        {
        //DO_IF_DEBUG(cout << "Doing relative cutoff\n";)
        Real maxLogNum = -200;
        Foreach(const ITensor& t, rho.blocks())
            maxLogNum = max(maxLogNum,t.scale().logNum());
        refNorm_ = LogNumber(maxLogNum,1);
        }
    //DO_IF_DEBUG(cout << "refNorm = " << refNorm_ << endl; )
    //cout << "WARNING - SETTING REFNORM TO 10\n";
    //refNorm_ = LogNumber(log(10),-1);
    //else DO_IF_DEBUG(cout << "Not doing relative cutoff\n";);

    //cerr << boost::format("refNorm = %.1E (lognum = %f, sign = %d)\n\n")
    //%Real(refNorm)%refNorm.logNum()%refNorm.sign();

    //1. Diagonalize each ITensor within rho.
    //   Store results in mmatrix and mvector.
    int itenind = 0;
    for(IQTensor::const_iten_it it = rho.const_iten_begin(); 
        it != rho.const_iten_end(); ++it)
        {
        const ITensor& t = *it;
        if(!t.index(1).noprime_equals(t.index(2)))
            { 
            Print(rho); 
            Print(t); 
            Error("Non-symmetric ITensor in density matrix, perhaps QNs not conserved?");
            }

        t.scaleTo(refNorm_);

        Matrix &UU = mmatrix.at(itenind);
        Vector &d =  mvector.at(itenind);

        //Diag ITensors within rho
        int n = t.index(1).m();
        Matrix M(n,n);
        t.toMatrix11NoScale(t.index(1),t.index(2),M);

        M *= -1;
        EigenValues(M,d,UU);
        d *= -1;

        d *= refNorm_.real();

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
            cerr << boost::format("\ndiff=%.2E\n")%Norm(Diff.TreatAsVector());
            Print(UU.t()*UU);
            Error("UU not unitary in diag_denmat");
            }
        
        if(fabs(d.sumels() + Trace(M)*refNorm_.real())/(fabs(d.sumels())+fabs(Trace(M)*refNorm_.real())) 
            > 1E-5)
            {
            cerr << boost::format("d.sumels() = %.10f, Trace(M)*refNorm_.real() = %.10f\n")
                     % d.sumels()        % (Trace(M)*refNorm_.real());
            Error("Total eigs != trace");
            }

#endif //STRONG_DEBUG
        }

    //2. Truncate eigenvalues

    //Sort all eigenvalues from smallest to largest
    //irrespective of quantum numbers
    sort(alleig.begin(),alleig.end());

    //Determine number of states to keep m
    int m = (int)alleig.size();
    svdtruncerr = 0;
    Real docut = -1;
    int mdisc = 0; 

    if(absoluteCutoff_)
        {
        while(m > maxm_ || (alleig[mdisc] < cutoff_ && m > minm_)
            && mdisc < (int)alleig.size())
            {
            if(alleig[mdisc] > 0)
                svdtruncerr += alleig[mdisc];
            else
                alleig[mdisc] = 0;

            ++mdisc;
            --m;
            }
        docut = (mdisc > 0 ? (alleig[mdisc-1] + alleig[mdisc])*0.5 : -1) + 1E-40;
        }
    else
	    {
	    Real scale = doRelCutoff_ ? alleig.back() : 1.0;
        while(m > maxm_ 
            || (svdtruncerr+alleig[mdisc] < cutoff_*scale && m > minm_)
            && mdisc < (int)alleig.size())
            {
            if(alleig[mdisc] > 0)
                svdtruncerr += alleig[mdisc];
            else
                alleig[mdisc] = 0;

            ++mdisc;
            --m;
            }
        docut = (mdisc > 0 
                ? (alleig[mdisc-1] + alleig[mdisc])*0.5 
                : -1) + 1E-40;
        svdtruncerr = (alleig.back() == 0 ? 0 : svdtruncerr/scale);
	    }

    if(showeigs_)
        {
        cout << endl;
        cout << boost::format("use_orig_m_ = %s")
                %(use_orig_m_?"true":"false")<<endl;
        cout << boost::format("Kept %d, discarded %d states in diag_denmat")
                                % m % mdisc << endl;
        cout << boost::format("svdtruncerr = %.2E")%svdtruncerr << endl;
        cout << boost::format("docut = %.2E")%docut << endl;
        cout << boost::format("cutoff=%.2E, minm=%d, maxm=%d")
                %cutoff_%minm_%maxm_ << endl;
        cout << "doRelCutoff is " << (doRelCutoff_ ? "true" : "false") << endl;
        cout << "absoluteCutoff is " << (absoluteCutoff_ ? "true" : "false") << endl;
        cout << "refNorm is " << refNorm_ << endl;
        int s = alleig.size();
        const int max_show = 20;
        int stop = s-min(s,max_show);
        cout << "Eigs: ";
        for(int j = s-1; j >= stop; --j)
            {
            cout << format(alleig[j] > 1E-3 ? ("%.3f") : ("%.3E")) 
                                % alleig[j];
            cout << ((j != stop) ? ", " : endl);
            }
        }

    assert(m <= maxm_); 
    assert(m < 20000);

    IQIndex active = (rho.index(1).primeLevel() == 0 ? rho.index(1)
                                                     : rho.index(2));

    //Truncate denmat eigenvalue vectors
    //Also form new Link index with appropriate m's for each block
    vector<inqn> iq; iq.reserve(rho.iten_size());
    itenind = 0;
    for(IQTensor::const_iten_it it = rho.const_iten_begin(); 
        it != rho.const_iten_end(); ++it)
        {
        const ITensor& t = *it;
        Vector& thisD = mvector.at(itenind);

        int this_m = 1;
        while(this_m <= thisD.Length() && thisD(this_m) > docut) 
            {
            if(thisD(this_m) < 0) thisD(this_m) = 0;
            ++this_m;
            }
        --this_m; //since the loop overshoots by 1

        if(m == 0 && thisD.Length() >= 1) // zero mps, just keep one arb state
            { this_m = 1; m = 1; docut = 1; }

        thisD.ReduceDimension(this_m);

        if(this_m == 0) { ++itenind; continue; }

        Index nm("qlink",this_m);
        Index act = deprimed(t.index(1));
        iq.push_back(inqn(nm,active.qn(act)));

        ++itenind;
        }
    if(iq.size() == 0)
        {
        Print(m);
        Print(docut);
        throw ResultIsZero("iq.size() == 0");
        }
    newmid = IQIndex("qlink",iq,active.dir()*Switch);
    } //void SVDWorker::diag_and_truncate

void SVDWorker::
buildUnitary(const IQTensor& rho, const vector<Matrix>& mmatrix, 
             const vector<Vector>& mvector,
             const IQIndex& newmid, IQTensor& U)
    {
    IQIndex active = (rho.index(1).primeLevel() == 0 ? rho.index(1)
                                                     : rho.index(2));

    // Construct orthogonalized IQTensor U
    vector<ITensor> terms; terms.reserve(rho.iten_size());
    int itenind = 0, kept_block = 0;
    int m = newmid.m();
    for(IQTensor::const_iten_it it = rho.const_iten_begin(); 
        it != rho.const_iten_end(); ++it)
        {
        const Vector& thisD = mvector.at(itenind);
        int this_m = thisD.Length();

        if(m == 0 && thisD.Length() >= 1) // zero mps, just keep one arb state
            { this_m = 1; m = 1; }

        if(this_m == 0) { ++itenind; continue; }

        const Index& nm = newmid.index(++kept_block);
        Index act = deprimed(it->index(1));

#ifdef DEBUG
        if(nm.m() != this_m)
            {
            Print(nm.m());
            Print(this_m);
            Error("Mismatched m");
            }
#endif

        Matrix Utrunc = mmatrix[itenind].Columns(1,this_m);

        ITensor term(act,nm); 
        term.fromMatrix11(act,nm,Utrunc); 
        terms.push_back(term);

        ++itenind;
        }
    U = IQTensor(active,newmid);
    for(size_t j = 0; j < terms.size(); ++j)
        U += terms[j];
    }

void SVDWorker::
buildCenter(const IQTensor& rho, const vector<Matrix>& mmatrix, 
            const vector<Vector>& mvector,
            const IQIndex& newmid, IQTensor& C)
    {
    vector<ITensor> terms; terms.reserve(rho.iten_size());
    int itenind = 0, kept_block = 0;
    int m = newmid.m();
    for(IQTensor::const_iten_it it = rho.const_iten_begin(); 
        it != rho.const_iten_end(); ++it)
        {
        const Vector& thisD = mvector.at(itenind);
        int this_m = thisD.Length();

        if(m == 0 && thisD.Length() >= 1) // zero mps, just keep one arb state
            { this_m = 1; m = 1; }

        if(this_m == 0) { ++itenind; continue; }

        const Index& nm = newmid.index(++kept_block);

        ITensor term(nm,sqrt(thisD)); 
        terms.push_back(term);

        ++itenind;
        }
    C = IQTensor(conj(newmid));
    for(size_t j = 0; j < terms.size(); ++j)
        C += terms[j];
    }


Real SVDWorker::diag_denmat(const IQTensor& rho, Vector& D, IQIndex& newmid, 
                            IQTensor& U)
    {
    vector<Matrix> mmatrix;
    vector<Vector> mvector;
    vector<Real> alleig;
    Real svdtruncerr = 0;

    diag_and_truncate(rho,mmatrix,mvector,alleig,svdtruncerr,newmid);
    
    buildUnitary(rho,mmatrix,mvector,newmid,U);

    D.ReDimension(newmid.m());
    for(int i = 1; i <= newmid.m(); ++i) 
        D(i) = alleig.at(alleig.size()-i);
    Global::lastd() = D;
    return svdtruncerr;
    } //Real SVDWorker::diag_denmat

Real SVDWorker::diag_denmat(const IQTensor& rho, Vector& D, IQIndex& newmid, IQTensor& C, IQTensor& U)
    {
    vector<Matrix> mmatrix;
    vector<Vector> mvector;
    vector<Real> alleig;
    Real svdtruncerr = 0;

    diag_and_truncate(rho,mmatrix,mvector,alleig,svdtruncerr,newmid);
    
    buildUnitary(rho,mmatrix,mvector,newmid,U);

    buildCenter(rho,mmatrix,mvector,newmid,C);

    D.ReDimension(newmid.m());
    for(int i = 1; i <= newmid.m(); ++i) 
        D(i) = alleig.at(alleig.size()-i);
    Global::lastd() = D;
    return svdtruncerr;
    } //Real SVDWorker::diag_denmat

    */

/*
Real SVDWorker::diag_denmat_complex(const IQTensor& rho, Vector& D, IQIndex& newmid, IQTensor& U)
    {
    bool docomplex = false;
    IQIndex active;
    Global::printdat() = true;
    if(rho.r() == 3) 
	{
	docomplex = true;
	for(int i = 1; i <= 3; i++)
	    if(rho.index(i).dir() == Out)
		active = rho.index(i);
	}
    else
	active = rho.finddir(Out);

    assert(active.primeLevel() == 0);

    vector<Matrix> mmatrixre(rho.iten_size());
    vector<Matrix> mmatrixim(rho.iten_size());
    vector<Vector> mvector(rho.iten_size());
    vector<Real> alleig;

    if(doRelCutoff_)
	{
        Real maxLogNum = -200;
        Foreach(const ITensor& t, rho.blocks())
	    {
	    t.scaleOutNorm();
	    maxLogNum = max(maxLogNum,t.scale().logNum());
	    }
        refNorm_ = LogNumber(maxLogNum,1);
	}
    //1. Diagonalize each ITensor within rho
    int itenind = 0;
    for(IQTensor::const_iten_it it = rho.const_iten_begin(); it != rho.const_iten_end(); ++it)
        {
        const ITensor& t = *it;
        t.scaleTo(refNorm_);

        Matrix &UUre = GET(mmatrixre,itenind);
        Matrix &UUim = GET(mmatrixim,itenind);
        Vector &d = GET(mvector,itenind);

        //Diag ITensors within rho
	int ii = 1, jj = 2;
	if(docomplex)
	    {
	    if(t.index(1) == IQIndex::IndReIm())
		ii = 2, jj = 3;
	    else if(t.index(2) == IQIndex::IndReIm())
		ii = 1, jj = 3;
	    else if(t.index(3) == IQIndex::IndReIm())
		ii = 1, jj = 2;
	    else 
		Error("bad IndReIm");
	    }
	//cout << "nontrivial indices are " << t.index(ii) SP t.index(jj) << endl;
	//cout << "t is " << t << endl;
	//cout << "refNorm_ is " << refNorm_ << endl;
        int n = t.index(ii).m();
        Matrix Mre(n,n), Mim(n,n);

	if(docomplex)
	    {
        ITensor tre(realPart(t)),
                tim(imagPart(t));
	    tre.scaleTo(refNorm_);
	    tim.scaleTo(refNorm_);
	    tre.toMatrix11NoScale(t.index(ii),t.index(jj),Mre);
	    tim.toMatrix11NoScale(t.index(ii),t.index(jj),Mim);
	    Mre *= -1;
	    Mim *= -1;
	    HermitianEigenvalues(Mre,Mim,d,UUre,UUim);
	    }
	else
	    {
	    t.toMatrix11NoScale(t.index(ii),t.index(jj),Mre);
	    Mre *= -1;
	    EigenValues(Mre,d,UUre);
	    }
	d *= -1;
        for(int j = 1; j <= n; ++j) 
            { alleig.push_back(d(j)); }
        ++itenind;
        }

    //2. Truncate eigenvalues

    //Sort all eigenvalues from smallest to largest
    //irrespective of quantum numbers
    sort(alleig.begin(),alleig.end());

    //Truncate
    Real docut = -1;
    Real svdtruncerr = 0;
    Real e1 = max(alleig.back(),1.0e-60);
    int mdisc = 0, m = (int)alleig.size();
    if(absoluteCutoff_)
	{
	//Sort all eigenvalues from largest to smallest
	//irrespective of quantum numbers
	reverse(alleig.begin(),alleig.end());
	m = minm_;
	while(m < maxm_ && m < (int)alleig.size() && alleig[m-1] > cutoff_ ) 
	    svdtruncerr += alleig[m++ - 1];
	reverse(alleig.begin(),alleig.end());
	mdisc = (int)alleig.size() - m;
	docut = (mdisc > 0 ?  (alleig.at(mdisc-1) + alleig[mdisc])*0.5 : -1);
	}
    else
	if(m > minm_)
	    {
	    Real sca = doRelCutoff_ ? e1 : 1.0;
	    for(; mdisc < (int)alleig.size(); mdisc++, m--)
		{
		if(((svdtruncerr += GET(alleig,mdisc)/sca) > cutoff_ && m <= maxm_) 
			   || m <= minm_)
		    { 
		    docut = (mdisc > 0 ? (alleig[mdisc-1] + alleig[mdisc])*0.5 : -1);
		    //Overshot by one, correct truncerr
		    svdtruncerr -= alleig[mdisc]/sca;
		    break; 
		    }
		}
	    }
    if(showeigs_)
	{
        cout << endl;
        cout << boost::format("use_orig_m_ = %s")%(use_orig_m_?"true":"false")<<endl;
        cout << boost::format("Kept %d, discarded %d states in diag_denmat_complex line 1310")
                                     % m % mdisc << endl;
        cout << boost::format("svdtruncerr = %.2E")%svdtruncerr << endl;
        cout << boost::format("docut = %.2E")%docut << endl;
        cout << boost::format("cutoff=%.2E, minm=%d, maxm=%d")%cutoff_%minm_%maxm_ << endl;
        cout << "doRelCutoff is " << (doRelCutoff_ ? "true" : "false") << endl;
        cout << "absoluteCutoff is " << (absoluteCutoff_ ? "true" : "false") << endl;
        cout << "refNorm is " << refNorm_ << endl;
        int s = alleig.size();
        const int max_show = 20;
        int stop = s-min(s,max_show);
        cout << "Eigs: ";
        for(int j = s-1; j >= stop; --j)
	    {
            cout << format(alleig[j] > 1E-3 ? ("%.3f") : ("%.3E")) 
                                % alleig[j];
            cout << ((j != stop) ? ", " : "\n");
	    }
        cout << endl;
	}

    assert(m <= maxm_); 
    assert(m < 20000);

    //3. Construct orthogonalized IQTensor U
    vector<ITensor> terms; terms.reserve(rho.iten_size());
    vector<inqn> iq; iq.reserve(rho.iten_size());
    itenind = 0;
    for(IQTensor::const_iten_it it = rho.const_iten_begin(); it != rho.const_iten_end(); ++it)
	{
        const ITensor& t = *it;
        const Vector& thisD = GET(mvector,itenind);

        int this_m = 1;
        for(; this_m <= thisD.Length(); ++this_m)
            if(thisD(this_m) < docut) 
		break;
        --this_m; //since for loop overshoots by 1

        if(m == 0 && thisD.Length() >= 1) // zero mps, just keep one arb state
            { this_m = 1; m = 1; docut = 1; }

        if(this_m == 0) 
	    { 
	    ++itenind; 
	    continue; 
	    }

        Index nm("qlink",this_m);
        Index act = deprimed(t.index(1));
	if(docomplex && act == Index::IndReIm())
	    act = deprimed(t.index(2));
        iq.push_back(inqn(nm,active.qn(act)));

        Matrix Utruncre = GET(mmatrixre,itenind).Columns(1,this_m);
        Matrix Utruncim;
	if(docomplex)
	    Utruncim = GET(mmatrixim,itenind).Columns(1,this_m);

        ITensor termre(act,nm),termim(act,nm); 
        termre.fromMatrix11(act,nm,Utruncre); 
	if(docomplex)
	    {
	    termim.fromMatrix11(act,nm,Utruncim); 
	    termre += ITensor::Complex_i() * termim;
	    }
        terms.push_back(termre);
        ++itenind;
	}
    newmid = IQIndex("qlink",iq,In);
    U = IQTensor(active,newmid);
    if(docomplex)
	U *= IQTensor::Complex_1();
    Foreach(const ITensor& t, terms) 
	U += t;
    D.ReDimension(m);
    for(int i = 1; i <= m; ++i) 
        D(i) = GET(alleig,alleig.size()-i);
    Global::lastd() = D;
    return svdtruncerr;
    } //Real SVDWorker::diag_denmat
*/
