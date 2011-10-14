#include "svdworker.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::pair;

Real SVDWorker::diag_denmat(const ITensor& rho, Vector& D, ITensor& U)
    {
    assert(rho.r() == 2);
    Index active = rho.index(1); active.noprime();

    if(!doRelCutoff_) rho.scaleTo(refNorm_);

    //Do the diagonalization
    Index ri = rho.index(1); ri.noprime();
    Matrix R,UU; rho.toMatrix11NoScale(ri,ri.primed(),R);
    R *= -1.0; EigenValues(R,D,UU); D *= -1.0;

    //Truncate
    Real svdtruncerr = 0.0;
    int mp = D.Length();
    Real sca = doRelCutoff_ ? D(1) : 1.0;
    if(absoluteCutoff_)
	{
	mp = minm_;
	while(mp < maxm_ && D(mp) < cutoff_ ) 
	    svdtruncerr += D(mp++);
	}
    else
	while(mp > maxm_ || (svdtruncerr+D(mp) < cutoff_*sca && mp > minm_)) 
	    svdtruncerr += D(mp--);
    if(!absoluteCutoff_)
	svdtruncerr = (D(1) == 0 ? 0 : svdtruncerr/sca);
    D.ReduceDimension(mp); 
    if(showeigs_)
	{
	cout << endl;
	cout << boost::format("truncate_ = %s")%(truncate_?"true":"false")<<endl;
	cout << boost::format("Kept %d states in diag_denmat\n")% mp;
	cout << boost::format("svdtruncerr = %.2E\n")%svdtruncerr;
	//cout << "doRelCutoff is " << doRelCutoff_ << endl;
	//cout << "refNorm is " << refNorm_ << endl;
	int stop = min(D.Length(),10);
	cout << "Eigs: ";
	for(int j = 1; j <= stop; ++j)
	    {
	    cout << boost::format(D(j) > 1E-3 ? ("%.3f") : ("%.3E")) % D(j);
	    cout << ((j != stop) ? ", " : "\n");
	    }
	}
    Index newmid(active.rawname(),mp,active.type());
    U = ITensor(active,newmid,UU.Columns(1,mp));
    lastd = D;
    return svdtruncerr;
    }

Real SVDWorker::diag_denmat_complex(const ITensor& rho, Vector& D, ITensor& U)
    {
    ITensor rhore,rhoim;
    rho.SplitReIm(rhore,rhoim);		// Need to fix this to put in Hermitian case!
    return diag_denmat(rhore,D,U);
    }

Real SVDWorker::diag_denmat(const IQTensor& rho, Vector& D, IQTensor& U)
    {
    assert(rho.r() == 2);
    IQIndex active = rho.finddir(Out);
    assert(active.primeLevel() == 0);

    vector<Matrix> mmatrix(rho.iten_size());
    vector<Vector> mvector(rho.iten_size());
    vector<Real> alleig;

    if(doRelCutoff_)
	{
        //DO_IF_DEBUG(cout << "Doing relative cutoff\n";)
        Real maxLogNum = -200;
        foreach(const ITensor& t, rho.itensors())
	    maxLogNum = max(maxLogNum,t.scale().logNum());
        refNorm_ = LogNumber(maxLogNum,1);
        //DO_IF_DEBUG(cout << "refNorm = " << refNorm << endl; )
	}
    //else DO_IF_DEBUG(cout << "Not doing relative cutoff\n";);

    //cerr << boost::format("refNorm = %.1E (lognum = %f, sign = %d)\n\n")
    //%Real(refNorm)%refNorm.logNum()%refNorm.sign();


    //1. Diagonalize each ITensor within rho
    int itenind = 0;
    for(IQTensor::const_iten_it it = rho.const_iten_begin(); it != rho.const_iten_end(); ++it)
	{
        const ITensor& t = *it;
        if(!t.index(1).noprime_equals(t.index(2)))
	    { 
	    Print(rho); 
	    Print(t); 
	    Error("Non-symmetric ITensor in density matrix, perhaps QNs not conserved?");
	    }

        t.scaleTo(refNorm_);

        Matrix &UU = GET(mmatrix,itenind);
        Vector &d = GET(mvector,itenind);

        //Diag ITensors within rho
        int n = t.index(1).m();
        Matrix M(n,n);
        t.toMatrix11NoScale(t.index(1),t.index(2),M);

#ifdef STRONG_DEBUG
        for(int r = 1; r <= n; ++r)
	    for(int c = r+1; c <= n; ++c)
		{
		if(fabs(M(r,c)-M(c,r)) > 1E-15)
		    {
		    Print(M);
		    Error("M not symmetric in diag_denmat");
		    }
		}
#endif //STRONG_DEBUG

        M *= -1;
        EigenValues(M,d,UU);
        d *= -1;

#ifdef STRONG_DEBUG
        Matrix Id(UU.Nrows(),UU.Nrows()); Id = 1;
        Matrix Diff = Id-(UU.t()*UU);
        if(Norm(Diff.TreatAsVector()) > 1E-12)
	    {
	    cerr << boost::format("\ndiff=%.2E\n")%Norm(Diff.TreatAsVector());
	    Print(UU.t()*UU);
	    Error("UU not unitary in diag_denmat");
	    }
        
        if(fabs(d.sumels() + Trace(M))/(fabs(d.sumels())+fabs(Trace(M))) > 1E-5)
	    {
	    cerr << boost::format("d.sumels() = %.10f, Trace(M) = %.10f\n")
				 % d.sumels()        % Trace(M);
	    Error("Total eigs != trace");
	    }

        /*
        Matrix DD(n,n); DD.TreatAsVector() = 0;
        for(int j = 1; j <= n; ++j) DD(j,j) = -d(j);
        Matrix nM = UU*DD*UU.t();
        for(int r = 1; r <= n; ++r)
        for(int c = r+1; c <= n; ++c)
        {
            if(fabs(M(r,c)) < 1E-16) continue;
            if(fabs(nM(r,c)-M(r,c))/(fabs(nM(r,c))+fabs(M(r,c))) > 1E-3)
            {
                Print(M);
                Print(nM);
                cerr << boost::format("nM(r,c)=%.2E\n")%nM(r,c);
                cerr << boost::format(" M(r,c)=%.2E\n")%M(r,c);
                Error("Inaccurate diag");
            }
        }
        */
#endif //STRONG_DEBUG

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
        docut = (mdisc > 0 ?  (alleig[mdisc-1] + alleig[mdisc])*0.5 : -1);
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
		    docut = (mdisc > 0 ?  (alleig[mdisc-1] + alleig[mdisc])*0.5 : -1);
		    //Overshot by one, correct truncerr
		    svdtruncerr -= alleig[mdisc]/sca;
		    break; 
		    }
		}
	    }
    if(showeigs_)
	{
        cout << endl;
        cout << boost::format("truncate_ = %s")%(truncate_?"true":"false")<<endl;
        cout << boost::format("Kept %d, discarded %d states in diag_denmat")
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
            cout << boost::format(alleig[j] > 1E-3 ? ("%.3f") : ("%.3E")) 
                                % alleig[j];
            cout << ((j != stop) ? ", " : "\n");
	    }
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

        if(this_m == 0) { ++itenind; continue; }

        Index nm("qlink",this_m);
        Index act = t.index(1).deprimed();
        iq.push_back(inqn(nm,active.qn(act)));

        Matrix Utrunc = GET(mmatrix,itenind).Columns(1,this_m);

        ITensor term(act,nm); 
        term.fromMatrix11(act,nm,Utrunc); 
        terms.push_back(term);

        ++itenind;
	}
    IQIndex newmid("qlink",iq,In);
    U = IQTensor(active,newmid);
    foreach(const ITensor& t, terms) 
	U += t;
    D.ReDimension(m);
    for(int i = 1; i <= m; ++i) 
        D(i) = GET(alleig,alleig.size()-i);
    lastd = D;
    return svdtruncerr;
    } //Real SVDWorker::diag_denmat

Real SVDWorker::diag_denmat_complex(const IQTensor& rho, Vector& D, IQTensor& U)
    {
    bool docomplex = false;
    IQIndex active;
    printdat = true;
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
        foreach(const ITensor& t, rho.itensors())
	    maxLogNum = max(maxLogNum,t.scale().logNum());
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
	ITensor tre,tim;
	if(docomplex)
	    {
	    t.SplitReIm(tre,tim);
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
        cout << boost::format("truncate_ = %s")%(truncate_?"true":"false")<<endl;
        cout << boost::format("Kept %d, discarded %d states in diag_denmat")
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
            cout << boost::format(alleig[j] > 1E-3 ? ("%.3f") : ("%.3E")) 
                                % alleig[j];
            cout << ((j != stop) ? ", " : "\n");
	    }
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
        Index act = t.index(1).deprimed();
	if(docomplex && act == Index::IndReIm())
	    act = t.index(2).deprimed();
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
    IQIndex newmid("qlink",iq,In);
    U = IQTensor(active,newmid);
    if(docomplex)
	U *= IQTensor::Complex_1();
    foreach(const ITensor& t, terms) 
	U += t;
    D.ReDimension(m);
    for(int i = 1; i <= m; ++i) 
        D(i) = GET(alleig,alleig.size()-i);
    lastd = D;
    return svdtruncerr;
    } //Real SVDWorker::diag_denmat

