#include "mps.h"

void plussers(const Index& l1, const Index& l2, Index& sumind, 
              ITensor& first, ITensor& second)
{
    sumind = Index(sumind.rawname(),l1.m()+l2.m(),sumind.type());
    first = ITensor(l1,sumind,1);
    second = ITensor(l2,sumind);
    for(int i = 1; i <= l2.m(); ++i) second(l2(i),sumind(l1.m()+i)) = 1;
}

void plussers(const IQIndex& l1, const IQIndex& l2, IQIndex& sumind, 
              IQTensor& first, IQTensor& second)
{
    map<Index,Index> l1map, l2map;
    vector<inqn> iq;
    foreach(const inqn& x, l1.iq())
	{
        Index ii = x.index;
        Index jj(ii.name(),ii.m(),ii.type());
        l1map[ii] = jj;
        iq.push_back(inqn(jj,x.qn));
	}
    foreach(const inqn& x, l2.iq())
	{
        Index ii = x.index;
        Index jj(ii.name(),ii.m(),ii.type());
        l2map[ii] = jj;
        iq.push_back(inqn(jj,x.qn));
	}
    sumind = IQIndex(sumind,iq);
    first = IQTensor(conj(l1),sumind);
    foreach(const inqn& x, l1.iq())
	{
        Index il1 = x.index;
        Index s1 = l1map[il1];
        ITensor t(il1,s1,1.0);
        first += t;
	}
    second = IQTensor(conj(l2),sumind);
    foreach(const inqn& x, l2.iq())
	{
        Index il2 = x.index;
        Index s2 = l2map[il2];
        ITensor t(il2,s2,1.0);
        second += t;
	}
}

template <class Tensor>
MPSt<Tensor>& MPSt<Tensor>::operator+=(const MPSt<Tensor>& other)
{
    primelinks(0,4);

    vector<Tensor> first(N), second(N);
    for(int i = 1; i < N; ++i)
    {
        IndexT l1 = this->RightLinkInd(i);
        IndexT l2 = other.RightLinkInd(i);
        IndexT r(l1.rawname());
        plussers(l1,l2,r,first[i],second[i]);
    }

    AAnc(1) = AA(1) * first[1] + other.AA(1) * second[1];
    for(int i = 2; i < N; ++i)
    {
        AAnc(i) = conj(first[i-1]) * AA(i) * first[i] 
                  + conj(second[i-1]) * other.AA(i) * second[i];
    }
    AAnc(N) = conj(first[N-1]) * AA(N) + conj(second[N-1]) * other.AA(N);

    noprimelink();

    //cerr << "WARNING: skipping orthogonalize in operator+=\n";
    orthogonalize();

    return *this;
}
template
MPSt<ITensor>& MPSt<ITensor>::operator+=(const MPSt<ITensor>& other);
template
MPSt<IQTensor>& MPSt<IQTensor>::operator+=(const MPSt<IQTensor>& other);

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
    while(mp > maxm_ || (svdtruncerr+D(mp) < cutoff_*D(1) && mp > minm_)) 
        { svdtruncerr += D(mp--); }
    svdtruncerr = (D(1) == 0 ? 0 : svdtruncerr/D(1));
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
        //assert(t.index(1).noprime_equals(t.index(2)));
        if(!t.index(1).noprime_equals(t.index(2)))
        { Print(rho); Print(t); Error("Non-symmetric ITensor in density matrix"); }

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
    if(m > minm_)
    for(; mdisc < (int)alleig.size(); mdisc++, m--){
    if(((svdtruncerr += GET(alleig,mdisc)/e1) > cutoff_ && m <= maxm_) 
       || m <= minm_)
    { 
        if(mdisc > 0)
             { docut = (alleig.at(mdisc-1) + alleig[mdisc])*0.5; }
        else { docut = -1; }

        //Overshot by one, correct truncerr
        svdtruncerr -= alleig[mdisc]/e1;

        break; 
    }}
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
            if(thisD(this_m) < docut) { break; }
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
    foreach(const ITensor& t, terms) U += t;

    D.ReDimension(m);
    for(int i = 1; i <= m; ++i) 
        { D(i) = GET(alleig,alleig.size()-i); }
    lastd = D;
    return svdtruncerr;
} //Real SVDWorker::diag_denmat


/* getCenterMatrix:
 * 
 *                    s                   s
 *                    |                   |
 * Decomposes A = -<--A-->- bond into -<--U-<-- -<--Lambda-->- bond
 *
 * A is replaced with the unitary U and Lambda is diagonal.
 * If A is the OC of an MPS, Lambda will contain the Schmidt weights. 
 *
 */
 /*
void getCenterMatrix(ITensor& A, const Index& bond, Real cutoff,int minm, int maxm, ITensor& Lambda, string newbondname = "")
{
    //Create combiner
    Combiner comb;
    foreach(const Index& i, A.indexn())
    if(!(i == bond || i == IndReIm || i.type() == Virtual))
    { 
        comb.addleft(i); 
    }
    foreach(const Index& i, A.index1())
    if(!(i == bond || i == IndReIm || i.type() == Virtual))
    { 
        comb.addleft(i); 
    }
    comb.init("combined");
    Index active = comb.right();

    //Apply combiner....
    //comb.init(active);
    ITensor Ac = comb * A;

    ITensor rho;
    if(Ac.is_complex())
    {
        ITensor re,im;
        Ac.SplitReIm(re,im);
        ITensor rec = conj(re), imc = conj(im);
        rec.primeind(active);
        rho = re * rec;
        imc.primeind(active);
        rho += im * imc;
    }
    else rho = Ac * primeind(Ac,active);
    assert(rho.r() == 2);

    //Diagonalize & truncate the density matrix
    ITensor Uc; Vector D; diag_denmat(rho,cutoff,minm,maxm,Uc,D);

    Lambda = conj(Uc) * Ac;
    A = Uc * comb; //should be conj(comb) with arrows

}
*/

int collapseCols(const Vector& Diag, Matrix& M)
{
    int nr = Diag.Length(), nc = int(Diag.sumels());
    assert(nr != 0);
    if(nc == 0) return nc;
    M = Matrix(nr,nc); M = 0;
    int c = 0;
    for(int r = 1; r <= nr; ++r)
    if(Diag(r) == 1) { M(r,++c) = 1; }
    return nc;
}

void convertToIQ(const BaseModel& model, const vector<ITensor>& A, vector<IQTensor>& qA, QN totalq, Real cut)
{
    const int N = A.size()-1;
    qA.resize(A.size());
    const bool is_mpo = A[1].hasindex(model.si(1).primed());
    const int Dim = model.dim();
    const int PDim = (is_mpo ? Dim : 1);

    vector<IQIndex> linkind(N);

    typedef map<QN,Vector>::value_type qD_vt;
    map<QN,Vector> qD; //Diags of compressor matrices by QN
    typedef map<QN,vector<ITensor> >::value_type qt_vt;
    map<QN,vector<ITensor> > qt; //ITensor blocks by QN
    typedef map<QN,ITensor>::value_type qC_vt;
    map<QN,ITensor> qC; //Compressor ITensors by QN
    ITensor block;
    vector<ITensor> nblock;
    vector<inqn> iq;

    QN q;

    qC[totalq] = ITensor(); //Represents Virtual index
    //First value of prev_q below set to totalq

    const int show_s = 0;

    Index bond, prev_bond;
    for(int s = 1; s <= N; ++s)
    {
        qD.clear(); qt.clear();
        if(s > 1) prev_bond = index_in_common(A[s-1],A[s],Link);
        if(s < N) bond = index_in_common(A[s],A[s+1],Link);

        if(s == show_s) { PrintDat(A[s]); }

        foreach(const qC_vt& x, qC) {
        const QN& prev_q = x.first; const ITensor& comp = x.second; 
        for(int n = 1; n <= Dim;  ++n)
        for(int u = 1; u <= PDim; ++u)
        {
            q = (is_mpo ? prev_q+model.si(s).qn(n)-model.si(s).qn(u) : prev_q-model.si(s).qn(n));

            //For the last site, only keep blocks 
            //compatible with specified totalq i.e. q=0 here
            if(s == N && q != QN()) continue;

            //Set Site indices of A[s] and its previous Link Index
            block = A[s];
            if(s != 1) block *= conj(comp);
            block *= model.si(s)(n);
            if(is_mpo) block *= model.siP(s)(u);

            //Initialize D Vector (D records which values of
            //the right Link Index to keep for the current QN q)
            int count = qD.count(q);
            Vector& D = qD[q];
            if(count == 0) { D.ReDimension(bond.m()); D = 0; }

            if(s == show_s)
            {
                cerr << boost::format("For n = %d\n")%n;
                cerr << boost::format("Got a block with norm %.10f\n")%block.norm();
                cerr << boost::format("bond.m() = %d\n")%bond.m();
                PrintDat(block);
                if(s != 1) PrintDat(comp);
            }

            bool keep_block = false;
            if(s == N) keep_block = true;
            else
            {
                if(bond.m() == 1 && block.norm() != 0) { D = 1; keep_block = true; }
                else
                {
                    ITensor summed_block;
                    if(s==1) summed_block = block;
                    else
                    {
                        //Here we sum over the previous link index
                        //which is already ok, analyze the one to the right
                        assert(comp.r()==2);
                        Index new_ind = (comp.index(1)==prev_bond ? comp.index(2) : comp.index(1));
                        summed_block = ITensor(new_ind,1) * block;
                    }
                    //cerr << boost::format("s = %d, bond=")%s << bond << "\n";
                    //summed_block.print("summed_block");

                    Real rel_cut = -1;
                    for(int j = 1; j <= bond.m(); ++j)
                    { rel_cut = max(fabs(summed_block.val1(j)),rel_cut); }
                    assert(rel_cut >= 0);
                    //Real rel_cut = summed_block.norm()/summed_block.vec_size();
                    rel_cut *= cut;
                    //cerr << "rel_cut == " << rel_cut << "\n";

                    if(rel_cut > 0)
                    for(int j = 1; j <= bond.m(); ++j)
                    if(fabs(summed_block.val1(j)) > rel_cut) 
                    { D(j) = 1; keep_block = true; }
                }
            } //else (s != N)

            //if(!keep_block && q == totalq) { D(1) = 1; keep_block = true; }

            if(keep_block)
            {
                qD[q] = D;

                if(is_mpo) 
                {
                block.addindex1(conj(model.si(s)(n).index()));
                block.addindex1(model.siP(s)(u).index());
                }
                else { block.addindex1(model.si(s)(n).index()); }

                qt[q].push_back(block);

                if(s==show_s)
                {
                block.print("Kept block",ShowData);
                cerr << "D = " << D << "\n";
                }
            }
        }}

        qC.clear();

        foreach(const qt_vt& x, qt)
        {
            const vector<ITensor>& blks = x.second;
            if(blks.size() != 0)
            {
                q = x.first; 
                if(s == N) 
                { foreach(const ITensor& t, blks) nblock.push_back(t); }
                else
                {
                    Matrix M; int mm = collapseCols(qD[q],M);
                    if(s==show_s)
                    {
                        cerr << boost::format("Adding block, mm = %d\n")%mm;
                        q.print("q");
                        cerr << "qD[q] = " << qD[q] << "\n";
                        cerr << "M = \n" << M << "\n";
                        int count = 0;
                        foreach(const ITensor& t, blks) 
                        t.print((boost::format("t%02d")%(++count)).str(),ShowData);
                    }
                    //string qname = (boost::format("ql%d(%+d:%d:%s)")%s%q.sz()%q.Nf()%(q.Nfp() == 0 ? "+" : "-")).str();
                    string qname = (boost::format("ql%d(%+d:%d)")%s%q.sz()%q.Nf()).str();
                    Index qbond(qname,mm);
                    ITensor compressor(bond,qbond,M);
                    foreach(const ITensor& t, blks) nblock.push_back(t * compressor);
                    iq.push_back(inqn(qbond,q));
                    qC[q] = compressor;
                }
            }
        }

        if(s != N) 
        { 
            if(iq.empty()) 
            {
                cerr << "At site " << s << "\n";
                Error("convertToIQ: no compatible QNs to put into Link.");
            }
            linkind[s] = IQIndex(nameint("qL",s),iq); iq.clear(); 
        }
        if(s == 1)
            qA[s] = (is_mpo ? IQTensor(conj(model.si(s)),model.siP(s),linkind[s]) : IQTensor(model.si(s),linkind[s]));
        else if(s == N)
            qA[s] = (is_mpo ? IQTensor(conj(linkind[s-1]),conj(model.si(s)),model.siP(s)) 
                                    : IQTensor(conj(linkind[s-1]),model.si(s)));
        else
            qA[s] = (is_mpo ? IQTensor(conj(linkind[s-1]),conj(model.si(s)),model.siP(s),linkind[s]) 
                                    : IQTensor(conj(linkind[s-1]),model.si(s),linkind[s]));

        foreach(const ITensor& nb, nblock) { qA[s] += nb; } nblock.clear();

        if(s==show_s)
        {
        qA[s].print((boost::format("qA[%d]")%s).str(),ShowData);
        Error("Stopping");
        }

    } //for loop over s

    IQIndex Center("Center",Index("center",1,Virtual),totalq,In);
    qA[1].addindex1(Center);
}

/*
template <class Tensor> 
template <class IQMPSType> 
void MPSt<Tensor>::convertToIQ(IQMPSType& iqpsi, QN totalq, Real cut) const
{
    assert(model_ != 0);
    const ModelT& sst = *model_;

    iqpsi = IQMPSType(sst,maxm,cutoff);

    if(!A[1].hasindex(si(1))) Error("convertToIQ: incorrect primelevel for conversion");
    bool is_mpo = A[1].hasindex(si(1).primed());
    const int Dim = si(1).m();
    const int PDim = (is_mpo ? Dim : 1);

    vector<IQIndex> linkind(N);

    typedef map<QN,Vector>::value_type qD_vt;
    map<QN,Vector> qD; //Diags of compressor matrices by QN
    typedef map<QN,vector<ITensor> >::value_type qt_vt;
    map<QN,vector<ITensor> > qt; //ITensor blocks by QN
    typedef map<QN,ITensor>::value_type qC_vt;
    map<QN,ITensor> qC; //Compressor ITensors by QN
    ITensor block;
    vector<ITensor> nblock;
    vector<inqn> iq;

    QN q;

    qC[totalq] = ITensor(); //Represents Virtual index
    //First value of prev_q below set to totalq

    const int show_s = 0;

    Index bond, prev_bond;
    for(int s = 1; s <= N; ++s)
    {
        if(debug3) cerr << "s = " << s << "\n";

        qD.clear(); qt.clear();
        if(s > 1) prev_bond = LinkInd(s-1); 
        if(s < N) bond = LinkInd(s);

        if(s == show_s) 
        {
            PrintDat(A[s]);
        }

        foreach(const qC_vt& x, qC) {
        const QN& prev_q = x.first; const ITensor& comp = x.second; 
        for(int n = 1; n <= Dim;  ++n)
        for(int u = 1; u <= PDim; ++u)
        {
            q = (is_mpo ? prev_q+si(s).qn(n)-si(s).qn(u) : prev_q-si(s).qn(n));

            //For the last site, only keep blocks 
            //compatible with specified totalq i.e. q=0 here
            if(s == N && q != QN()) continue;

            //Set Site indices of A[s] and its previous Link Index
            block = A[s];
            if(s != 1) block *= conj(comp);
            block *= si(s)(n);
            if(is_mpo) block *= siP(s)(u);

            //Initialize D Vector (D records which values of
            //the right Link Index to keep for the current QN q)
            int count = qD.count(q);
            Vector& D = qD[q];
            if(count == 0) { D.ReDimension(bond.m()); D = 0; }

            if(s == show_s)
            {
                cerr << boost::format("For n = %d\n")%n;
                cerr << boost::format("Got a block with norm %.10f\n")%block.norm();
                cerr << boost::format("bond.m() = %d\n")%bond.m();
                PrintDat(block);
                if(s != 1) PrintDat(comp);
            }

            bool keep_block = false;
            if(s == N) keep_block = true;
            else
            {
                if(bond.m() == 1 && block.norm() != 0) { D = 1; keep_block = true; }
                else
                {
                    ITensor summed_block;
                    if(s==1) summed_block = block;
                    else
                    {
                        //Here we sum over the previous link index
                        //which is already ok, analyze the one to the right
                        assert(comp.r()==2);
                        Index new_ind = (comp.index(1)==prev_bond ? comp.index(2) : comp.index(1));
                        summed_block = ITensor(new_ind,1) * block;
                    }
                    //cerr << boost::format("s = %d, bond=")%s << bond << "\n";
                    //summed_block.print("summed_block");

                    Real rel_cut = -1;
                    for(int j = 1; j <= bond.m(); ++j)
                    { rel_cut = max(fabs(summed_block.val1(j)),rel_cut); }
                    assert(rel_cut >= 0);
                    //Real rel_cut = summed_block.norm()/summed_block.vec_size();
                    rel_cut *= cut;
                    //cerr << "rel_cut == " << rel_cut << "\n";

                    if(rel_cut > 0)
                    for(int j = 1; j <= bond.m(); ++j)
                    if(fabs(summed_block.val1(j)) > rel_cut) 
                    { D(j) = 1; keep_block = true; }
                }
            } //else (s != N)

            //if(!keep_block && q == totalq) { D(1) = 1; keep_block = true; }

            if(keep_block)
            {
                qD[q] = D;

                if(is_mpo) 
                {
                block.addindex1(conj(si(s)(n).index()));
                block.addindex1(siP(s)(u).index());
                }
                else { block.addindex1(si(s)(n).index()); }

                qt[q].push_back(block);

                if(s==show_s)
                {
                block.print("Kept block",ShowData);
                cerr << "D = " << D << "\n";
                }
            }
        }}

        qC.clear();

        foreach(const qt_vt& x, qt)
        {
            const vector<ITensor>& blks = x.second;
            if(blks.size() != 0)
            {
                q = x.first; 
                if(s == N) 
                { foreach(const ITensor& t, blks) nblock.push_back(t); }
                else
                {
                    Matrix M; int mm = collapseCols(qD[q],M);
                    if(s==show_s)
                    {
                        cerr << boost::format("Adding block, mm = %d\n")%mm;
                        q.print("q");
                        cerr << "qD[q] = " << qD[q] << "\n";
                        cerr << "M = \n" << M << "\n";
                        int count = 0;
                        foreach(const ITensor& t, blks) 
                        t.print((boost::format("t%02d")%(++count)).str(),ShowData);
                    }
                    //string qname = (boost::format("ql%d(%+d:%d:%s)")%s%q.sz()%q.Nf()%(q.Nfp() == 0 ? "+" : "-")).str();
                    string qname = (boost::format("ql%d(%+d:%d)")%s%q.sz()%q.Nf()).str();
                    Index qbond(qname,mm);
                    ITensor compressor(bond,qbond,M);
                    foreach(const ITensor& t, blks) nblock.push_back(t * compressor);
                    iq.push_back(inqn(qbond,q));
                    qC[q] = compressor;
                }
            }
        }

        if(s != N) 
        { 
            if(iq.empty()) 
            {
                cerr << "At site " << s << "\n";
                Error("convertToIQ: no compatible QNs to put into Link.");
            }
            linkind[s] = IQIndex(nameint("qL",s),iq); iq.clear(); 
        }
        if(s == 1)
        {
            iqpsi.AAnc(s) = (is_mpo ? IQTensor(conj(si(s)),siP(s),linkind[s]) : IQTensor(si(s),linkind[s]));
            IQIndex Center("Center",Index("center",1,Virtual),totalq,In);
            iqpsi.AAnc(1).addindex1(Center);
        }
        else if(s == N)
        {
            iqpsi.AAnc(s) = (is_mpo ? IQTensor(conj(linkind[s-1]),conj(si(s)),siP(s)) 
                                    : IQTensor(conj(linkind[s-1]),si(s)));
        }
        else
        {
            iqpsi.AAnc(s) = (is_mpo ? IQTensor(conj(linkind[s-1]),conj(si(s)),siP(s),linkind[s]) 
                                    : IQTensor(conj(linkind[s-1]),si(s),linkind[s]));
        }

        foreach(const ITensor& nb, nblock) { iqpsi.AAnc(s) += nb; } nblock.clear();

        if(0) //try to get this working ideally
        if(!is_mpo && s > 1) 
        {
            IQTensor AA = iqpsi.bondTensor(s-1);
            iqpsi.doSVD(s-1,AA,Fromleft);
        }

        if(s==show_s)
        {
        iqpsi.AA(s).print((boost::format("qA[%d]")%s).str(),ShowData);
        Error("Stopping");
        }

    } //for loop over s

    assert(check_QNs(iqpsi));

} //void convertToIQ(IQMPSType& iqpsi) const
*/

