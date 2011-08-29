#include "mps.h"

void diag_denmat(const ITensor& rho, Real cutoff, int minm, int maxm, ITensor& nU, Vector& D)
{
    assert(rho.r() == 2);
    Index active = rho.index(1); active.noprime();

    //Do the diagonalization
    Index ri = rho.index(1); ri.noprime();
    Matrix R,U; rho.toMatrix11(ri,ri.primed(),R);
    R *= -1.0; EigenValues(R,D,U); D *= -1.0;

    //Truncate
    Real err = 0.0;
    int mp = D.Length();
    while(mp > maxm || (err+D(mp) < cutoff*D(1) && mp > minm)) err += D(mp--);
    svdtruncerr = (D(1) == 0 ? 0.0 : err/D(1));
    D.ReduceDimension(mp); lastd = D;
    Index newmid(active.rawname(),mp,active.type());
    nU = ITensor(active,newmid,U.Columns(1,mp));
}

Vector do_denmat_Real(const ITensor& AA, ITensor& A, ITensor& B, Real cutoff,int minm, int maxm, Direction dir)
{
    Index mid = index_in_common(A,B,Link);
    if(mid.is_null()) mid = Index("mid");

    ITensor& to_orth = (dir==Fromleft ? A : B);
    ITensor& newoc   = (dir==Fromleft ? B : A);

    //Create combiner
    Combiner comb;

    int unique_link = 0;
    foreach(const Index& i, to_orth.indexn())
    if(!(newoc.hasindex(i) || i == IndReIm || i.type() == Virtual))
    { 
        if(i.type() == Link) ++unique_link; 
        comb.addleft(i); 
    }
    foreach(const Index& i, to_orth.index1())
    if(!(newoc.hasindex(i) || i == IndReIm || i.type() == Virtual))
    { 
        if(i.type() == Link) ++unique_link; 
        comb.addleft(i); 
    }

    //Init combiner
    //Index active(mid.rawname());
    //comb.init(active);
    comb.init(mid.rawname());
    Index active = comb.right();

    //Print(AA);
    //Print(to_orth);
    //Print(newoc);
    //Print(unique_link);

    //Check if we're at the edge
    if(unique_link == 0)
    {
        //Handle the right-edge/Fromright and left-edge/Fromleft
        //cases by simply turning the appropriate Combiner into
        //an ITensor and using it as the new edge tensor

        newoc = comb * AA;
        //comb.toITensor(to_orth); to_orth = conj(to_orth);
        to_orth = comb; to_orth.conj();

        Vector eigs_kept(active.m()); eigs_kept = 1.0/active.m();
        return eigs_kept; 
    }

    //Apply combiner....
    ITensor AAc = AA * comb;

    ITensor rho;
    if(AAc.is_complex())
    {
        ITensor re,im;
        AAc.SplitReIm(re,im);
        ITensor rec = conj(re), imc = conj(im);
        rec.primeind(active);
        rho = re * rec;
        imc.primeind(active);
        rho += im * imc;
    }
    else
    {
        ITensor AAcconj = AAc;
        AAcconj.primeind(active);
        rho = AAc * AAcconj;
    }
    assert(rho.r() == 2);

    //Diagonalize the density matrix
    //and form unitary ITensor nU
    ITensor U; Vector D;
    diag_denmat(rho,cutoff,minm,maxm,U,D);

    to_orth = U * comb; //should be conj(comb) with arrows
    newoc   = AAc * conj(U);

    return D;
}

void diag_denmat(const IQTensor& rho, Real cutoff, int minm, int maxm, IQTensor& nU, Vector& eigs_kept)
{
    IQIndex active = rho.finddir(Out);
    assert(active.primelevel == 0);

    vector<Matrix> mmatrix(rho.iten_size());
    vector<Vector> mvector(rho.iten_size());
    vector<Real> alleig;

    Real maxlogfac = -1.0e20;
    const bool do_relative_cutoff = true;
    if(do_relative_cutoff)
    {
        foreach(const ITensor& t, rho.itensors())
        { maxlogfac = max(maxlogfac,t.logfac()); }
    }
    else
	{
        const Real lognormref = 1.3e-14;
        maxlogfac = 2.0*lognormref;
	}

    //1. Diagonalize each ITensor within rho
    int itenind = 0;
    foreach(const ITensor& t, rho.itensors())
	{
        assert(t.index(1).noprime_equals(t.index(2)));
        //if(!t.index(1).noprime_equals(t.index(2)))
        //{ Print(t); Error("Non-symmetric ITensor in density matrix"); }

        t.normlogto(maxlogfac); //Changes the logfac but preserves tensor

        Matrix &U = mmatrix[itenind];
        Vector &d = mvector[itenind];

        //Diag ITensors within rho
        int n = t.index(1).m();
        Matrix M(n,n);
        t.toMatrix11(t.index(2),t.index(1),M);

        M *= -1.0;
        EigenValues(M,d,U);
        d *= -1.0;

        for(int j = 1; j <= n; ++j) 
        { alleig.push_back(d(j)); }
        ++itenind;
	}

    //2. Truncate eigenvalues

    //Sort all eigenvalues from smallest to largest
    //irrespective of quantum numbers
    sort(alleig.begin(),alleig.end());

    //Truncate
    Real e1 = max(alleig.back(),1.0e-60), docut = 0;
    svdtruncerr = 0;
    int m = 0, mkeep = (int)alleig.size();
    if(mkeep > minm)
    for(; m < (int)alleig.size(); m++, mkeep--)
    if(((svdtruncerr += alleig[m]/e1) > cutoff && mkeep <= maxm) || mkeep <= minm)
    { 
        docut = (m > 0 ?  (alleig[m-1] + alleig[m]) * 0.5 : 0);
        svdtruncerr -= alleig[m]/e1;
        break; 
    }
    //cerr << "\nDiscarded " << m << " states in diag_denmat\n";
    m = (int)alleig.size()-m;

    assert(m <= maxm); 
    assert(m < 20000);

    //3. Construct orthogonalized IQTensor nU
    vector<ITensor> terms; terms.reserve(rho.iten_size());
    vector<inqn> iq; iq.reserve(rho.iten_size());
    itenind = 0;
    foreach(const ITensor& t, rho.itensors())
	{
        const Vector& thisD = mvector[itenind];
        int this_m = 1;
        for(; this_m <= thisD.Length(); ++this_m)
        if(thisD(this_m) < docut) { break; }
        --this_m; //since for loop overshoots by 1

        if(this_m == 0) { ++itenind; continue; }

        if(mkeep == 0 && thisD.Length() >= 1) // zero mps, just keep one arb state
        { this_m = 2; mkeep = 1; docut = 1; }

        Index nm("qlink",this_m);
        Index act = t.index(1).deprimed();
        iq.push_back(inqn(nm,active.qn(act)));

        Matrix UU = mmatrix[itenind].Columns(1,this_m);

        assert(act.primelevel == 0);
        assert(active.hasindex(act));
        assert(act.m() == UU.Nrows());

        terms.push_back(ITensor(act,nm,UU));
        ++itenind;
	}
    IQIndex newmid("qlink",iq,In);
    nU = IQTensor(active,newmid);
    foreach(const ITensor& t, terms) nU += t;

    eigs_kept.ReDimension(m);
    for(int i = 1; i <= m; ++i) eigs_kept(i) = alleig[alleig.size()-i];
    lastd = eigs_kept;
} //void diag_denmat

Vector do_denmat_Real(const IQTensor& nA, IQTensor& A, IQTensor& B, Real cutoff, int minm,int maxm, Direction dir)
{
    IQIndex mid = index_in_common(A,B,Link);
    if(mid.is_null()) mid = IQIndex("mid");

    if(nA.iten_size() == 0) Error("zero size in do_denmat_Real(IQTensor)");

    IQTensor& to_orth = (dir==Fromleft ? A : B);
    IQTensor& newoc   = (dir==Fromleft ? B : A);

    int unique_link = 0;
    //Create combiner
    IQCombiner comb;
    foreach(const IQIndex& i, to_orth.iqinds())
    if(!(newoc.hasindex(i) || i == IQIndReIm || i.type() == Virtual))
    { 
        //if(i.type() == Link) { ++unique_link; cerr << "found unique link = " << i << "\n"; }
        if(i.type() == Link) ++unique_link; 
        comb.addleft(i); 
    }

    //Init combiner
    //IQIndex newmid(mid.rawname(),Link,Out);
    //comb.init(newmid);
    comb.doCondense(false);
    comb.init(mid.rawname());
    IQIndex newmid = comb.right();

    //Check if we're at the edge
    //bool edge_case = (to_orth.num_index(Link) <= 1 ? true : false);
    bool edge_case = (unique_link == 0 ? true : false);

    if(edge_case)
    {
        //Handle the right-edge/Fromright and left-edge/Fromleft
        //cases by simply turning the appropriate IQCombiner into
        //an IQTensor and using it as the new edge tensor
        newoc = comb * nA;
        to_orth = comb; to_orth.conj();

        Vector eigs_kept(newmid.m()); eigs_kept = 1.0/newmid.m();
        return eigs_kept; 
    }


    //Apply combiner
    IQTensor nnA;
    nnA = nA * comb;

    //Apply condenser
    IQIndex active("Condensed");
    Condenser cond(newmid,active);
    IQTensor ncA;
    ncA = nnA * cond;

    IQIndex activep(primeBoth,active,4);
    IQTensor rho;
    if(ncA.hasindex(IQIndReIm))
	{
        //cout << "doing complex/Real denmat!" << endl;
        IQTensor cre,cim;
        ncA.SplitReIm(cre,cim);
        IQTensor creconj = conj(cre), cimconj = conj(cim);
        creconj.ind_inc_prime(active,4);
        rho = cre * creconj;
        cimconj.ind_inc_prime(active,4);
        rho += cim * cimconj;
	}
    else
	{
        //conj reverses arrow directions
        IQTensor ncAconj = conj(ncA);
        ncAconj.ind_inc_prime(active,4);
        rho = ncA * ncAconj;
	}

    IQTensor nU; Vector eigs_kept;
    diag_denmat(rho,cutoff,minm,maxm,nU,eigs_kept);

    to_orth = nU * cond * conj(comb);
    newoc  = conj(nU) * ncA;

    return eigs_kept;
} //Vector do_denmat_Real(const IQTensor& nA,...)

Vector do_denmat_Real(const vector<IQTensor>& nA, const IQTensor& A, const IQTensor& B, IQTensor& U,
	Real cutoff, int minm,int maxm, Direction dir, bool donormalize, bool do_relative_cutoff)
{
    // Make a density matrix that is summed over the nA
    IQIndex mid = index_in_common(A,B,Link);
    if(mid.is_null()) mid = IQIndex("mid");

    int num_states = nA.size();
    if(num_states == 0) Error("zero size in do_denmat_Real(vector<IQTensor>)");

    const IQTensor& to_orth = (dir==Fromleft ? A : B);
    const IQTensor& newoc   = (dir==Fromleft ? B : A);

    int unique_link = 0;
    //Create combiner
    IQCombiner comb;
    foreach(const IQIndex& i, to_orth.iqinds())
    if(!(newoc.hasindex(i) || i == IQIndReIm || i.type() == Virtual))
    { 
        if(i.type() == Link) ++unique_link;
        comb.addleft(i); 
    }

    //Check if we're at the edge
    //bool edge_case = (to_orth.num_index(Link) <= 1 ? true : false);
    bool edge_case = (unique_link == 0 ? true : false);

    if(edge_case)
    {
        //Handle the right-edge/Fromright and left-edge/Fromleft
        //cases by simply turning the appropriate IQCombiner into
        //an IQTensor and using it as the new edge tensor

        //IQIndex newmid(mid.rawname(),Link,Out);

        comb.init(mid.rawname());

        U = comb; U.conj(); 

        Vector eigs_kept(comb.right().m()); eigs_kept = 1.0/comb.right().m();
        return eigs_kept; 
    }

    //Combine
    //IQIndex c(mid.rawname());
    //comb.init(c);
    comb.doCondense(false);
    comb.init(mid.rawname());
    IQIndex c = comb.right();
    vector<IQTensor> nnA; nnA.reserve(nA.size());
    foreach(const IQTensor& iqt, nA) nnA.push_back(iqt * comb);

    //Condense
    vector<IQTensor> ncA; ncA.reserve(nnA.size());
    IQIndex active("Condensed");
    Condenser cond(c,active);
    foreach(const IQTensor& iqt, nnA) ncA.push_back(iqt * cond);

    IQIndex activep(primeBoth,active,4);
    IQTensor rho;
    if(ncA.front().hasindex(IQIndReIm))
    {
        //Error("not doing complex 873249827");
        foreach(const IQTensor& iqt, ncA)
        {
            IQTensor iqtre,iqtim;
            iqt.SplitReIm(iqtre,iqtim);
            IQTensor iqtreconj = conj(iqtre), iqtimconj = conj(iqtim);
            iqtreconj.ind_inc_prime(active,4);
            //cout << "k = " << k << endl;
            //cout << "ncA[k] is " << ncA[k] << endl;
            //cout << "ncAconj is " << ncAconj << endl;
            IQTensor r = iqtre * iqtreconj;
            iqtimconj.ind_inc_prime(active,4);
            r += iqtim * iqtimconj;
            //cout << "r is " << r << endl;
            if(num_states == 1) rho = r;
            else                rho += r;
            //cout << "now rho is " << rho << endl;
        }
    }
    else
    {
        foreach(const IQTensor& iqt, ncA)
        {
            IQTensor iqtconj = conj(iqt);
            iqtconj.ind_inc_prime(active,4);
            IQTensor r = iqt * iqtconj;
            if(num_states == 1) rho = r;
            else                rho += r;
        }
    }
    if(rho.iten_size() == 0)
	{
        Error("rho iten_size is 0!!!");
	}

    rho *= 1.0/num_states;

    IQTensor nU; Vector eigs_kept;
    diag_denmat(rho,cutoff,minm,maxm,nU,eigs_kept);

    U = nU * cond * conj(comb);

    return eigs_kept;

} //void do_denmat_Real(const vector<IQTensor>& nA, ... )

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
void getCenterMatrix(ITensor& A, const Index& bond, Real cutoff,int minm, int maxm, ITensor& Lambda, string newbondname)
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

    assert(A.checkDim());
    assert(Lambda.checkDim());
}


void psiHphi(const MPS& psi, const MPO& H, const ITensor& LB, const ITensor& RB, const MPS& phi, Real& re, Real& im) //<psi|H|phi>
{
    int N = psi.NN();
    if(N != phi.NN() || H.NN() < N) Error("mismatched N in psiHphi");
    MPS psiconj(psi);
    for(int i = 1; i <= N; ++i) psiconj.AAnc(i) = conj(primed(psi.AA(i)));
    ITensor L = (LB.is_null() ? phi.AA(1) : LB * phi.AA(1));
    L *= H.AA(1); L *= psiconj.AA(1);
    for(int i = 2; i <= N; ++i)
	{ L *= phi.AA(i); L *= H.AA(i); L *= psiconj.AA(i); }
    if(!RB.is_null()) L *= RB;
    if(L.is_complex())
    {
        if(L.Length() != 1) Error("Non-scalar result in psiHphi.");
        re = L(1) * exp(L.logfac());
        im = L(2) * exp(L.logfac());
    }
    else 
    {
        if(L.Length() != 1) Error("Non-scalar result in psiHphi.");
        re = L.norm();
        im = 0;
    }
}
Real psiHphi(const MPS& psi, const MPO& H, const ITensor& LB, const ITensor& RB, const MPS& phi) //Re[<psi|H|phi>]
{
    Real re,im; psiHphi(psi,H,LB,RB,phi,re,im);
    if(im != 0) cerr << "Real psiHphi: WARNING, dropping non-zero imaginary part of expectation value.\n";
    return re;
}

void plussers(const Index& l1, const Index& l2, Index& sumind, ITensor& first, ITensor& second)
{
    sumind = Index(sumind.rawname(),l1.m()+l2.m(),sumind.type());
    first = ITensor(l1,sumind,1);
    second = ITensor(l2,sumind);
    for(int i = 1; i <= l2.m(); ++i) second(l2(i),sumind(l1.m()+i)) = 1;
}

void plussers(const IQIndex& l1, const IQIndex& l2, IQIndex& sumind, IQTensor& first, IQTensor& second)
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
    first = IQTensor(l1,sumind);
    foreach(const inqn& x, l1.iq())
	{
        Index il1 = x.index;
        Index s1 = l1map[il1];
        ITensor t(il1,s1,1.0);
        first += t;
	}
    second = IQTensor(l2,sumind);
    foreach(const inqn& x, l2.iq())
	{
        Index il2 = x.index;
        Index s2 = l2map[il2];
        ITensor t(il2,s2,1.0);
        second += t;
	}
}


