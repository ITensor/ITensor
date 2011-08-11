#include "mps.h"

void diag_denmat(const ITensor& rho, Real cutoff, int minm, int maxm, Matrix& U, Vector& D)
{
    assert(rho.r() == 2);

    //Do the diagonalization
    Index ri = rho.index(1); ri.noprime();
    Matrix R; rho.toMatrix11(ri,ri.primed(),R);
    R *= -1.0; EigenValues(R,D,U); D *= -1.0;

    //Truncate
    Real err = 0.0;
    int mp = D.Length();
    while(mp > maxm || (err+D(mp) < cutoff*D(1) && mp > minm)) err += D(mp--);
    svdtruncerr = (D(1) == 0 ? 0.0 : err/D(1));
    D.ReduceDimension(mp); lastd = D;
    U = U.Columns(1,mp);
}

Vector do_denmat_Real(const ITensor& AA, ITensor& A, ITensor& B, Real cutoff,int minm, int maxm, Direction dir)
{
    Index mid = index_in_common(A,B,Link);
    if(mid.is_null()) mid = Index("mid",1);

    ITensor& to_orth = (dir==Fromleft ? A : B);
    ITensor& newoc   = (dir==Fromleft ? B : A);

    //Create combiner
    Index active("combined");
    Combiner comb(active);

    int unique_link = 0;
    foreach(const Index& i, to_orth.indexn())
    if(!(newoc.hasindex(i) || i == IndReIm || i.type() == Virtual))
    { 
        if(i.type() == Link) ++unique_link; 
        comb.addleft(active,i); 
    }
    foreach(const Index& i, to_orth.index1())
    if(!(newoc.hasindex(i) || i == IndReIm || i.type() == Virtual))
    { 
        if(i.type() == Link) ++unique_link; 
        comb.addleft(active,i); 
    }

    //Check if we're at the edge
    if(unique_link == 0)
    {
        //Handle the right-edge/Fromright and left-edge/Fromleft
        //cases by simply turning the appropriate Combiner into
        //an ITensor and using it as the new edge tensor

        active.setname(mid.rawname());

        newoc = comb * AA;
        comb.toITensor(to_orth); to_orth = conj(to_orth);

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
    Matrix U_; Vector D; diag_denmat(rho,cutoff,minm,maxm,U_,D);

    //Form unitary ITensor nU
    mid = Index(mid.rawname(),D.Length(),mid.type());
    ITensor U(active,mid,U_);

    to_orth = U * comb; //should be conj(comb) with arrows
    newoc   = AAc * conj(U);

    return D;
}

void diag_denmat(IQTensor& rho, Real cutoff, int minm, int maxm, IQTensor& nU, Vector& eigs_kept, bool do_relative_cutoff)
{
    IQIndex active = rho.finddir(Out);
    assert(active.primelevel == 0);

    vector<Matrix> mmatrix(rho.iten_size());
    vector<Vector> mvector(rho.iten_size());
    vector<Real> alleig;

    Real maxlogfac = -1.0e20;
    for(IQTensor::const_iten_it i = rho.const_iten_begin(); i != rho.const_iten_end(); ++i)
    {
        maxlogfac = max(maxlogfac,i->logfac());
    }

    if(!do_relative_cutoff) 
	{
        const Real lognormref = 1.3e-14;
        maxlogfac = 2.0*lognormref;
	}

    //printdat = true; cerr << "rho = " << endl << rho << endl; printdat = false; //DEBUG

    int itenind = 0;
    for(IQTensor::iten_it i = rho.iten_begin(); i != rho.iten_end(); ++i)
	{
        i->normlogto(maxlogfac); //Changes the logfac but preserves tensor

        //Check that all ITensors in rho have indices
        //identical up to their primelevel
        //if(0)
        //if(!i->has_symmetric_indices()) //i->index(1).noprime_equals(i->index(2)))
        //{
        //    cout << "rho is " << rho << endl;
        //    cout << "*i is " << *i << endl << endl;
        //    Error("ITensor failed to have symmetric indices. Perhaps quantum numbers weren't conserved?");
        //}
        //assert(i->has_symmetric_indices());

        //Diag ITensors within rho
        int n = i->index(1).m();
        Matrix M(n,n), U;
        Vector d;
        i->toMatrix11(i->index(2),i->index(1),M);
        if(i->index(2).m() != i->index(1).m())
        {
            cout << "rho is " << rho << endl;
            cout << "*i is " << *i << endl << endl;
            cerr << "i->index(2) = " << i->index(2) << "\n"; 
            cerr << "i->index(1) = " << i->index(1) << "\n"; 
            Error("diag_denmat: density matrix not square.");
        }
    
#ifndef NDEBUG
        //double debug_norm = Norm(Matrix(Matrix(M)-Matrix(M.t())).TreatAsVector());
        //if(debug_norm >= 1E-4)
        //{
            //cerr << "rho = " << endl << rho << endl;   
            //printdat = true;
            //cerr << "i = " << endl << *i << endl;
            //printdat = false;
        //}
#endif

        M *= -1.0;
        EigenValues(M,d,U);
        d *= -1.0;

        for(int j = 1; j <= n; j++) 
        { alleig.push_back(d(j)); }
        mmatrix[itenind] = U;
        mvector[itenind] = d;
        itenind++;
	}

    //Sort all eigenvalues from smallest to largest
    //irrespective of quantum numbers
    sort(alleig.begin(),alleig.end());

    //Truncate eigenvalues
    Real e1 = max(alleig.back(),1.0e-60), su = 0.0, docut = 0.0;
    int m = 0, mkeep = (int)alleig.size();
    if(mkeep > minm)
    for(; m < (int)alleig.size(); m++, mkeep--)
    {
        if(((su += alleig[m]/e1) > cutoff && mkeep <= maxm) || mkeep <= minm)
        //if(((su += alleig[m]/e1) > cutoff && mkeep <= maxm))
        { 
            docut = (m > 0 ?  (alleig[m-1] + alleig[m]) * 0.5 : 0.0); 
            su -= alleig[m]/e1;
            break; 
        }
    }
    m = (int)alleig.size()-m;

    if(m > maxm) 
    {
        cerr << format("minm = %d, maxm = %d\n") % minm % maxm;
        Error("bad m, too big");
    }
    if(m > 20000) Error("bad m, > 20000");
    //if(0 && mkeep == 0) cout << "mkeep = 0" << endl;
    //if(showeigs) { cout << "m, maxm, mkeep, docut, cutoff, e1  are " << m SP maxm SP mkeep SP docut SP cutoff SP e1 << endl; }


    //Construct ITensors for orthogonalized IQTensor (i.e. nU)
    list<ITensor> terms;
    itenind = 0;
    int totkept = 0;
    vector<inqn> iq;
    for(IQTensor::iten_it i = rho.iten_begin(); i != rho.iten_end(); ++i)
	{
        int j = 1;
        for( ; j <= mvector[itenind].Length(); j++)
        { if(mvector[itenind](j) < docut) break; }

        if(mkeep == 0 && mvector[itenind].Length() >= 1)	// zero mps, just keep one arb state
        {
            j = 2;
            mkeep = 1;
            docut = 1.0;
        }
        j -= 1;
        totkept += j;
        if(j == 0) { itenind++; continue; }

        Index nm("qlink",j);
        Index act = i->index(1).deprimed();
        iq.push_back(inqn(nm,active.qn(act)));

        Matrix UU = mmatrix[itenind].Columns(1,j);

        assert(act.primelevel == 0);
        assert(active.hasindex(act));
        assert(act.m() == UU.Nrows());

        ITensor tU(act,nm);
        tU.fromMatrix11(act,nm,UU);
        terms.push_back(tU);
        ++itenind;
	}

    IQIndex newmid("qlink",iq,In);

    nU = IQTensor(active,newmid);
    foreach(const ITensor& t, terms) nU += t;

    eigs_kept.ReDimension(m);
    for(int i = 1; i <= m; ++i) eigs_kept(i) = alleig[alleig.size()-i];
    lastd = eigs_kept;
} //void diag_denmat

Vector do_denmat_Real(const IQTensor& nA, IQTensor& A, IQTensor& B,
    Real cutoff, int minm,int maxm, Direction dir)
{
    const bool do_relative_cutoff = true;
    if(nA.iten_size() == 0) Error("zero size in do_denmat_Real(IQTensor)");

    IQTensor& to_orth = (dir==Fromleft ? A : B);
    IQTensor& newoc   = (dir==Fromleft ? B : A);

    int unique_link = 0;
    //Create combiner
    IQCombiner comb;
    foreach(const IQIndex& i, to_orth.iqindex)
    if(!(newoc.hasindex(i) || i == IQIndReIm || i.type() == Virtual))
    { 
        //if(i.type() == Link) { ++unique_link; cerr << "found unique link = " << i << "\n"; }
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

        IQIndex mid = index_in_common(A,B,Link);
        IQIndex newmid(mid.rawname(),Link,Out);

        comb.init(newmid);

        newoc = comb * nA;
        to_orth = conj(comb.toIQTensor());

        Vector eigs_kept(newmid.m()); eigs_kept = 1.0/newmid.m();
        return eigs_kept; 
    }

    //Init combiner
    IQIndex c("Combined");
    comb.init(c);

    /*
#ifndef NDEBUG
    //do_denmat_Real should always be used to 
    //orthogonalize toward the orthogonality center
    foreach(IQIndex I,comb.left)
    {
        if(I.type() == Link && I.dir() != In)
        {
            cerr << "I = " << I << endl;
            Error("do_denmat_Real direction was not toward orthogonality center.");
        }
    }
#endif
    */

    //Apply combiner
    IQTensor nnA;
    nnA = nA * comb;

    //Apply condenser
    IQIndex active("Condensed");
    Condenser cond(c,active,nnA);
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
    diag_denmat(rho,cutoff,minm,maxm,nU,eigs_kept,do_relative_cutoff);

    to_orth = nU * cond * conj(comb);
    newoc  = conj(nU) * ncA;

    return eigs_kept;
} //Vector do_denmat_Real(const IQTensor& nA,...)

Vector do_denmat_Real(const vector<IQTensor>& nA, const IQTensor& A, const IQTensor& B, IQTensor& U,
	Real cutoff, int minm,int maxm, Direction dir, bool donormalize, bool do_relative_cutoff)
{
    // Make a density matrix that is summed over the nA

    int num_states = nA.size();
    if(num_states == 0) Error("zero size in do_denmat_Real(vector<IQTensor>)");

    const IQTensor& to_orth = (dir==Fromleft ? A : B);
    const IQTensor& newoc   = (dir==Fromleft ? B : A);

    int unique_link = 0;
    //Create combiner
    IQCombiner comb;
    foreach(const IQIndex& i, to_orth.iqindex)
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

        IQIndex mid = index_in_common(A,B,Link);
        IQIndex newmid(mid.rawname(),Link,Out);

        comb.init(newmid);

        U = conj(comb.toIQTensor());

        Vector eigs_kept(newmid.m()); eigs_kept = 1.0/newmid.m();
        return eigs_kept; 
    }

    //Combine
    IQIndex c("Combined");
    comb.init(c);
    list<IQTensor> nnA;
    foreach(const IQTensor& iqt, nA) nnA.push_back(iqt * comb);

    //Condense
    list<IQTensor> ncA;
    IQIndex active("Condensed");
    Condenser cond(c,active,nnA.front());
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
    diag_denmat(rho,cutoff,minm,maxm,nU,eigs_kept,do_relative_cutoff);

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
    Index active("combined");
    Combiner comb(active);

    foreach(const Index& i, A.indexn())
    if(!(i == bond || i == IndReIm || i.type() == Virtual))
    { 
        comb.addleft(active,i); 
    }
    foreach(const Index& i, A.index1())
    if(!(i == bond || i == IndReIm || i.type() == Virtual))
    { 
        comb.addleft(active,i); 
    }

    //Apply combiner....
    ITensor Ac = A * comb;

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
    Matrix U; Vector D; diag_denmat(rho,cutoff,minm,maxm,U,D);

    //Form unitary ITensor nU
    //Index nb(bond.type(),(newbondname == "" ? "c" : newbondname),U.Ncols());
    Index nb("c",D.Length(),bond.type());
    ITensor Uc(active,nb); Uc.fromMatrix11(active,nb,U);

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


