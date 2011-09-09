#include "mps.h"

void diag_denmat(const ITensor& rho, Real cutoff, int minm, int maxm, ITensor& nU, Vector& D, Real refScale = DefaultRefScale)
{
    assert(rho.r() == 2);
    Index active = rho.index(1); active.noprime();

    //if(refScale != DefaultRefScale) rho.scaleTo(refScale);

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

void diag_denmat(const IQTensor& rho, Real cutoff, int minm, int maxm, IQTensor& nU, Vector& eigs_kept, LogNumber refScale = DefaultRefScale)
{
    IQIndex active = rho.finddir(Out);
    assert(active.primelevel == 0);

    vector<Matrix> mmatrix(rho.iten_size());
    vector<Vector> mvector(rho.iten_size());
    vector<Real> alleig;

    if(refScale == DefaultRefScale) 
    {
        Real maxLogNum = -50;
        foreach(const ITensor& t, rho.itensors())
        { maxLogNum = max(maxLogNum,t.scale().logNum()); }
        assert(maxLogNum > -50);
        assert(maxLogNum <  50);
        refScale = LogNumber(maxLogNum,1);
    }

    //cerr << format("refScale = %.1E (lognum = %f, sign = %d)\n")
    //%Real(refScale)%refScale.logNum()%refScale.sign();

//#define USE_REFSCALE

    //1. Diagonalize each ITensor within rho
    int itenind = 0;
    for(IQTensor::const_iten_it it = rho.const_iten_begin(); it != rho.const_iten_end(); ++it)
	{
        const ITensor& t = *it;
        assert(t.index(1).noprime_equals(t.index(2)));
        //if(!t.index(1).noprime_equals(t.index(2)))
        //{ Print(rho); Print(t); Error("Non-symmetric ITensor in density matrix"); }

        t.scaleTo(refNorm);

        Matrix &U = GET(mmatrix,itenind);
        Vector &d = GET(mvector,itenind);

        //Diag ITensors within rho
        int n = t.index(1).m();
        Matrix M(n,n);
#ifdef USE_REFSCALE
        t.toMatrix11NoScale(t.index(2),t.index(1),M);
#else
        t.toMatrix11(t.index(2),t.index(1),M);
#endif

        M *= -1;
        EigenValues(M,d,U);
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
    Real e1 = max(alleig.back(),1.0e-60), docut = 0;
    svdtruncerr = 0;
    int m = 0, mkeep = (int)alleig.size();
    if(mkeep > minm)
    for(; m < (int)alleig.size(); m++, mkeep--){
    if(((svdtruncerr += GET(alleig,m)/e1) > cutoff && mkeep <= maxm) || mkeep <= minm)
    { 
        docut = (m > 0 ?  (GET(alleig,m-1) + GET(alleig,m))*0.5 : 0);
        svdtruncerr -= GET(alleig,m)/e1;
        break; 
    }}
    //cerr << "\nDiscarded " << m << " states in diag_denmat\n";
    m = (int)alleig.size()-m;

    assert(m <= maxm); 
    assert(m < 20000);

    //3. Construct orthogonalized IQTensor nU
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

        if(mkeep == 0 && thisD.Length() >= 1) // zero mps, just keep one arb state
        { this_m = 2; mkeep = 1; docut = 1; }
        --this_m; //since for loop overshoots by 1

        if(this_m == 0) { ++itenind; continue; }

        Index nm("qlink",this_m);
        Index act = t.index(1).deprimed();
        iq.push_back(inqn(nm,active.qn(act)));

        Matrix UU = GET(mmatrix,itenind).Columns(1,this_m);

        ITensor term(act,nm); term.fromMatrix11(act,nm,UU); 
#ifdef USE_REFSCALE
        term *= refScale;
#endif
        terms.push_back(term);

        ++itenind;
	}
    IQIndex newmid("qlink",iq,In);
    nU = IQTensor(active,newmid);
    foreach(const ITensor& t, terms) nU += t;

    eigs_kept.ReDimension(m);
    for(int i = 1; i <= m; ++i) eigs_kept(i) = GET(alleig,alleig.size()-i);
    lastd = eigs_kept;
} //void diag_denmat

template<class Tensor>
Vector tensorSVD(const Tensor& AA, Tensor& A, Tensor& B, 
Real cutoff, int minm, int maxm, Direction dir, LogNumber refScale)
{
    typedef typename Tensor::IndexT IndexT;
    typedef typename Tensor::CombinerT CombinerT;

    if(AA.vec_size() == 0) 
    {
        A *= 0;
        B *= 0;
        Vector eigs(1); eigs = 1;
        return eigs;
        //Error("tensorSVD(Tensor): input tensor had zero size.");
    }

    IndexT mid = index_in_common(A,B,Link);
    if(mid.is_null()) mid = IndexT("mid");

    Tensor& to_orth = (dir==Fromleft ? A : B);
    Tensor& newoc   = (dir==Fromleft ? B : A);

    CombinerT comb;

    int unique_link = 0; //number of Links unique to to_orth
    for(int j = 1; j <= to_orth.r(); ++j) 
    { 
        const IndexT& I = to_orth.index(j);
        //if(I.type() == Link) debug2 = true;
        if(!(newoc.hasindex(I) || I == Tensor::ReImIndex || I.type() == Virtual))
        {
            if(I.type() == Link) ++unique_link;
            //cerr << "Adding left index " << I << "\n";
            comb.addleft(I);
        }
        //debug2 = false;
    }

    //Check if we're at the edge
    if(unique_link == 0)
    {
        comb.init(mid.rawname());
        assert(comb.check_init());
        comb.product(AA,newoc);
        to_orth = comb; to_orth.conj();
        Vector eigs_kept(comb.right().m()); eigs_kept = 1.0/comb.right().m();
        return eigs_kept; 
    }

    //Apply combiner
    comb.doCondense(true);
    comb.init(mid.rawname());
    Tensor AAc; comb.product(AA,AAc);

    /*
    cerr << "WARNING: checking norms -- slow\n";

    if(fabs(AA.norm()-AAc.norm())/AA.norm() > 1E-5)
    {
        cerr << format("AA.norm() = %.10f\n")%AA.norm();
        cerr << format("AAc.norm() = %.10f\n")%AAc.norm();
        cerr << format("rel diff = %.10f\n")%(fabs(AA.norm()-AAc.norm())/AA.norm());
        Error("Incorrect norm for combined tensor.");
    }

    if(fabs(AA.sumels()-AAc.sumels())/fabs(AA.sumels()) > 1E-5)
    {
        Print(AA.sumels());
        Print(AAc.sumels());
        cerr << format("rel diff = %.3E\n")%(fabs(AA.sumels()-AAc.sumels())/fabs(AA.sumels()));
        Error("Incorrect total for condensed tensor.");
    }
    */

    assert(LogNumber(AAc.norm()).isFinite());

    if(LogNumber(AAc.norm()).isNan())
    {
        Print(LogNumber(AAc.norm()));
        Print(AAc.norm());
        PrintDat(AAc);
        Error("Norm was nan");
    }

    const IndexT& active = comb.right();

    Tensor rho;
    if(AAc.is_complex())
    {
        Tensor re,im;
        AAc.SplitReIm(re,im);
        rho = re; rho.conj(); rho.primeind(active);
        rho *= re;
        im *= conj(primeind(im,active));
        rho += im;
    }
    else 
    { 
        Tensor AAcc = conj(AAc); 
        AAcc.primeind(active); 
        rho = AAc*AAcc; 
    }
    assert(rho.r() == 2);

    Tensor U; Vector eigs_kept;
    diag_denmat(rho,cutoff,minm,maxm,U,eigs_kept,refScale);

    comb.conj();
    comb.product(U,to_orth);
    newoc = conj(U) * AAc;

    return eigs_kept;
}
template Vector tensorSVD<ITensor>(const ITensor& AA, ITensor& A, ITensor& B, Real cutoff, int minm, int maxm, Direction dir, LogNumber refScale);
template Vector tensorSVD<IQTensor>(const IQTensor& AA, IQTensor& A, IQTensor& B, Real cutoff, int minm, int maxm, Direction dir, LogNumber refScale);


/*
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
*/

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

