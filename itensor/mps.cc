#include "mps.h"

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
        AAnc(i) = conj(first[i-1]) * AA(i) * first[i] + conj(second[i-1]) * other.AA(i) * second[i];
    }
    AAnc(N) = conj(first[N-1]) * AA(N) + conj(second[N-1]) * other.AA(N);

    noprimelink();

    position(N);
    position(1);

    return *this;
}
template
MPSt<ITensor>& MPSt<ITensor>::operator+=(const MPSt<ITensor>& other);
template
MPSt<IQTensor>& MPSt<IQTensor>::operator+=(const MPSt<IQTensor>& other);

void diag_denmat(const ITensor& rho, Real cutoff, int minm, int maxm, 
	    ITensor& nU, Vector& D, bool doRelCutoff, LogNumber refNorm)
{
    assert(rho.r() == 2);
    Index active = rho.index(1); active.noprime();

    if(!doRelCutoff) rho.scaleTo(refNorm);

    //Do the diagonalization
    Index ri = rho.index(1); ri.noprime();
    Matrix R,U; rho.toMatrix11(ri,ri.primed(),R);
    R *= -1.0; EigenValues(R,D,U); D *= -1.0;

    //Truncate
    Real err = 0.0;
    int mp = D.Length();
    while(mp > maxm || (err+D(mp) < cutoff*D(1) && mp > minm)) err += D(mp--);
    svdtruncerr = (D(1) == 0 ? 0.0 : err/D(1));
    if(showeigs)
	{
	cout << "doRelCutoff is " << doRelCutoff << endl;
	cout << "refNorm is " << refNorm << endl;
        cout << format("\nKept %d states in diag_denmat\n")% mp;
        cout << format("svdtruncerr = %.2E\n")%svdtruncerr;
        int stop = min(D.Length(),10);
        cout << "Eigs: ";
        for(int j = 1; j <= stop; ++j)
        {
            cout << format(D(j) > 1E-3 ? ("%.3f") : ("%.3E")) % D(j);
            cout << ((j != stop) ? ", " : "\n");
        }
    }

    D.ReduceDimension(mp); lastd = D;
    Index newmid(active.rawname(),mp,active.type());
    nU = ITensor(active,newmid,U.Columns(1,mp));
}

void diag_denmat(const IQTensor& rho, Real cutoff, int minm, int maxm, 
IQTensor& nU, Vector& eigs_kept, 
bool doRelCutoff = false, LogNumber refNorm = 1)
{
    IQIndex active = rho.finddir(Out);
    assert(active.primelevel == 0);

    vector<Matrix> mmatrix(rho.iten_size());
    vector<Vector> mvector(rho.iten_size());
    vector<Real> alleig;

    if(doRelCutoff)
	{
        DO_IF_DEBUG(cout << "Doing relative cutoff\n";)

        Real maxLogNum = -200;
        foreach(const ITensor& t, rho.itensors())
	    maxLogNum = max(maxLogNum,t.scale().logNum());
        assert(maxLogNum > -200);
        assert(maxLogNum <  200);
        refNorm = LogNumber(maxLogNum,1);
        DO_IF_DEBUG(cout << "refNorm = " << refNorm << endl; )
	}
    else
        DO_IF_DEBUG(cout << "Not doing relative cutoff\n";);

    //cerr << format("refNorm = %.1E (lognum = %f, sign = %d)\n\n")
    //%Real(refNorm)%refNorm.logNum()%refNorm.sign();


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
        t.toMatrix11NoScale(t.index(1),t.index(2),M);

        M *= -1;
        EigenValues(M,d,U);
        d *= -1;

        //if(itenind == 1) Print(U.t()*U);
        //Print(d);

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
    //cerr << format("e1 = %.10f\n")%e1;
    svdtruncerr = 0;
    int mdisc = 0, m = (int)alleig.size();
    if(m > minm)
    for(; mdisc < (int)alleig.size(); mdisc++, m--){
    if(((svdtruncerr += GET(alleig,mdisc)/e1) > cutoff && m <= maxm) 
       || m <= minm)
    { 
        if(mdisc > 0)
             { docut = (GET(alleig,mdisc-1) + GET(alleig,mdisc))*0.5; }
        else { docut = 0; }

        //Overshot by one, correct truncerr
        svdtruncerr -= GET(alleig,mdisc)/e1;

        break; 
    }}
    if(showeigs)
    {
        cout << format("\nKept %d, discarded %d states in diag_denmat\n")
                        % m % mdisc;
        cout << format("svdtruncerr = %.2E\n")%svdtruncerr;
        cout << format("docut = %.2E\n")%docut;
        int s = alleig.size();
        int stop = s-min(s,10);
        cout << "Eigs: ";
        for(int j = s-1; j >= stop; --j)
        {
            cout << format(alleig[j] > 1E-3 ? ("%.3f") : ("%.3E")) % alleig[j];
            cout << ((j != stop) ? ", " : "\n");
        }
    }

    assert(m <= maxm); 
    assert(m < 20000);

    //3. Construct orthogonalized IQTensor nU
    vector<ITensor> terms; terms.reserve(rho.iten_size());
    vector<inqn> iq; iq.reserve(rho.iten_size());
    itenind = 0;
    for(IQTensor::const_iten_it it = rho.const_iten_begin(); it != rho.const_iten_end(); ++it)
	{
        const ITensor& t = *it;
        //PrintDat(t);
        const Vector& thisD = GET(mvector,itenind);

        //cerr << "thisD = "; 
        //for(int j = 1; j <= thisD.Length(); ++j)
        //    cerr << format("%.2E ")%thisD(j);
        //cerr << "\n";
        //Print(thisD.Length());

        int this_m = 1;
        for(; this_m <= thisD.Length(); ++this_m)
            if(thisD(this_m) < docut) { break; }
        --this_m; //since for loop overshoots by 1

        if(m == 0 && thisD.Length() >= 1) // zero mps, just keep one arb state
            { this_m = 1; m = 1; docut = 1; }

        if(this_m == 0) { ++itenind; continue; }

        //Print(this_m);

        Index nm("qlink",this_m);
        Index act = t.index(1).deprimed();
        iq.push_back(inqn(nm,active.qn(act)));

        Matrix UU = GET(mmatrix,itenind).Columns(1,this_m);

        //if(itenind == 1)
        //    Print(UU.t()*UU);

        ITensor term(act,nm); term.fromMatrix11(act,nm,UU); 
        //PrintDat(term);
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
Real cutoff, int minm, int maxm, Direction dir, 
bool doRelCutoff, LogNumber refNorm)
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
        if(!(newoc.hasindex(I) || I == Tensor::ReImIndex 
             || I.type() == Virtual))
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
        Vector eigs_kept(comb.right().m()); 
        eigs_kept = 1.0/comb.right().m();
        return eigs_kept; 
    }

    //Apply combiner
    comb.doCondense(true);
    comb.init(mid.rawname());
    Tensor AAc; comb.product(AA,AAc);

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
    diag_denmat(rho,cutoff,minm,maxm,U,eigs_kept,doRelCutoff,refNorm);

    comb.conj();
    comb.product(U,to_orth);
    newoc = conj(U) * AAc;

    return eigs_kept;
}
template Vector 
tensorSVD<ITensor>(const ITensor& AA, ITensor& A, ITensor& B, 
                   Real cutoff, int minm, int maxm, Direction dir, 
                   bool doRelCutoff, LogNumber refScale);
template Vector 
tensorSVD<IQTensor>(const IQTensor& AA, IQTensor& A, IQTensor& B, 
                   Real cutoff, int minm, int maxm, Direction dir, 
                   bool doRelCutoff, LogNumber refScale);

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

