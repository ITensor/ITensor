#include "mpo.h"

void nmultMPO(const IQMPO& Aorig, const IQMPO& Borig, IQMPO& res,Real cut, int maxm)
{
    if(Aorig.NN() != Borig.NN()) Error("nmultMPO(IQMPO): Mismatched N");
    int N = Borig.NN();
    IQMPO A(Aorig), B(Borig);

    A.position(1);
    B.position(1);
    B.primeall();

    res=A;
    res.primelinks(0,4);
    res.mapprime(1,2,primeSite);

    IQTensor clust,nfork;
    vector<int> midsize(N);
    for(int i = 1; i < N; ++i)
	{
        if(i == 1) { clust = A.AA(i) * B.AA(i); }
        else       { clust = nfork * A.AA(i) * B.AA(i); }
        if(i == N-1) break;

        IQIndex oldmid = res.RightLinkInd(i);
        nfork = IQTensor(A.RightLinkInd(i),B.RightLinkInd(i),oldmid);
        if(clust.iten_size() == 0)	// this product gives 0 !!
        { 
            cerr << boost::format("WARNING: clust.iten_size()==0 in nmultMPO (i=%d).\n")%i; 
            res *= 0;
            return; 
        }
        tensorSVD(clust, res.AAnc(i), nfork,cut,1,maxm,Fromleft,A.doRelCutoff(),A.refNorm());
        IQIndex mid = index_in_common(res.AA(i),nfork,Link);
        assert(mid.dir() == In);
        mid.conj();
        midsize[i] = mid.m();
        assert(res.RightLinkInd(i+1).dir() == Out);
        assert(res.si(i+1).dir() == Out);
        res.AAnc(i+1) = IQTensor(mid,conj(res.si(i+1)),res.si(i+1).primed().primed(),res.RightLinkInd(i+1));
	}

    nfork = clust * A.AA(N) * B.AA(N);
    if(nfork.iten_size() == 0)	// this product gives 0 !!
    { 
        cerr << "WARNING: nfork.iten_size()==0 in nmultMPO\n"; 
        res *= 0;
        return; 
    }

    res.doSVD(N-1,nfork,Fromright,false);
    res.noprimelink();
    res.mapprime(2,1,primeSite);
    res.cutoff = cut;
    res.position(N);
    res.position(1);

}//void nmultMPO(const IQMPO& Aorig, const IQMPO& Borig, IQMPO& res,Real cut, int maxm)

void napplyMPO(const IQMPS& x, const IQMPO& K, IQMPS& res, Real cutoff, int maxm)
{
    if(cutoff < 0) cutoff = x.cutoff;
    if(maxm < 0) maxm = x.maxm;
    int N = x.NN();
    if(K.NN() != N) Error("Mismatched N in napplyMPO");
    if(x.right_lim() > 3)
    {
        cerr << "x is " << endl << x << endl;
        Error("bad right_lim for x");
    }
    if(K.right_lim() > 3)
    {
        //cerr << "K is " << endl << K << endl;
        Error("bad right_lim for K");
    }

    res = x; res.maxm = maxm; res.cutoff = cutoff;
    res.primelinks(0,4);
    res.mapprime(0,1,primeSite);

    IQTensor clust,nfork;
    vector<int> midsize(N);
    int maxdim = 1;
    for(int i = 1; i < N; i++)
	{
        if(i == 1) { clust = x.AA(i) * K.AA(i); }
        else { clust = nfork * (x.AA(i) * K.AA(i)); }
        if(i == N-1) break; //No need to SVD for i == N-1

        IQIndex oldmid = res.RightLinkInd(i); assert(oldmid.dir() == Out);
        nfork = IQTensor(x.RightLinkInd(i),K.RightLinkInd(i),oldmid);
        if(clust.iten_size() == 0)	// this product gives 0 !!
        { res *= 0; return; }
        tensorSVD(clust, res.AAnc(i), nfork,
                  cutoff,1,maxm,Fromleft,K.doRelCutoff(),K.refNorm());
        IQIndex mid = index_in_common(res.AA(i),nfork,Link);
        assert(mid.dir() == In);
        mid.conj();
        midsize[i] = mid.m();
        maxdim = max(midsize[i],maxdim);
        assert(res.RightLinkInd(i+1).dir() == Out);
        res.AAnc(i+1) = IQTensor(mid,res.si(i+1).primed(),res.RightLinkInd(i+1));
	}
    nfork = clust * x.AA(N) * K.AA(N);
    if(nfork.iten_size() == 0)	// this product gives 0 !!
	{ res *= 0; return; }

    res.doSVD(N-1,nfork,Fromright,false);
    res.noprimelink();
    res.mapprime(1,0,primeSite);
    res.position(1);
    res.maxm = x.maxm; res.cutoff = x.cutoff;

} //void napplyMPO

//Expensive: scales as m^3 k^3!
void exact_applyMPO(const IQMPS& x, const IQMPO& K, IQMPS& res)
{
    int N = x.NN();
    if(K.NN() != N) Error("Mismatched N in exact_applyMPO");

    res = x;
    res.position(1);

    res.AAnc(1) = x.AA(1) * K.AA(1);
    for(int j = 1; j < N; ++j)
	{
        //cerr << boost::format("exact_applyMPO: step %d\n") % j;
        //Compute product of MPS tensor and MPO tensor
        res.AAnc(j+1) = x.AA(j+1) * K.AA(j+1); //m^2 k^2 d^2

        //Add common IQIndices to IQCombiner
        IQCombiner comb; comb.doCondense(false);
        foreach(const IQIndex& I, res.AA(j).iqinds())
        if(res.AA(j+1).hasindex(I) && I != IQIndReIm && I.type() != Virtual)
        { assert(I.dir() == Out); comb.addleft(I);}
        comb.init(nameint("a",j));

        //Apply combiner to product tensors
        res.AAnc(j) = res.AA(j) * comb; //m^3 k^3 d
        res.AAnc(j+1) = conj(comb) * res.AA(j+1); //m^3 k^3 d
	}
    res.mapprime(1,0,primeSite);
    //res.position(1);
} //void exact_applyMPO
