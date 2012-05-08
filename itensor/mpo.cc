//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#include "mpo.h"
using namespace std;
using boost::format;

template<class Tensor> 
void MPOt<Tensor>::
position(int i, bool preserve_shape)
    {
    if(isNull()) Error("position: MPS is null");
    while(l_orth_lim_ < i-1)
        {
        if(l_orth_lim_ < 0) l_orth_lim_ = 0;
        Tensor WF = AA(l_orth_lim_+1) * AA(l_orth_lim_+2);
        svdBond(l_orth_lim_+1,WF,Fromleft,preserve_shape);
        }
    while(r_orth_lim_ > i+1)
        {
        if(r_orth_lim_ > N+1) r_orth_lim_ = N+1;
        Tensor WF = AA(r_orth_lim_-2) * AA(r_orth_lim_-1);
        svdBond(r_orth_lim_-2,WF,Fromright,preserve_shape);
        }
    }
template void MPOt<ITensor>::
position(int b, bool preserve_shape);
template void MPOt<IQTensor>::
position(int b, bool preserve_shape);

template <class Tensor>
void MPOt<Tensor>::
svdBond(int b, const Tensor& AA, Direction dir, bool preserve_shape)
    {
    if(preserve_shape)
        {
        //The idea of the preserve_shape flag is to 
        //leave any external indices of the MPO on the
        //tensors they originally belong to
        Error("preserve_shape not currently implemented");
        }

    if(dir == Fromleft && b-1 > l_orth_lim_)
        {
        std::cout << boost::format("b=%d, l_orth_lim_=%d")
                %b%l_orth_lim_ << std::endl;
        Error("b-1 > l_orth_lim_");
        }
    if(dir == Fromright && b+2 < r_orth_lim_)
        {
        std::cout << boost::format("b=%d, r_orth_lim_=%d")
                %b%r_orth_lim_ << std::endl;
        Error("b+2 < r_orth_lim_");
        }

    SparseT D;
    svd_.svd(b,AA,A[b],D,A[b+1]);

    //Push singular values/amplitudes
    //to the right or left as requested
    //and update orth_lims
    if(dir == Fromleft)
        {
        A[b+1] *= D;

        l_orth_lim_ = b;
        if(r_orth_lim_ < b+2) r_orth_lim_ = b+2;
        }
    else //dir == Fromright
        {
        A[b] *= D;

        if(l_orth_lim_ > b-1) l_orth_lim_ = b-1;
        r_orth_lim_ = b+1;
        }
    }
template void MPOt<ITensor>::
svdBond(int b, const ITensor& AA, Direction dir, bool preserve_shape);
template void MPOt<IQTensor>::
svdBond(int b, const IQTensor& AA, Direction dir, bool preserve_shape);

int 
findCenter(const IQMPO& psi)
    {
    for(int j = 1; j <= psi.NN(); ++j) 
        {
        const IQTensor& A = psi.AA(j);
        if(A.r() == 0) Error("Zero rank tensor in IQMPO");
        bool allOut = true;
        for(int i = 1; i <= A.r(); ++i)
            {
            //Only look at Link IQIndices
            if(A.index(i).type() != Link) continue;

            if(A.index(i).dir() != Out)
                {
                allOut = false;
                break;
                }
            }

        //Found the ortho. center
        if(allOut) return j;
        }
    return -1;
    }

void
checkQNs(const IQMPO& psi)
    {
    const int N = psi.NN();

    QN zero;

    int center = findCenter(psi);
    //std::cerr << boost::format("Found the OC at %d\n") % center;
    if(center == -1)
        {
        Error("Did not find an ortho. center");
        }

    //Check that all IQTensors have zero div
    //including the ortho. center
    for(int i = 1; i <= N; ++i) 
        {
        if(psi.AA(i).isNull())
            {
            std::cerr << boost::format("AA(%d) null, QNs not well defined\n")%i;
            Error("QNs not well defined");
            }
        try {
            checkQNs(psi.AA(i));
            }
        catch(const ITError& e)
            {
            std::cerr << "At i = " << i << "\n";
            throw e;
            }
        if(psi.AA(i).div() != zero)
            {
            std::cerr << "At i = " << i << "\n";
            Print(psi.AA(i));
            Error("Non-zero div IQTensor in IQMPO");
            }
        }

    //Check arrows from left edge
    for(int i = 1; i < center; ++i)
        {
        if(psi.RightLinkInd(i).dir() != In) 
            {
            std::cerr << boost::format("checkQNs: At site %d to the left of the OC, Right side Link not pointing In\n")%i;
            Error("Incorrect Arrow in IQMPO");
            }
        if(i > 1)
            {
            if(psi.LeftLinkInd(i).dir() != Out) 
                {
                std::cerr << boost::format("checkQNs: At site %d to the left of the OC, Left side Link not pointing Out\n")%i;
                Error("Incorrect Arrow in IQMPO");
                }
            }
        }

    //Check arrows from right edge
    for(int i = N; i > center; --i)
        {
        if(i < N)
        if(psi.RightLinkInd(i).dir() != Out) 
            {
            std::cerr << boost::format("checkQNs: At site %d to the right of the OC, Right side Link not pointing Out\n")%i;
            Error("Incorrect Arrow in IQMPO");
            }
        if(psi.LeftLinkInd(i).dir() != In) 
            {
            std::cerr << boost::format("checkQNs: At site %d to the right of the OC, Left side Link not pointing In\n")%i;
            Error("Incorrect Arrow in IQMPO");
            }
        }
    }

template <class MPOType>
void 
nmultMPO(const MPOType& Aorig, const MPOType& Borig, MPOType& res,Real cut, int maxm)
    {
    typedef typename MPOType::TensorT Tensor;
    typedef typename MPOType::IndexT IndexT;
    if(Aorig.NN() != Borig.NN()) Error("nmultMPO(MPOType): Mismatched N");
    int N = Borig.NN();
    MPOType A(Aorig), B(Borig);

    SVDWorker svd = A.svd();
    svd.cutoff(cut);
    svd.maxm(maxm);

    A.position(1);
    B.position(1);
    B.primeall();

    res=A;
    res.primelinks(0,4);
    res.mapprime(1,2,primeSite);

    Tensor clust,nfork;
    vector<int> midsize(N);
    for(int i = 1; i < N; ++i)
        {
        if(i == 1) 
            { 
            clust = A.AA(i) * B.AA(i); 
            }
        else       
            { 
            clust = nfork * A.AA(i) * B.AA(i); 
            }

        if(i == N-1) break;

        IndexT oldmid = res.RightLinkInd(i);
        nfork = Tensor(A.RightLinkInd(i),B.RightLinkInd(i),oldmid);

        /*
        if(clust.norm() == 0) // this product gives 0 !!
            { 
            cerr << boost::format("WARNING: clust.norm()==0 in nmultMPO (i=%d).\n")%i; 
            res *= 0;
            return; 
            }
            */

        svd.denmatDecomp(i,clust, res.AAnc(i), nfork,Fromleft);

        IndexT mid = index_in_common(res.AA(i),nfork,Link);
        mid.conj();
        midsize[i] = mid.m();
        res.AAnc(i+1) = Tensor(mid,conj(res.si(i+1)),primed(res.si(i+1),2),res.RightLinkInd(i+1));
        }

    nfork = clust * A.AA(N) * B.AA(N);

    /*
    if(nfork.norm() == 0) // this product gives 0 !!
        { 
        cerr << "WARNING: nfork.norm()==0 in nmultMPO\n"; 
        res *= 0;
        return; 
        }
        */

    res.doSVD(N-1,nfork,Fromright);
    res.noprimelink();
    res.mapprime(2,1,primeSite);
    res.cutoff(cut);
    res.orthogonalize();

    }//void nmultMPO(const MPOType& Aorig, const IQMPO& Borig, IQMPO& res,Real cut, int maxm)
template
void nmultMPO(const MPO& Aorig, const MPO& Borig, MPO& res,Real cut, int maxm);
template
void nmultMPO(const IQMPO& Aorig, const IQMPO& Borig, IQMPO& res,Real cut, int maxm);

void 
napplyMPO(const IQMPS& x, const IQMPO& K, IQMPS& res, Real cutoff, int maxm, bool allow_arb_position)
    {
    if(cutoff < 0) cutoff = x.cutoff();
    if(maxm < 0) maxm = x.maxm();
    int N = x.NN();
    if(K.NN() != N) Error("Mismatched N in napplyMPO");
    if(x.rightLim() > 3)
        {
        Error("bad rightLim for x");
        }
    if(!allow_arb_position && K.rightLim() > 3)
        {
        //cerr << "K is " << endl << K << endl;
        Error("bad rightLim for K");
        }

    SVDWorker svd = K.svd();
    svd.cutoff(cutoff);
    svd.maxm(maxm);

    res = x; 
    res.maxm(maxm); 
    res.cutoff(cutoff);
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
	    throw ResultIsZero("clust.iten size == 0");
        svd.denmatDecomp(i,clust, res.AAnc(i), nfork,Fromleft);
        IQIndex mid = index_in_common(res.AA(i),nfork,Link);
        assert(mid.dir() == In);
        mid.conj();
        midsize[i] = mid.m();
        maxdim = max(midsize[i],maxdim);
        assert(res.RightLinkInd(i+1).dir() == Out);
        res.AAnc(i+1) = IQTensor(mid,primed(res.si(i+1)),res.RightLinkInd(i+1));
        }
    nfork = clust * x.AA(N) * K.AA(N);
    if(nfork.iten_size() == 0)	// this product gives 0 !!
	throw ResultIsZero("nfork.iten size == 0");

    res.doSVD(N-1,nfork,Fromright);
    res.noprimelink();
    res.mapprime(1,0,primeSite);
    res.position(1);
    res.maxm(x.maxm()); 
    res.cutoff(x.cutoff());

    } //void napplyMPO

//Expensive: scales as m^3 k^3!
void 
exact_applyMPO(const IQMPS& x, const IQMPO& K, IQMPS& res)
    {
    int N = x.NN();
    if(K.NN() != N) Error("Mismatched N in exact_applyMPO");

    res = x;

    res.AAnc(1) = x.AA(1) * K.AA(1);
    for(int j = 1; j < N; ++j)
        {
        //cerr << boost::format("exact_applyMPO: step %d\n") % j;
        //Compute product of MPS tensor and MPO tensor
        res.AAnc(j+1) = x.AA(j+1) * K.AA(j+1); //m^2 k^2 d^2

        //Add common IQIndices to IQCombiner
        IQCombiner comb; comb.doCondense(false);
        Foreach(const IQIndex& I, res.AA(j).iqinds())
        if(res.AA(j+1).hasindex(I) && I != IQIndex::IndReIm())
            { assert(I.dir() == Out); comb.addleft(I);}
        comb.init(nameint("a",j));

        //Apply combiner to product tensors
        res.AAnc(j) = res.AA(j) * comb; //m^3 k^3 d
        res.AAnc(j+1) = conj(comb) * res.AA(j+1); //m^3 k^3 d
        }
    res.mapprime(1,0,primeSite);
    //res.orthogonalize();
    } //void exact_applyMPO

