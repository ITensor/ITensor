//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#include "hambuilder.h"
#include "sweeps.h"

namespace itensor {

using std::istream;
using std::ostream;
using std::cout;
using std::endl;
using std::vector;
using std::find;
using std::pair;
using std::make_pair;
using std::string;

template <class Tensor>
MPOt<Tensor>::
MPOt() 
    : 
    Parent(),
    logrefNorm_(DefaultLogRefScale)
    { 
    }
template MPOt<ITensor>::MPOt();
template MPOt<IQTensor>::MPOt();

template <class Tensor>
MPOt<Tensor>::
MPOt(const SiteSet& sites,
     Real _logrefNorm) 
    : 
    Parent(sites)
    { 
    // Norm of psi^2 = 1 = norm = sum of denmat evals. 
    // This translates to Tr{Adag A} = norm.  
    // Ref. norm is Tr{1} = d^N, d = 2 S=1/2, d = 4 for Hubbard, etc
    if(_logrefNorm == DefaultLogRefScale) logrefNorm_ = sites.N();

    //Set all tensors to identity ops
    for(int j = 1; j <= N(); ++j)
        {
        Anc(j) = sites.op("Id",j);
        }
    putMPOLinks(*this);
    }
template
MPOt<ITensor>::
MPOt(const SiteSet& sites, Real _logrefNorm);
template
MPOt<IQTensor>::
MPOt(const SiteSet& sites, Real _logrefNorm);

/*
template<class Tensor> 
void MPOt<Tensor>::
position(int i, const OptSet& opts)
    {
    if(isNull()) Error("position: MPS is null");

    while(l_orth_lim_ < i-1)
        {
        if(l_orth_lim_ < 0) l_orth_lim_ = 0;
        Tensor WF = A(l_orth_lim_+1) * A(l_orth_lim_+2);
        svdBond(l_orth_lim_+1,WF,Fromleft,opts);
        }
    while(r_orth_lim_ > i+1)
        {
        if(r_orth_lim_ > N_+1) r_orth_lim_ = N_+1;
        Tensor WF = A(r_orth_lim_-2) * A(r_orth_lim_-1);
        svdBond(r_orth_lim_-2,WF,Fromright,opts);
        }

    is_ortho_ = true;
    }
template void MPOt<ITensor>::
position(int b, const OptSet& opts);
template void MPOt<IQTensor>::
position(int b, const OptSet& opts);
*/

/*
template <class Tensor>
void MPOt<Tensor>::
orthogonalize(const OptSet& opts)
    {
    //Do a half-sweep to the right, orthogonalizing each bond
    //but do not truncate since the basis to the right might not
    //be ortho (i.e. use the current m).
    //svd_.useOrigM(true);
    int orig_maxm = maxm();
    Real orig_cutoff = cutoff();
    Foreach(Spectrum& spec, spectrum_)
        {
        spec.maxm(MAX_M);
        spec.cutoff(MIN_CUT);
        }

    position(1);
    position(N_);

    //Now basis is ortho, ok to truncate
    Foreach(Spectrum& spec, spectrum_)
        {
        spec.useOrigM(false);
        spec.maxm(orig_maxm);
        spec.cutoff(orig_cutoff);
        }
    position(1);

    is_ortho_ = true;
    }
template
void MPOt<ITensor>::orthogonalize(const OptSet& opts);
template
void MPOt<IQTensor>::orthogonalize(const OptSet& opts);
*/


template <class Tensor>
MPOt<Tensor>& MPOt<Tensor>::
plusEq(const MPOt<Tensor>& other_,
       const OptSet& opts)
    {
    if(doWrite())
        Error("operator+= not supported if doWrite(true)");

    //cout << "calling new orthog in sum" << endl;
    if(!this->isOrtho())
        {
        try { 
            orthogonalize(); 
            }
        catch(const ResultIsZero& rz) 
            { 
            *this = other_;
            return *this;
            }
        }

    if(!other_.isOrtho())
        {
        MPOt<Tensor> other(other_);
        try { 
            other.orthogonalize(); 
            }
        catch(const ResultIsZero& rz) 
            { 
            return *this;
            }
        return addAssumeOrth(other,opts);
        }

    return addAssumeOrth(other_,opts);
    }
template
MPOt<ITensor>& MPOt<ITensor>::plusEq(const MPOt<ITensor>& other, const OptSet&);
template
MPOt<IQTensor>& MPOt<IQTensor>::plusEq(const MPOt<IQTensor>& other, const OptSet&);

int 
findCenter(const IQMPO& psi)
    {
    for(int j = 1; j <= psi.N(); ++j) 
        {
        const IQTensor& A = psi.A(j);
        if(A.r() == 0) Error("Zero rank tensor in IQMPO");
        bool allOut = true;
        Foreach(const IQIndex& I, A.indices())
            {
            //Only look at Link IQIndices
            if(I.type() != Link) continue;

            if(I.dir() != Out)
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
checkQNs(const IQMPO& H)
    {
    const int N = H.N();

    const QN Zero;

    int center = findCenter(H);
    if(center == -1)
        {
        Error("Did not find an ortho. center");
        }

    //Check that all IQTensors have zero div
    //including the ortho. center
    for(int i = 1; i <= N; ++i) 
        {
        if(H.A(i).isNull())
            {
            println("A(",i,") null, QNs not well defined");
            Error("QNs not well defined");
            }
        if(div(H.A(i)) != Zero)
            {
            cout << "At i = " << i << endl;
            Print(H.A(i));
            Error("Non-zero div IQTensor in IQMPO");
            }
        }

    //Check arrows from left edge
    for(int i = 1; i < center; ++i)
        {
        if(rightLinkInd(H,i).dir() != In) 
            {
            println("checkQNs: At site ",i," to the left of the OC, Right side Link not pointing In");
            Error("Incorrect Arrow in IQMPO");
            }
        if(i > 1)
            {
            if(leftLinkInd(H,i).dir() != Out) 
                {
                println("checkQNs: At site ",i," to the left of the OC, Left side Link not pointing Out");
                Error("Incorrect Arrow in IQMPO");
                }
            }
        }

    //Check arrows from right edge
    for(int i = N; i > center; --i)
        {
        if(i < N)
        if(rightLinkInd(H,i).dir() != Out) 
            {
            println("checkQNs: At site ",i," to the right of the OC, Right side Link not pointing Out");
            Error("Incorrect Arrow in IQMPO");
            }
        if(leftLinkInd(H,i).dir() != In) 
            {
            println("checkQNs: At site ",i," to the right of the OC, Left side Link not pointing In");
            Error("Incorrect Arrow in IQMPO");
            }
        }
    }

template <class MPOType>
void 
nmultMPO(const MPOType& Aorig, const MPOType& Borig, MPOType& res,
         const OptSet& opts)
    {
    typedef typename MPOType::TensorT Tensor;
    typedef typename MPOType::IndexT IndexT;
    if(Aorig.N() != Borig.N()) Error("nmultMPO(MPOType): Mismatched N");
    const int N = Borig.N();

    MPOType A(Aorig);
    A.position(1);

    MPOType B;
    if(&Borig == &Aorig)
        {
        B = A;
        }
    else
        {
        B = Borig;
        B.position(1);
        }

    B.primeall();

    res=A;
    res.primelinks(0,4);
    res.mapprime(1,2,Site);

    Tensor clust,nfork;
    vector<int> midsize(N);
    for(int i = 1; i < N; ++i)
        {
        if(i == 1) 
            { 
            clust = A.A(i) * B.A(i); 
            }
        else       
            { 
            clust = nfork * A.A(i) * B.A(i); 
            }

        if(i == N-1) break;

        IndexT oldmid = rightLinkInd(res,i);
        nfork = Tensor(rightLinkInd(A,i),rightLinkInd(B,i),oldmid);

        /*
        if(clust.norm() == 0) // this product gives 0 !!
            { 
            cout << "WARNING: clust.norm()==0 in nmultMPO i=" << i << endl;
            res *= 0;
            return; 
            }
            */

        denmatDecomp(clust, res.Anc(i), nfork,Fromleft,opts);

        IndexT mid = commonIndex(res.A(i),nfork,Link);
        mid.conj();
        midsize[i] = mid.m();
        res.Anc(i+1) = Tensor(mid,conj(res.si(i+1)),prime(res.si(i+1),2),rightLinkInd(res,i+1));
        }

    nfork = clust * A.A(N) * B.A(N);

    /*
    if(nfork.norm() == 0) // this product gives 0 !!
        { 
        cerr << "WARNING: nfork.norm()==0 in nmultMPO\n"; 
        res *= 0;
        return; 
        }
        */

    res.svdBond(N-1,nfork,Fromright);
    res.noprimelink();
    res.mapprime(2,1,Site);
    res.orthogonalize();

    }//void nmultMPO(const MPOType& Aorig, const IQMPO& Borig, IQMPO& res,Real cut, int maxm)
template
void nmultMPO(const MPO& Aorig, const MPO& Borig, MPO& res, const OptSet&);
template
void nmultMPO(const IQMPO& Aorig, const IQMPO& Borig, IQMPO& res,const OptSet& );


template <class Tensor>
void 
zipUpApplyMPO(const MPSt<Tensor>& psi, 
              const MPOt<Tensor>& K, 
              MPSt<Tensor>& res,
              const OptSet& opts)
    {
    typedef typename Tensor::IndexT
    IndexT;

    const
    bool allow_arb_position = opts.getBool("AllowArbPosition",false);

    if(&psi == &res)
        Error("psi and res must be different MPS instances");

    //Real cutoff = opts.getReal("Cutoff",psi.cutoff());
    //int maxm = opts.getInt("Maxm",psi.maxm());

    const int N = psi.N();
    if(K.N() != N) 
        Error("Mismatched N in zipUpApplyMPO");

    if(!psi.isOrtho() || psi.orthoCenter() != 1)
        Error("Ortho center of psi must be site 1");

    if(!allow_arb_position && (!K.isOrtho() || K.orthoCenter() != 1))
        Error("Ortho center of K must be site 1");

#ifdef DEBUG
    checkQNs(psi);
    checkQNs(K);
    /*
    cout << "Checking divergence in zip" << endl;
    for(int i = 1; i <= N; i++)
	div(psi.A(i));
    for(int i = 1; i <= N; i++)
	div(K.A(i));
    cout << "Done Checking divergence in zip" << endl;
    */
#endif

    res = psi; 
    res.primelinks(0,4);
    res.mapprime(0,1,Site);

    Tensor clust,nfork;
    vector<int> midsize(N);
    int maxdim = 1;
    for(int i = 1; i < N; i++)
        {
        if(i == 1) { clust = psi.A(i) * K.A(i); }
        else { clust = nfork * (psi.A(i) * K.A(i)); }
        if(i == N-1) break; //No need to SVD for i == N-1

        IndexT oldmid = rightLinkInd(res,i); assert(oldmid.dir() == Out);
        nfork = Tensor(rightLinkInd(psi,i),rightLinkInd(K,i),oldmid);
        //if(clust.iten_size() == 0)	// this product gives 0 !!
	    //throw ResultIsZero("clust.iten size == 0");
        denmatDecomp(clust, res.Anc(i), nfork,Fromleft,opts);
        IndexT mid = commonIndex(res.A(i),nfork,Link);
        //assert(mid.dir() == In);
        mid.conj();
        midsize[i] = mid.m();
        maxdim = max(midsize[i],maxdim);
        assert(rightLinkInd(res,i+1).dir() == Out);
        res.Anc(i+1) = Tensor(mid,prime(res.si(i+1)),rightLinkInd(res,i+1));
        }
    nfork = clust * psi.A(N) * K.A(N);
    //if(nfork.iten_size() == 0)	// this product gives 0 !!
	//throw ResultIsZero("nfork.iten size == 0");

    res.svdBond(N-1,nfork,Fromright,opts);
    res.noprimelink();
    res.mapprime(1,0,Site);
    res.position(1);
    } //void zipUpApplyMPO
template
void 
zipUpApplyMPO(const MPS& x, const MPO& K, MPS& res, const OptSet& opts);
template
void 
zipUpApplyMPO(const IQMPS& x, const IQMPO& K, IQMPS& res, const OptSet& opts);

//Expensive: scales as m^3 k^3!
template<class Tensor>
void 
exactApplyMPO(const MPSt<Tensor>& x, 
              const MPOt<Tensor>& K, 
              MPSt<Tensor>& res,
              const OptSet& opts)
    {
    typedef typename Tensor::IndexT
    IndexT;
    typedef typename Tensor::CombinerT
    CombinerT;

    int N = x.N();
    if(K.N() != N) Error("Mismatched N in exactApplyMPO");

    if(&res != &x)
        res = x;

    res.Anc(1) = x.A(1) * K.A(1);
    for(int j = 1; j < N; ++j)
        {
        //cout << "exact_applyMPO: step " << j << endl;
        //Compute product of MPS tensor and MPO tensor
        res.Anc(j+1) = x.A(j+1) * K.A(j+1); //m^2 k^2 d^2

        //Add common IQIndices to IQCombiner
        CombinerT comb; 
        Foreach(const IndexT& I, res.A(j).indices())
            {
            if(hasindex(res.A(j+1),I))
                comb.addleft(I);
            }
        comb.init(nameint("a",j));

        //Apply combiner to product tensors
        res.Anc(j) = res.A(j) * comb; //m^3 k^3 d
        res.Anc(j+1) = conj(comb) * res.A(j+1); //m^3 k^3 d
        }
    res.mapprime(1,0,Site);
    res.orthogonalize(opts);
    } //void exact_applyMPO
template
void 
exactApplyMPO(const MPS& x, const MPO& K, MPS& res, const OptSet&);
template
void 
exactApplyMPO(const IQMPS& x, const IQMPO& K, IQMPS& res, const OptSet&);


template<class Tensor>
void
fitApplyMPO(const MPSt<Tensor>& psi,
            const MPOt<Tensor>& K,
            MPSt<Tensor>& res,
            const OptSet& opts)
    {
    fitApplyMPO(1.,psi,K,res,opts);
    }
template
void fitApplyMPO(const MPSt<ITensor>& psi, const MPOt<ITensor>& K, MPSt<ITensor>& res, const OptSet& opts);
template
void fitApplyMPO(const MPSt<IQTensor>& psi, const MPOt<IQTensor>& K, MPSt<IQTensor>& res, const OptSet& opts);

template<class Tensor>
void
fitApplyMPO(Real fac,
            const MPSt<Tensor>& psi,
            const MPOt<Tensor>& K,
            MPSt<Tensor>& res,
            const OptSet& opts)
    {
    const int N = psi.N();
    const int nsweep = opts.getInt("Nsweep",1);
    const bool verbose = opts.getBool("Verbose",false);
    const bool normalize = opts.getBool("Normalize",true);

    const MPSt<Tensor> origPsi(psi);

    vector<Tensor> BK(N+2);

    BK.at(N) = origPsi.A(N)*K.A(N)*conj(prime(res.A(N)));
    for(int n = N-1; n > 2; --n)
        {
        BK.at(n) = BK.at(n+1)*origPsi.A(n)*K.A(n)*conj(prime(res.A(n)));
        }

    res.position(1);

    for(int sw = 1; sw <= nsweep; ++sw)
        {
        for(int b = 1, ha = 1; ha <= 2; sweepnext(b,ha,N))
            {
            if(verbose)
                {
                println("Sweep=",sw,", HS=",ha,", Bond=(",b,",",b+1,")");
                }

            Tensor lwfK = (BK.at(b-1).isNull() ? origPsi.A(b) : BK.at(b-1)*origPsi.A(b));
            lwfK *= K.A(b);
            Tensor rwfK = (BK.at(b+2).isNull() ? origPsi.A(b+1) : BK.at(b+2)*origPsi.A(b+1));
            rwfK *= K.A(b+1);

            Tensor wfK = lwfK*rwfK;
            wfK.noprime();
            wfK *= fac;

            if(normalize) wfK /= wfK.norm();
            res.svdBond(b,wfK,(ha==1?Fromleft:Fromright),opts&Opt("UseSVD",true));

            if(verbose)
                {
                printfln("    Trunc. err=%.1E, States kept=%s",
                         res.spectrum(b).truncerr(),
                         showm(linkInd(res,b)) );
                }

            if(ha == 1)
                BK.at(b) = lwfK * conj(prime(res.A(b)));
            else
                BK.at(b+1) = rwfK * conj(prime(res.A(b+1)));
            }
        }
    }
template
void
fitApplyMPO(Real fac,const MPSt<ITensor>& psi,const MPOt<ITensor>& K,MPSt<ITensor>& res,const OptSet& opts);
template
void
fitApplyMPO(Real fac,const MPSt<IQTensor>& psi,const MPOt<IQTensor>& K,MPSt<IQTensor>& res,const OptSet& opts);

template<class Tensor>
Real
fitApplyMPO(const MPSt<Tensor>& psiA, 
            Real mpofac,
            const MPSt<Tensor>& psiB,
            const MPOt<Tensor>& K,
            MPSt<Tensor>& res,
            const OptSet& opts)
    {
    return fitApplyMPO(1.,psiA,mpofac,psiB,K,res,opts);
    }
template
Real
fitApplyMPO(const MPSt<ITensor>& psiA, Real mpofac,const MPSt<ITensor>& psiB,const MPOt<ITensor>& K,MPSt<ITensor>& res,const OptSet& opts);
template
Real
fitApplyMPO(const MPSt<IQTensor>& psiA, Real mpofac,const MPSt<IQTensor>& psiB,const MPOt<IQTensor>& K,MPSt<IQTensor>& res,const OptSet& opts);

template<class Tensor>
Real
fitApplyMPO(Real mpsfac,
            const MPSt<Tensor>& psiA, 
            Real mpofac,
            const MPSt<Tensor>& psiB,
            const MPOt<Tensor>& K,
            MPSt<Tensor>& res,
            const OptSet& opts)
    {
    if(&psiA == &res || &psiB == &res)
        {
        Error("fitApplyMPO: Result MPS cannot be same as an input MPS");
        }
    const int N = psiA.N();
    const int nsweep = opts.getInt("Nsweep",1);

    vector<Tensor> B(N+2),
                   BK(N+2);

    B.at(N) = psiA.A(N)*conj(prime(res.A(N),Link));
    BK.at(N) = psiB.A(N)*K.A(N)*conj(prime(res.A(N)));
    for(int n = N-1; n > 2; --n)
        {
        B.at(n) = B.at(n+1)*psiA.A(n)*conj(prime(res.A(n),Link));
        BK.at(n) = BK.at(n+1)*psiB.A(n)*K.A(n)*conj(prime(res.A(n)));
        }

    res.position(1);

    for(int sw = 1; sw <= nsweep; ++sw)
        {
        for(int b = 1, ha = 1; ha <= 2; sweepnext(b,ha,N))
            {
            Tensor lwf = (B.at(b-1).isNull() ? psiA.A(b) : B.at(b-1)*psiA.A(b));
            Tensor rwf = (B.at(b+2).isNull() ? psiA.A(b+1) : psiA.A(b+1)*B.at(b+2));

            Tensor lwfK = (BK.at(b-1).isNull() ? psiB.A(b) : BK.at(b-1)*psiB.A(b));
            lwfK *= K.A(b);
            Tensor rwfK = (BK.at(b+2).isNull() ? psiB.A(b+1) : BK.at(b+2)*psiB.A(b+1));
            rwfK *= K.A(b+1);

            Tensor wf = mpsfac*noprime(lwf*rwf) + mpofac*noprime(lwfK*rwfK);
            wf.noprime();

            res.svdBond(b,wf,(ha==1?Fromleft:Fromright),opts&Opt("UseSVD",true));

            if(ha == 1)
                {
                B.at(b) = lwf * conj(prime(res.A(b),Link));
                BK.at(b) = lwfK * conj(prime(res.A(b)));
                }
            else
                {
                B.at(b+1) = rwf * conj(prime(res.A(b+1),Link));
                BK.at(b+1) = rwfK * conj(prime(res.A(b+1)));
                }
            }
        }

    Tensor olp = B.at(3);
    olp *= psiA.A(2);
    olp *= conj(prime(res.A(2),Link));
    olp *= psiA.A(1);
    olp *= conj(prime(res.A(1),Link));

    return olp.toComplex().real();
    }
template
Real
fitApplyMPO(Real mpsfac,const MPSt<ITensor>& psiA, Real mpofac,const MPSt<ITensor>& psiB,const MPOt<ITensor>& K,MPSt<ITensor>& res,const OptSet& opts);
template
Real
fitApplyMPO(Real mpsfac,const MPSt<IQTensor>& psiA, Real mpofac,const MPSt<IQTensor>& psiB,const MPOt<IQTensor>& K,MPSt<IQTensor>& res,const OptSet& opts);

template<class Tensor>
void 
expsmallH(const MPOt<Tensor>& H, 
          MPOt<Tensor>& K, 
          Real tau, 
          Real Etot, 
          Real Kcutoff,
          OptSet opts)
    {
    const int ord = opts.getInt("ExpHOrder",50);
    const bool verbose = opts.getBool("Verbose",false);
    opts.add("Cutoff",MIN_CUT);
    opts.add("Maxm",MAX_M);

    MPOt<Tensor> Hshift(H.sites());
    Hshift.Anc(1) *= -Etot;
    Hshift.plusEq(H,opts);
    Hshift.Anc(1) *= -tau;

    vector<MPOt<Tensor> > xx(2);
    xx.at(0) = MPOt<Tensor>(H.sites());
    xx.at(1) = Hshift;

    //
    // Exponentiate by building up a Taylor series in reverse:
    //      o=1    o=2      o=3      o=4  
    // K = 1-t*H*(1-t*H/2*(1-t*H/3*(1-t*H/4*(...))))
    //
    if(verbose) cout << "Exponentiating H, order: " << endl;
    for(int o = ord; o >= 1; --o)
        {
        if(verbose) 
            {
            cout << o << " "; 
            cout.flush();
            }
        if(o > 1) xx[1].Anc(1) *= 1.0 / o;

        K = sum(xx,opts);
        if(o > 1)
            nmultMPO(K,Hshift,xx[1],opts);
        }
    if(verbose) cout << endl;
    }
template
void 
expsmallH(const MPO& H, MPO& K, Real tau, Real Etot, Real Kcutoff, OptSet opts);
template
void 
expsmallH(const IQMPO& H, IQMPO& K, Real tau, Real Etot, Real Kcutoff, OptSet opts);

template<class Tensor>
void 
expH(const MPOt<Tensor>& H, MPOt<Tensor>& K, 
     Real tau, 
     Real Etot,
     Real Kcutoff, 
     int ndoub,
     OptSet opts)
    {
    const bool verbose = opts.getBool("Verbose",false);
    Real ttau = tau / pow(2.0,ndoub);
    //cout << "ttau in expH is " << ttau << endl;

    Real smallcut = 0.1*Kcutoff*pow(0.25,ndoub);
    expsmallH(H, K, ttau,Etot,smallcut,opts);

    if(verbose) cout << "Starting doubling in expH" << endl;
    for(int doub = 1; doub <= ndoub; ++doub)
        {
        //cout << " Double step " << doub << endl;
        if(doub == ndoub) 
            opts.add("Cutoff",Kcutoff);
        else
            opts.add("Cutoff",0.1 * Kcutoff * pow(0.25,ndoub-doub));
        //cout << "in expH, K.cutoff is " << K.cutoff << endl;
        MPOt<Tensor> KK;
        nmultMPO(K,K,KK,opts);
        K = KK;
        /*
        if(doub == ndoub)
            {
            cout << "step " << doub << ", K is " << endl;
            cout << "K.cutoff, K.maxm are " << K.cutoff SP K.maxm << endl;
            for(int i = 1; i <= N; i++)
                cout << i SP K.A[i];
            }
        */
        }
    }
template
void 
expH(const MPO& H, MPO& K, Real tau, Real Etot,Real Kcutoff, int ndoub, OptSet);
template
void 
expH(const IQMPO& H, IQMPO& K, Real tau, Real Etot,Real Kcutoff, int ndoub, OptSet);

template<class Tensor>
void
applyExpH(const MPSt<Tensor>& psi, 
          const MPOt<Tensor>& H, 
          Real tau, 
          MPSt<Tensor>& res, 
          const OptSet& opts)
    {
    typedef typename Tensor::IndexT
    IndexT;

    typedef MPSt<Tensor>
    MPST;

    if(&psi == &res) Error("Must pass distinct MPS arguments to applyExpH");

    const int order = opts.getInt("Order",10);

    const int N = res.N();
    const int nsweep = opts.getInt("Nsweep",1);

    res.position(1);

    vector<Tensor> lastB(N+2),
                   B(N+2),
                   BH(N+2);

    B.at(N) = psi.A(N)*conj(prime(psi.A(N),Link));
    BH.at(N) = psi.A(N)*H.A(N)*conj(prime(psi.A(N)));
    for(int n = N-1; n > 2; --n)
        {
        B.at(n) = B.at(n+1)*psi.A(n)*conj(prime(psi.A(n),Link));
        BH.at(n) = BH.at(n+1)*psi.A(n)*H.A(n)*conj(prime(psi.A(n)));
        }

    lastB = B;

    MPST last(psi);

    bool up = true;

    for(int ord = order, n = 0; ord >= 1; --ord, ++n)
        {
        const Real mpofac = -tau/(1.*ord);

        if(n > 0) lastB.swap(B);

        for(int sw = 1; sw <= nsweep; ++sw)
            {
            for(int b = 1, ha = 1; ha <= 2; sweepnext(b,ha,N))
                {
                Tensor lwf,rwf,
                       lwfH,rwfH;

                if(up)
                    {
                    lwf = (B.at(b-1).isNull() ? psi.A(b) : B.at(b-1)*psi.A(b));
                    rwf = (B.at(b+2).isNull() ? psi.A(b+1) : B.at(b+2)*psi.A(b+1));

                    lwfH = (BH.at(b-1).isNull() ? last.A(b) : BH.at(b-1)*last.A(b));
                    lwfH *= H.A(b);
                    rwfH = (BH.at(b+2).isNull() ? last.A(b+1) : BH.at(b+2)*last.A(b+1));
                    rwfH *= H.A(b+1);
                    }
                else //dn
                    {
                    lwf = (B.at(b-1).isNull() ? conj(prime(psi.A(b),Link)) : B.at(b-1)*conj(prime(psi.A(b),Link)));
                    rwf = (B.at(b+2).isNull() ? conj(prime(psi.A(b+1),Link)) : B.at(b+2)*conj(prime(psi.A(b+1),Link)));

                    lwfH = (BH.at(b-1).isNull() ? conj(prime(last.A(b))) : BH.at(b-1)*conj(prime(last.A(b))));
                    lwfH *= H.A(b);
                    rwfH = (BH.at(b+2).isNull() ? conj(prime(last.A(b+1))) : BH.at(b+2)*conj(prime(last.A(b+1))));
                    rwfH *= H.A(b+1);
                    }

                Tensor wf = noprime(lwf*rwf) + mpofac*noprime(lwfH*rwfH);
                if(!up) wf.conj();

                res.svdBond(b,wf,(ha==1?Fromleft:Fromright),opts&Opt("UseSVD",true));

                if(up)
                    {
                    if(ha == 1)
                        {
                        B.at(b) = lwf * conj(prime(res.A(b),Link));
                        BH.at(b) = lwfH * conj(prime(res.A(b)));
                        }
                    else
                        {
                        B.at(b+1) = rwf * conj(prime(res.A(b+1),Link));
                        BH.at(b+1) = rwfH * conj(prime(res.A(b+1)));
                        }
                    }
                else //dn
                    {
                    if(ha == 1)
                        {
                        B.at(b) = lwf * res.A(b);
                        BH.at(b) = lwfH * res.A(b);
                        }
                    else
                        {
                        B.at(b+1) = rwf * res.A(b+1);
                        BH.at(b+1) = rwfH * res.A(b+1);
                        }
                    }
                }
            }

        last = res;

        up = !up;

        } // for ord

    }
template
void
applyExpH(const MPSt<ITensor>& psi, const MPOt<ITensor>& H, Real tau, MPSt<ITensor>& res, const OptSet& opts);
template
void
applyExpH(const MPSt<IQTensor>& psi, const MPOt<IQTensor>& H, Real tau, MPSt<IQTensor>& res, const OptSet& opts);

void
putMPOLinks(MPO& W, const OptSet& opts)
    {
    const string pfix = opts.getString("Prefix","l");
    vector<Index> links(W.N());
    for(int b = 1; b < W.N(); ++b)
        {
        links.at(b) = Index(format("%s%d",pfix,b));
        }
    W.Anc(1) *= links.at(1)(1);
    for(int b = 2; b < W.N(); ++b)
        {
        W.Anc(b) *= links.at(b-1)(1);
        W.Anc(b) *= links.at(b)(1);
        }
    W.Anc(W.N()) *= links.at(W.N()-1)(1);
    }

void
putMPOLinks(IQMPO& W, const OptSet& opts)
    {
    QN q;
    const int N = W.N();
    const string pfix = opts.getString("Prefix","l");

    vector<IQIndex> links(N);
    for(int b = 1; b < N; ++b)
        {
        string nm = format("%s%d",pfix,b);
               
        q += div(W.A(b),Opt("Fast"));
        links.at(b) = IQIndex(nm,Index(nm),q);
        }

    W.Anc(1) *= links.at(1)(1);
    for(int b = 2; b < N; ++b)
        {
        W.Anc(b) *= conj(links.at(b-1)(1));
        W.Anc(b) *= links.at(b)(1);
        }
    W.Anc(N) *= conj(links.at(N-1)(1));
    }

}; //namespace itensor
