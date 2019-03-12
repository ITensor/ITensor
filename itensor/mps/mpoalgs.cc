//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/util/print_macro.h"
#include "itensor/mps/mpo.h"
#include "itensor/mps/localop.h"

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

template<class MPOType>
void 
nmultMPO(MPOType const& Aorig, 
         MPOType const& Borig, 
         MPOType& res,
         Args args)
    {
    using Tensor = typename MPOType::TensorT;

    if(!args.defined("Cutoff")) args.add("Cutoff",1E-14);

    if(Aorig.N() != Borig.N()) Error("nmultMPO(MPOType): Mismatched N");
    const int N = Borig.N();

    auto A = Aorig;
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
    auto siA = uniqueIndex(A.A(1),B.A(1),A.A(2));
    auto siB = uniqueIndex(B.A(1),A.A(1),B.A(2));
    res.Aref(1) = Tensor(siA,siB,linkInd(A,1));

    //println("\n\n");
    //Print(A);
    //Print(B);

    Tensor clust,nfork;
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

        nfork = Tensor(linkInd(A,i),linkInd(B,i),linkInd(res,i));

        //printfln("Before denmatDecomp i=%d:",i);
        //Print(clust);
        //Print(res.A(i));
        //Print(nfork);

        denmatDecomp(clust,res.Aref(i),nfork,Fromleft,args);

        //printfln("After denmatDecomp i=%d:",i);
        //Print(res.A(i));
        //Print(nfork);

        //println("<>--<>--<>--<>");
        //PAUSE;

        auto mid = commonIndex(res.A(i),nfork,Link);
        mid.dag();
        auto siA = uniqueIndex(A.A(i+1),A.A(i),A.A(i+2),B.A(i+1));
        auto siB = uniqueIndex(B.A(i+1),B.A(i),B.A(i+2),A.A(i+1));
        //Print(siA);
        //Print(siB);
        res.Aref(i+1) = Tensor(mid,siA,siB,rightLinkInd(res,i+1));
        }

    nfork = clust * A.A(N) * B.A(N);

    res.svdBond(N-1,nfork,Fromright, args);
    for(auto i : range1(N))
        {
        if(i < N)
            {
            auto l = linkInd(res,i);
            res.Aref(i).mapprime(l,l.primeLevel(),0);
            res.Aref(i+1).mapprime(l,l.primeLevel(),0);
            }
        res.Aref(i).mapprime(2,1);
        }
    res.orthogonalize();
    }
template
void nmultMPO(const MPO& Aorig, const MPO& Borig, MPO& res, Args);
template
void nmultMPO(const IQMPO& Aorig, const IQMPO& Borig, IQMPO& res,Args);

template<class Tensor>
MPSt<Tensor>
applyMPO(MPOt<Tensor> const& K,
         MPSt<Tensor> const& x,
         Args const& args)
    {
    auto method = args.getString("Method","DensityMatrix");

    //This is done here because fitApplyMPO() has a different default behavior
    //(for backwards compatability)
    auto normalize = args.getBool("Normalize",false);
    auto argsp = args;
    argsp.add("Normalize=",normalize);

    MPSt<Tensor> res;
    if(method == "DensityMatrix")
        res = exactApplyMPO(K,x,argsp);
    else if(method == "Fit")
        res = fitApplyMPO(x,K,argsp);
    else
        Error("applyMPO currently supports the following methods: 'DensityMatrix' (previously called with exactApplyMPO), 'Fit' (previously called with fitApplyMPO)");

    return res;
    }
template
MPS
applyMPO(MPO const& K, MPS const& x, Args const&);
template
IQMPS
applyMPO(IQMPO const& K, IQMPS const& x, Args const&);

template<class Tensor>
MPSt<Tensor>
applyMPO(MPOt<Tensor> const& K,
         MPSt<Tensor> const& x,
         MPSt<Tensor> const& x0,
         Args const& args)
    {
    auto method = args.getString("Method","Fit");

    //This is done here because fitApplyMPO() has a different default behavior
    //(for backwards compatability)
    auto normalize = args.getBool("Normalize",false);
    auto argsp = args;
    argsp.add("Normalize=",normalize);

    MPSt<Tensor> res = x0;
    if(method == "DensityMatrix")
        Error("applyMPO method 'DensityMatrix' does not accept an input MPS");
    else if(method == "Fit")
        fitApplyMPO(x,K,res,argsp);
    else
        Error("applyMPO currently supports the following methods: 'DensityMatrix' (previously called with exactApplyMPO), 'Fit' (previously called with fitApplyMPO)");

    return res;
    }
template
MPS
applyMPO(MPO const& K, MPS const& x, MPS const& x0, Args const&);
template
IQMPS
applyMPO(IQMPO const& K, IQMPS const& x, IQMPS const& x0, Args const&);

template<class Tensor>
MPSt<Tensor>
exactApplyMPO(MPOt<Tensor> const& K,
              MPSt<Tensor> const& psi,
              Args const& args)
    {
    auto cutoff = args.getReal("Cutoff",1E-13);
    auto ignore_degeneracy = args.getBool("IgnoreDegeneracy",true);
    auto maxm_set = args.defined("Maxm");
    auto verbose = args.getBool("Verbose",false);
    auto normalize = args.getBool("Normalize",false);
    auto siteType = getIndexType(args,"SiteType",Site);
    auto linkType = getIndexType(args,"LinkType",Link);

    auto dargs = Args{"Cutoff",cutoff,"IgnoreDegeneracy",ignore_degeneracy,"Verbose",verbose};
    if(maxm_set) dargs.add("Maxm",args.getInt("Maxm"));

    if(noprime(findtype(K.A(1),Site)) != findtype(psi.A(1),Site))
        {
        Error("MPS and MPO have different site indices in exactApplyMPO");
        }

    int plev = 14741;

    auto res = psi;

    auto N = psi.N();

    //Set up conjugate psi and K
    auto psic = psi;
    auto Kc = K;
    for(auto j : range1(N)) 
        {
        //Modify prime levels of psic and Kc
        if(j == 1)
            {
            auto ci = commonIndex(psi.A(1),psi.A(2),linkType);
            psic.Aref(j) = dag(mapprime(psi.A(j),siteType,0,2,ci,0,plev));
            ci = commonIndex(Kc.A(1),Kc.A(2),linkType);
            Kc.Aref(j) = dag(mapprime(K.A(j),siteType,0,2,ci,0,plev));
            }
        else
            {
            psic.Aref(j) = dag(mapprime(psi.A(j),siteType,0,2,linkType,0,plev));
            Kc.Aref(j) = dag(mapprime(K.A(j),siteType,0,2,linkType,0,plev));
            }
        }

    //Build environment tensors from the left
    if(verbose) print("Building environment tensors...");
    auto E = std::vector<Tensor>(N+1);
    E.at(1) = psi.A(1)*K.A(1)*Kc.A(1)*psic.A(1);
    for(int j = 2; j < N; ++j)
        {
        E.at(j) = E.at(j-1)*psi.A(j)*K.A(j)*Kc.A(j)*psic.A(j);
        //assert(rank(E[j])==4);
        }
    if(verbose) println("done");

    //O is the representation of the product of K*psi in the new MPS basis
    auto O = psi.A(N)*K.A(N);
    O.noprime(siteType);

    auto rho = E.at(N-1) * O * dag(prime(O,plev));
    Tensor U,D;
    dargs.add("IndexName=",nameint("a",N));
    auto spec = diagHermitian(rho,U,D,dargs);
    if(verbose) printfln("  j=%02d truncerr=%.2E m=%d",N-1,spec.truncerr(),commonIndex(U,D).m());

    res.Aref(N) = dag(U);

    O = O*U*psi.A(N-1)*K.A(N-1);
    O.noprime(siteType);

    for(int j = N-1; j > 1; --j)
        {
        if(not maxm_set)
            {
            //Infer maxm from bond dim of original MPS
            //times bond dim of MPO
            //i.e. upper bound on rank of rho
            auto cip = commonIndex(psi.A(j),E.at(j-1));
            auto ciw = commonIndex(K.A(j),E.at(j-1));
            auto maxm = (cip) ? cip.m() : 1l;
            maxm *= (ciw) ? ciw.m() : 1l;
            dargs.add("Maxm",maxm);
            }
        rho = E.at(j-1) * O * dag(prime(O,plev));
        dargs.add("IndexName=",nameint("a",j));
        auto spec = diagHermitian(rho,U,D,dargs);
        O = O*U*psi.A(j-1)*K.A(j-1);
        O.noprime(siteType);
        res.Aref(j) = dag(U);
        if(verbose) printfln("  j=%02d truncerr=%.2E m=%d",j,spec.truncerr(),commonIndex(U,D).m());
        }

    if(normalize) O /= norm(O);
    res.Aref(1) = O;
    res.leftLim(0);
    res.rightLim(2);

    return res;
    }
template
MPS
exactApplyMPO(MPO const& K, MPS const& x, Args const&);
template
IQMPS
exactApplyMPO(IQMPO const& K, IQMPS const& x, Args const&);

template<class Tensor>
MPSt<Tensor>
exactApplyMPO(MPSt<Tensor> const& x,
              MPOt<Tensor> const& K,
              Args const& args)
    {
    return exactApplyMPO(K,x,args);
    }
template
MPS
exactApplyMPO(MPS const& x, MPO const& K, Args const&);
template
IQMPS
exactApplyMPO(IQMPS const& x, IQMPO const& K, Args const&);

template<class Tensor>
void 
exactApplyMPO(MPSt<Tensor> const& x, 
              MPOt<Tensor> const& K, 
              MPSt<Tensor>      & res,
              Args const& args)
    {
    res = exactApplyMPO(K,x);
    }
template
void 
exactApplyMPO(MPS const& x, MPO const& K, MPS& res, Args const&);
template
void 
exactApplyMPO(IQMPS const& x, IQMPO const& K, IQMPS& res, Args const&);

template<class Tensor>
MPSt<Tensor>
fitApplyMPO(MPSt<Tensor> const& psi,
            MPOt<Tensor> const& K,
            Args const& args)
    {
    MPSt<Tensor> res;
    fitApplyMPO(1.,psi,K,res,args);
    return res;
    }
template
MPSt<ITensor> fitApplyMPO(MPSt<ITensor> const& psi, MPOt<ITensor> const& K, Args const& args);
template
MPSt<IQTensor> fitApplyMPO(MPSt<IQTensor> const& psi, MPOt<IQTensor> const& K, Args const& args);

template<class Tensor>
void
fitApplyMPO(MPSt<Tensor> const& psi,
            MPOt<Tensor> const& K,
            MPSt<Tensor>& res,
            Args const& args)
    {
    fitApplyMPO(1.,psi,K,res,args);
    }
template
void fitApplyMPO(const MPSt<ITensor>& psi, const MPOt<ITensor>& K, MPSt<ITensor>& res, const Args& args);
template
void fitApplyMPO(const MPSt<IQTensor>& psi, const MPOt<IQTensor>& K, MPSt<IQTensor>& res, const Args& args);

template<class Tensor>
void
fitApplyMPO(Real fac,
            MPSt<Tensor> const& psi,
            MPOt<Tensor> const& K,
            MPSt<Tensor>& res,
            Args const& args)
    {
    auto nsweep = args.getInt("Nsweep",1);
    Sweeps sweeps(nsweep);
    auto cutoff = args.getReal("Cutoff",-1);
    if(cutoff >= 0) sweeps.cutoff() = cutoff;
    auto maxm = args.getInt("Maxm",-1);
    if(maxm >= 1) sweeps.maxm() = maxm;
    fitApplyMPO(fac,psi,K,res,sweeps,args);
    }
template
void
fitApplyMPO(Real fac,
            MPSt<ITensor> const& psi,
            MPOt<ITensor> const& K,
            MPSt<ITensor>& res,
            Args const& args);
template
void
fitApplyMPO(Real fac,
            MPSt<IQTensor> const& psi,
            MPOt<IQTensor> const& K,
            MPSt<IQTensor>& res,
            Args const& args);

template<class Tensor>
void
fitApplyMPO(Real fac,
            MPSt<Tensor> const& psi,
            MPOt<Tensor> const& K,
            MPSt<Tensor>& res,
            Sweeps const& sweeps,
            Args args)
    {
    auto N = psi.N();
    auto verbose = args.getBool("Verbose",false);
    auto normalize = args.getBool("Normalize",true);

    const auto origPsi = psi;

    if(not res) res = origPsi;
    res.position(1);

    auto BK = vector<Tensor>(N+2);
    BK.at(N) = origPsi.A(N)*K.A(N)*dag(prime(res.A(N)));
    for(auto n = N-1; n > 2; --n)
        {
        BK.at(n) = BK.at(n+1)*origPsi.A(n)*K.A(n)*dag(prime(res.A(n)));
        }

    for(auto sw : range1(sweeps.nsweep()))
        {
        args.add("Sweep",sw);
        args.add("Cutoff",sweeps.cutoff(sw));
        args.add("Minm",sweeps.minm(sw));
        args.add("Maxm",sweeps.maxm(sw));
        args.add("Noise",sweeps.noise(sw));

        for(int b = 1, ha = 1; ha <= 2; sweepnext(b,ha,N))
            {
            if(verbose)
                {
                printfln("Sweep=%d, HS=%d, Bond=(%d,%d)",sw,ha,b,b+1);
                }

            auto lwfK = (BK.at(b-1) ? BK.at(b-1)*origPsi.A(b) : origPsi.A(b));
            lwfK *= K.A(b);
            auto rwfK = (BK.at(b+2) ? BK.at(b+2)*origPsi.A(b+1) : origPsi.A(b+1));
            rwfK *= K.A(b+1);

            auto wfK = lwfK*rwfK;
            wfK.noprime();
            wfK *= fac;

            if(normalize) wfK /= norm(wfK);
            auto PH = LocalOp<Tensor>(K.A(b),K.A(b+1),BK.at(b-1),BK.at(b+2));
            auto spec = res.svdBond(b,wfK,(ha==1?Fromleft:Fromright),PH,args);

            if(verbose)
                {
                printfln("    Trunc. err=%.1E, States kept=%s",
                         spec.truncerr(),
                         showm(linkInd(res,b)) );
                }

            if(ha == 1)
                {
                BK.at(b) = lwfK * dag(prime(res.A(b)));
                }
            else
                {
                BK.at(b+1) = rwfK * dag(prime(res.A(b+1)));
                }
            }
        }
    }
template
void
fitApplyMPO(Real fac,
            MPSt<ITensor> const& psi,
            MPOt<ITensor> const& K,
            MPSt<ITensor>& res, 
            Sweeps const&, 
            Args args);
template
void
fitApplyMPO(Real fac,
            MPSt<IQTensor> const& psi,
            MPOt<IQTensor> const& K,
            MPSt<IQTensor>& res, 
            Sweeps const&, 
            Args args);

template<class Tensor>
Real
fitApplyMPO(MPSt<Tensor> const& psiA, 
            Real mpofac,
            MPSt<Tensor> const& psiB,
            MPOt<Tensor> const& K,
            MPSt<Tensor>& res,
            Args const& args)
    {
    return fitApplyMPO(1.,psiA,mpofac,psiB,K,res,args);
    }
template
Real
fitApplyMPO(MPSt<ITensor> const& psiA, 
            Real mpofac,
            MPSt<ITensor> const& psiB,
            MPOt<ITensor> const& K,
            MPSt<ITensor> & res,
            Args const& args);
template
Real
fitApplyMPO(MPSt<IQTensor> const& psiA, 
            Real mpofac,
            MPSt<IQTensor> const& psiB,
            MPOt<IQTensor> const& K,
            MPSt<IQTensor> & res,
            Args const& args);

template<class Tensor>
Real
fitApplyMPO(Real mpsfac,
            MPSt<Tensor> const& psiA, 
            Real mpofac,
            MPSt<Tensor> const& psiB,
            MPOt<Tensor> const& K,
            MPSt<Tensor>& res,
            Args const& args)
    {
    if(&psiA == &res || &psiB == &res)
        {
        Error("fitApplyMPO: Result MPS cannot be same as an input MPS");
        }
    auto N = psiA.N();
    auto nsweep = args.getInt("Nsweep",1);

    res.position(1);

    vector<Tensor> B(N+2),
                   BK(N+2);

    B.at(N) = psiA.A(N)*dag(prime(res.A(N),Link));
    BK.at(N) = psiB.A(N)*K.A(N)*dag(prime(res.A(N)));
    for(int n = N-1; n > 2; --n)
        {
        B.at(n) = B.at(n+1)*psiA.A(n)*dag(prime(res.A(n),Link));
        BK.at(n) = BK.at(n+1)*psiB.A(n)*K.A(n)*dag(prime(res.A(n)));
        }


    for(int sw = 1; sw <= nsweep; ++sw)
        {
        for(int b = 1, ha = 1; ha <= 2; sweepnext(b,ha,N))
            {
            Tensor lwf = (B.at(b-1) ? B.at(b-1)*psiA.A(b) : psiA.A(b));
            Tensor rwf = (B.at(b+2) ? psiA.A(b+1)*B.at(b+2) : psiA.A(b+1));

            Tensor lwfK = (BK.at(b-1) ? BK.at(b-1)*psiB.A(b) : psiB.A(b));
            lwfK *= K.A(b);
            Tensor rwfK = (BK.at(b+2) ? BK.at(b+2)*psiB.A(b+1) : psiB.A(b+1));
            rwfK *= K.A(b+1);

            Tensor wf = mpsfac*noprime(lwf*rwf) + mpofac*noprime(lwfK*rwfK);
            wf.noprime();

            res.svdBond(b,wf,(ha==1?Fromleft:Fromright),args+Args("UseSVD",true));

            if(ha == 1)
                {
                B.at(b) = lwf * dag(prime(res.A(b),Link));
                BK.at(b) = lwfK * dag(prime(res.A(b)));
                }
            else
                {
                B.at(b+1) = rwf * dag(prime(res.A(b+1),Link));
                BK.at(b+1) = rwfK * dag(prime(res.A(b+1)));
                }
            }
        }

    auto olp = B.at(3);
    olp *= psiA.A(2);
    olp *= dag(prime(res.A(2),Link));
    olp *= psiA.A(1);
    olp *= dag(prime(res.A(1),Link));

    return olp.real();
    }
template
Real
fitApplyMPO(Real mpsfac,
            MPSt<ITensor> const& psiA, 
            Real mpofac,
            MPSt<ITensor> const& psiB,
            MPOt<ITensor> const& K,
            MPSt<ITensor>& res,
            Args const& args);
template
Real
fitApplyMPO(Real mpsfac,
            MPSt<IQTensor> const& psiA, 
            Real mpofac,
            MPSt<IQTensor> const& psiB,
            MPOt<IQTensor> const& K,
            MPSt<IQTensor> & res,
            Args const& args);

template<class Tensor>
void 
expsmallH(MPOt<Tensor> const& H, 
          MPOt<Tensor> & K, 
          Real tau, 
          Real Etot, 
          Real Kcutoff,
          Args args)
    {
    int ord = args.getInt("ExpHOrder",50);
    bool verbose = args.getBool("Verbose",false);
    args.add("Cutoff",MIN_CUT);
    args.add("Maxm",MAX_M);

    MPOt<Tensor> Hshift(H.sites());
    Hshift.Anc(1) *= -Etot;
    Hshift.plusEq(H,args);
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

        K = sum(xx,args);
        if(o > 1)
            nmultMPO(K,Hshift,xx[1],args);
        }
    if(verbose) cout << endl;
    }
template
void 
expsmallH(MPO const& H, MPO & K, Real tau, Real Etot, Real Kcutoff, Args args);
template
void 
expsmallH(IQMPO const& H, IQMPO & K, Real tau, Real Etot, Real Kcutoff, Args args);

template<class Tensor>
void 
expH(MPOt<Tensor> const& H, 
     MPOt<Tensor>& K, 
     Real tau, 
     Real Etot,
     Real Kcutoff, 
     int ndoub, 
     Args args)
    {
    const bool verbose = args.getBool("Verbose",false);
    Real ttau = tau / pow(2.0,ndoub);
    //cout << "ttau in expH is " << ttau << endl;

    Real smallcut = 0.1*Kcutoff*pow(0.25,ndoub);
    expsmallH(H, K, ttau,Etot,smallcut,args);

    if(verbose) cout << "Starting doubling in expH" << endl;
    for(int doub = 1; doub <= ndoub; ++doub)
        {
        //cout << " Double step " << doub << endl;
        if(doub == ndoub) 
            args.add("Cutoff",Kcutoff);
        else
            args.add("Cutoff",0.1 * Kcutoff * pow(0.25,ndoub-doub));
        //cout << "in expH, K.cutoff is " << K.cutoff << endl;
        MPOt<Tensor> KK;
        nmultMPO(K,K,KK,args);
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
expH(const MPO& H, MPO& K, Real tau, Real Etot,Real Kcutoff, int ndoub, Args);
template
void 
expH(const IQMPO& H, IQMPO& K, Real tau, Real Etot,Real Kcutoff, int ndoub, Args);

template<class Tensor>
void
applyExpH(MPSt<Tensor> const& psi, 
          MPOt<Tensor> const& H, 
          Real tau, 
          MPSt<Tensor>& res, 
          Args const& args)
    {
    //using IndexT = typename Tensor::index_type;
    using MPST = MPSt<Tensor>;

    if(&psi == &res) Error("Must pass distinct MPS arguments to applyExpH");

    const int order = args.getInt("Order",10);

    const int N = res.N();
    const int nsweep = args.getInt("Nsweep",1);

    res.position(1);

    vector<Tensor> lastB(N+2),
                   B(N+2),
                   BH(N+2);

    B.at(N) = psi.A(N)*dag(prime(psi.A(N),Link));
    BH.at(N) = psi.A(N)*H.A(N)*dag(prime(psi.A(N)));
    for(int n = N-1; n > 2; --n)
        {
        B.at(n) = B.at(n+1)*psi.A(n)*dag(prime(psi.A(n),Link));
        BH.at(n) = BH.at(n+1)*psi.A(n)*H.A(n)*dag(prime(psi.A(n)));
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
                    lwf = (B.at(b-1) ? B.at(b-1)*psi.A(b) : psi.A(b) );
                    rwf = (B.at(b+2) ? B.at(b+2)*psi.A(b+1) : psi.A(b+1));

                    lwfH = (BH.at(b-1) ? BH.at(b-1)*last.A(b) : last.A(b));
                    lwfH *= H.A(b);
                    rwfH = (BH.at(b+2) ? BH.at(b+2)*last.A(b+1) : last.A(b+1));
                    rwfH *= H.A(b+1);
                    }
                else //dn
                    {
                    lwf = (B.at(b-1) ? B.at(b-1)*dag(prime(psi.A(b),Link)) : dag(prime(psi.A(b),Link)));
                    rwf = (B.at(b+2) ? B.at(b+2)*dag(prime(psi.A(b+1),Link)) : dag(prime(psi.A(b+1),Link)));

                    lwfH = (BH.at(b-1) ? BH.at(b-1)*dag(prime(last.A(b))) : dag(prime(last.A(b))));
                    lwfH *= H.A(b);
                    rwfH = (BH.at(b+2) ? BH.at(b+2)*dag(prime(last.A(b+1))) : dag(prime(last.A(b+1))));
                    rwfH *= H.A(b+1);
                    }

                auto wf = noprime(lwf*rwf) + mpofac*noprime(lwfH*rwfH);
                if(!up) wf.dag();

                res.svdBond(b,wf,(ha==1?Fromleft:Fromright),args+Args("UseSVD",true));

                if(up)
                    {
                    if(ha == 1)
                        {
                        B.at(b) = lwf * dag(prime(res.A(b),Link));
                        BH.at(b) = lwfH * dag(prime(res.A(b)));
                        }
                    else
                        {
                        B.at(b+1) = rwf * dag(prime(res.A(b+1),Link));
                        BH.at(b+1) = rwfH * dag(prime(res.A(b+1)));
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
applyExpH(const MPSt<ITensor>& psi, const MPOt<ITensor>& H, Real tau, MPSt<ITensor>& res, const Args& args);
template
void
applyExpH(const MPSt<IQTensor>& psi, const MPOt<IQTensor>& H, Real tau, MPSt<IQTensor>& res, const Args& args);


//
// For now this is unsupported
//
template<class Tensor>
void 
zipUpApplyMPO(MPSt<Tensor> const& psi, 
              MPOt<Tensor> const& K, 
              MPSt<Tensor>& res, 
              Args const& args)
    {
    using IndexT = typename Tensor::index_type;

    const
    bool allow_arb_position = args.getBool("AllowArbPosition",false);

    if(&psi == &res)
        Error("psi and res must be different MPS instances");

    //Real cutoff = args.getReal("Cutoff",psi.cutoff());
    //int maxm = args.getInt("Maxm",psi.maxm());

    auto N = psi.N();
    if(K.N() != N) 
        Error("Mismatched N in zipUpApplyMPO");

    if(!itensor::isOrtho(psi) || itensor::orthoCenter(psi) != 1)
        Error("Ortho center of psi must be site 1");

    if(!allow_arb_position && (!itensor::isOrtho(K) || itensor::orthoCenter(K) != 1))
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
        denmatDecomp(clust, res.Anc(i), nfork,Fromleft,args);
        IndexT mid = commonIndex(res.A(i),nfork,Link);
        //assert(mid.dir() == In);
        mid.dag();
        midsize[i] = mid.m();
        maxdim = std::max(midsize[i],maxdim);
        assert(rightLinkInd(res,i+1).dir() == Out);
        res.Anc(i+1) = Tensor(mid,prime(res.sites()(i+1)),rightLinkInd(res,i+1));
        }
    nfork = clust * psi.A(N) * K.A(N);
    //if(nfork.iten_size() == 0)	// this product gives 0 !!
	//throw ResultIsZero("nfork.iten size == 0");

    res.svdBond(N-1,nfork,Fromright,args);
    res.noprimelink();
    res.mapprime(1,0,Site);
    res.position(1);
    } //void zipUpApplyMPO
template
void 
zipUpApplyMPO(const MPS& x, const MPO& K, MPS& res, const Args& args);
template
void 
zipUpApplyMPO(const IQMPS& x, const IQMPO& K, IQMPS& res, const Args& args);


} //namespace itensor
