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

void 
nmultMPO(MPO const& Aorig, 
         MPO const& Borig, 
         MPO& res,
         Args args)
    {
    if(!args.defined("Cutoff")) args.add("Cutoff",1E-14);

    if(length(Aorig) != length(Borig)) Error("nmultMPO(MPO): Mismatched MPO length");
    const int N = length(Borig);

    auto A = Aorig;
    A.position(1);

    MPO B;
    if(&Borig == &Aorig)
        {
        B = A;
        }
    else
        {
        B = Borig;
        B.position(1);
        }

    B.prime();

    res=A;
    auto siA = uniqueIndex(A(1),B(1),A(2));
    auto siB = uniqueIndex(B(1),A(1),B(2));
    res.ref(1) = ITensor(siA,siB,linkIndex(A,1));

    //Print(A);
    //Print(B);

    ITensor clust,nfork;
    for(int i = 1; i < N; ++i)
        {
        if(i == 1) 
            { 
            clust = A(i) * B(i); 
            }
        else       
            { 
            clust = nfork * A(i) * B(i); 
            }

        if(i == N-1) break;

        nfork = ITensor(linkIndex(A,i),linkIndex(B,i),linkIndex(res,i));

        denmatDecomp(clust,res.ref(i),nfork,Fromleft,args);

        auto mid = commonIndex(res(i),nfork,"Link");
        mid.dag();
        auto siA = uniqueIndex(A(i+1),A(i),A(i+2),B(i+1));
        auto siB = uniqueIndex(B(i+1),B(i),B(i+2),A(i+1));
        res.ref(i+1) = ITensor(mid,siA,siB,rightLinkIndex(res,i+1));
        }

    nfork = clust * A(N) * B(N);

    res.svdBond(N-1,nfork,Fromright, args);
    for(auto i : range1(N))
        {
        if(i < N)
            {
            auto l = linkIndex(res,i);
            res.ref(i).noPrime(l);
            res.ref(i+1).noPrime(l);
            }
        res.ref(i).replaceTags("2","1");
        }
    res.orthogonalize();

    }


MPS
applyMPO(MPO const& K,
         MPS const& x,
         Args const& args)
    {
    if(not x(1).store()) Error("Error in applyMPO, MPS is uninitialized.");
    if(not K(1).store()) Error("Error in applyMPO, MPO is uninitialized.");
    auto method = args.getString("Method","DensityMatrix");

    //This is done here because fitApplyMPO() has a different default behavior
    //(for backwards compatability)
    auto normalize = args.getBool("Normalize",false);
    auto argsp = args;
    argsp.add("Normalize=",normalize);

    MPS res;
    if(method == "DensityMatrix")
        res = exactApplyMPO(K,x,argsp);
    else if(method == "Fit")
        res = fitApplyMPO(x,K,argsp);
    else
        Error("applyMPO currently supports the following methods: 'DensityMatrix' (previously called with exactApplyMPO), 'Fit' (previously called with fitApplyMPO)");

    return res;
    }


MPS
applyMPO(MPO const& K,
         MPS const& x,
         MPS const& x0,
         Args const& args)
    {
    if(not x(1).store()) Error("Error in applyMPO, MPS is uninitialized.");
    if(not K(1).store()) Error("Error in applyMPO, MPO is uninitialized.");
    if(not x0(1).store()) Error("Error in applyMPO, guess MPS is uninitialized.");
    auto method = args.getString("Method","Fit");

    //This is done here because fitApplyMPO() has a different default behavior
    //(for backwards compatability)
    auto normalize = args.getBool("Normalize",false);
    auto argsp = args;
    argsp.add("Normalize=",normalize);

    MPS res = x0;
    if(method == "DensityMatrix")
        Error("applyMPO method 'DensityMatrix' does not accept an input MPS");
    else if(method == "Fit")
        fitApplyMPO(x,K,res,argsp);
    else
        Error("applyMPO currently supports the following methods: 'DensityMatrix' (previously called with exactApplyMPO), 'Fit' (previously called with fitApplyMPO)");

    return res;
    }


MPS
exactApplyMPO(MPO const& K,
              MPS const& psi,
              Args args)
    {
    if( args.defined("Maxm") )
      {
      if( args.defined("MaxDim") )
        {
        Global::warnDeprecated("Args Maxm and MaxDim are both defined. Maxm is deprecated in favor of MaxDim, MaxDim will be used.");
        }
      else
        {
        Global::warnDeprecated("Arg Maxm is deprecated in favor of MaxDim.");
        args.add("MaxDim",args.getInt("Maxm"));
        }
      }

    auto cutoff = args.getReal("Cutoff",1E-13);
    auto dargs = Args{"Cutoff",cutoff};
    auto maxdim_set = args.defined("MaxDim");
    if(maxdim_set) dargs.add("MaxDim",args.getInt("MaxDim"));
    auto verbose = args.getBool("Verbose",false);
    auto normalize = args.getBool("Normalize",false);
    //auto siteType = getIndexType(args,"SiteType",Site);
    //auto linkType = getIndexType(args,"LinkType",Link);
    auto siteTags = getTagSet(args,"SiteTags","Site");
    auto linkTags = getTagSet(args,"LinkTags","Link");

    if(not commonIndex(K(1),psi(1),siteTags))
        Error("MPS and MPO have different site indices in applyMPO method 'DensityMatrix'");

    auto plev = 14741;
    auto plevtag = format("%d",plev);

    auto res = psi;

    auto N = length(psi);

    //Set up conjugate psi and K
    auto psic = psi;
    auto Kc = K;
    for(auto j : range1(N)) 
        {
        //Modify prime levels of psic and Kc
        if(j == 1)
            {
            auto ci = commonIndex(psi(1),psi(2),linkTags);
            psic.ref(j) = dag(prime(replaceTags(psi(j),"0","2",siteTags),plev,ci));
            ci = commonIndex(Kc(1),Kc(2),linkTags);
            Kc.ref(j) = dag(prime(replaceTags(K(j),"0","2",siteTags),plev,ci));
            }
        else
            {
            psic.ref(j) = dag(replaceTags(replaceTags(psi(j),"0","2",siteTags),"0",plevtag,linkTags));
            Kc.ref(j) = dag(replaceTags(replaceTags(K(j),"0","2",siteTags),"0",plevtag,linkTags));
            }
        }

    //Build environment tensors from the left
    if(verbose) print("Building environment tensors...");
    auto E = std::vector<ITensor>(N+1);
    E.at(1) = psi(1)*K(1)*Kc(1)*psic(1);
    for(int j = 2; j < N; ++j)
        {
        E.at(j) = E.at(j-1)*psi(j)*K(j)*Kc(j)*psic(j);
        //assert(order(E[j])==4);
        }
    if(verbose) println("done");

    //O is the representation of the product of K*psi in the new MPS basis
    auto O = psi(N)*K(N);
    O.noPrime(siteTags);

    auto rho = E.at(N-1) * O * dag(prime(O,plev));
    ITensor U,D;
    dargs.add("Tags=",format("Link,l=%d",N-1));
    auto spec = diagHermitian(rho,U,D,dargs);
    if(verbose) printfln("  j=%02d truncerr=%.2E m=%d",N-1,spec.truncerr(),dim(commonIndex(U,D)));

    res.ref(N) = dag(U);

    O = O*U*psi(N-1)*K(N-1);
    O.noPrime(siteTags);

    for(int j = N-1; j > 1; --j)
        {
        if(not maxdim_set)
            {
            //Infer maxdim from bond dim of original MPS
            //times bond dim of MPO
            //i.e. upper bound on order of rho
            auto cip = commonIndex(psi(j),E.at(j-1));
            auto ciw = commonIndex(K(j),E.at(j-1));
            auto maxdim = (cip) ? dim(cip) : 1l;
            maxdim *= (ciw) ? dim(ciw) : 1l;
            dargs.add("MaxDim",maxdim);
            }
        rho = E.at(j-1) * O * dag(prime(O,plev));
        //TODO: make sure this tag convention is working
        dargs.add("Tags=",format("Link,l=%d",j-1));
        auto spec = diagHermitian(rho,U,D,dargs);
        O = O*U*psi(j-1)*K(j-1);
        O.noPrime(siteTags);
        res.ref(j) = dag(U);
        if(verbose) printfln("  j=%02d truncerr=%.2E m=%d",j,spec.truncerr(),dim(commonIndex(U,D)));
        }

    if(normalize) O /= norm(O);
    res.ref(1) = O;
    res.leftLim(0);
    res.rightLim(2);

    return res;
    }


MPS
exactApplyMPO(MPS const& x,
              MPO const& K,
              Args const& args)
    {
    return exactApplyMPO(K,x,args);
    }


void 
exactApplyMPO(MPS const& x, 
              MPO const& K, 
              MPS      & res,
              Args const& args)
    {
    res = exactApplyMPO(K,x);
    }


MPS
fitApplyMPO(MPS const& psi,
            MPO const& K,
            Args const& args)
    {
    MPS res;
    fitApplyMPO(1.,psi,K,res,args);
    return res;
    }


void
fitApplyMPO(MPS const& psi,
            MPO const& K,
            MPS& res,
            Args const& args)
    {
    fitApplyMPO(1.,psi,K,res,args);
    }


void
fitApplyMPO(Real fac,
            MPS const& psi,
            MPO const& K,
            MPS& res,
            Args args)
    {
    if( args.defined("Maxm") )
      {
      if( args.defined("MaxDim") )
        {
        Global::warnDeprecated("Args Maxm and MaxDim are both defined. Maxm is deprecated in favor of MaxDim, MaxDim will be used.");
        }
      else
        {
        Global::warnDeprecated("Arg Maxm is deprecated in favor of MaxDim.");
        args.add("MaxDim",args.getInt("Maxm"));
        }
      }

    auto nsweep = args.getInt("Nsweep",1);
    Sweeps sweeps(nsweep);
    auto cutoff = args.getReal("Cutoff",-1);
    if(cutoff >= 0) sweeps.cutoff() = cutoff;
    auto maxdim = args.getInt("MaxDim",-1);
    if(maxdim >= 1) sweeps.maxdim() = maxdim;
    fitApplyMPO(fac,psi,K,res,sweeps,args);
    }

void
fitApplyMPO(Real fac,
            MPS const& psi,
            MPO const& K,
            MPS& res,
            Sweeps const& sweeps,
            Args args)
    {
    auto N = length(psi);
    auto verbose = args.getBool("Verbose",false);
    auto normalize = args.getBool("Normalize",true);

    const auto origPsi = psi;

    if(not res) res = origPsi;
    res.position(1);

    auto BK = vector<ITensor>(N+2);
    BK.at(N) = origPsi(N)*K(N)*dag(prime(res(N)));
    for(auto n = N-1; n > 2; --n)
        {
        BK.at(n) = BK.at(n+1)*origPsi(n)*K(n)*dag(prime(res(n)));
        }

    for(auto sw : range1(sweeps.nsweep()))
        {
        args.add("Sweep",sw);
        args.add("Cutoff",sweeps.cutoff(sw));
        args.add("MinDim",sweeps.mindim(sw));
        args.add("MaxDim",sweeps.maxdim(sw));
        args.add("Noise",sweeps.noise(sw));

        for(int b = 1, ha = 1; ha <= 2; sweepnext(b,ha,N))
            {
            if(verbose)
                {
                printfln("Sweep=%d, HS=%d, Bond=(%d,%d)",sw,ha,b,b+1);
                }

            //TODO: does this tag the correct bond, independent of the sweep direction?
            args.add("Tags",format("Link,l=%d",b));

            auto lwfK = (BK.at(b-1) ? BK.at(b-1)*origPsi(b) : origPsi(b));
            lwfK *= K(b);
            auto rwfK = (BK.at(b+2) ? BK.at(b+2)*origPsi(b+1) : origPsi(b+1));
            rwfK *= K(b+1);

            auto wfK = lwfK*rwfK;
            wfK.noPrime();
            wfK *= fac;

            if(normalize) wfK /= norm(wfK);
            auto PH = LocalOp(K(b),K(b+1),BK.at(b-1),BK.at(b+2));
            auto spec = res.svdBond(b,wfK,(ha==1?Fromleft:Fromright),PH,args);

            if(verbose)
                {
                printfln("    Trunc. err=%.1E, States kept=%s",
                         spec.truncerr(),
                         showDim(linkIndex(res,b)) );
                }

            if(ha == 1)
                {
                BK.at(b) = lwfK * dag(prime(res(b)));
                }
            else
                {
                BK.at(b+1) = rwfK * dag(prime(res(b+1)));
                }
            }
        }
    }


Real
fitApplyMPO(MPS const& psiA, 
            Real mpofac,
            MPS const& psiB,
            MPO const& K,
            MPS& res,
            Args const& args)
    {
    return fitApplyMPO(1.,psiA,mpofac,psiB,K,res,args);
    }


Real
fitApplyMPO(Real mpsfac,
            MPS const& psiA, 
            Real mpofac,
            MPS const& psiB,
            MPO const& K,
            MPS& res,
            Args const& args)
    {
    if(&psiA == &res || &psiB == &res)
        {
        Error("fitApplyMPO: Result MPS cannot be same as an input MPS");
        }
    auto N = length(psiA);
    auto nsweep = args.getInt("Nsweep",1);

    res.position(1);

    vector<ITensor> B(N+2),
                   BK(N+2);

    B.at(N) = psiA(N)*dag(prime(res(N),"Link"));
    BK.at(N) = psiB(N)*K(N)*dag(prime(res(N)));
    for(int n = N-1; n > 2; --n)
        {
        B.at(n) = B.at(n+1)*psiA(n)*dag(prime(res(n),"Link"));
        BK.at(n) = BK.at(n+1)*psiB(n)*K(n)*dag(prime(res(n)));
        }


    for(int sw = 1; sw <= nsweep; ++sw)
        {
        for(int b = 1, ha = 1; ha <= 2; sweepnext(b,ha,N))
            {
            ITensor lwf = (B.at(b-1) ? B.at(b-1)*psiA(b) : psiA(b));
            ITensor rwf = (B.at(b+2) ? psiA(b+1)*B.at(b+2) : psiA(b+1));

            ITensor lwfK = (BK.at(b-1) ? BK.at(b-1)*psiB(b) : psiB(b));
            lwfK *= K(b);
            ITensor rwfK = (BK.at(b+2) ? BK.at(b+2)*psiB(b+1) : psiB(b+1));
            rwfK *= K(b+1);

            ITensor wf = mpsfac*noPrime(lwf*rwf) + mpofac*noPrime(lwfK*rwfK);
            wf.noPrime();

            res.svdBond(b,wf,(ha==1?Fromleft:Fromright),args+Args("UseSVD",true));

            if(ha == 1)
                {
                B.at(b) = lwf * dag(prime(res(b),"Link"));
                BK.at(b) = lwfK * dag(prime(res(b)));
                }
            else
                {
                B.at(b+1) = rwf * dag(prime(res(b+1),"Link"));
                BK.at(b+1) = rwfK * dag(prime(res(b+1)));
                }
            }
        }

    auto olp = B.at(3);
    olp *= psiA(2);
    olp *= dag(prime(res(2),"Link"));
    olp *= psiA(1);
    olp *= dag(prime(res(1),"Link"));

    return olp.elt();
    }

void
applyExpH(MPS const& psi, 
          MPO const& H, 
          Real tau, 
          MPS& res, 
          Args const& args)
    {
    if(&psi == &res) Error("Must pass distinct MPS arguments to applyExpH");

    const int order = args.getInt("Order",10);

    const int N = length(res);
    const int nsweep = args.getInt("Nsweep",1);

    res.position(1);

    vector<ITensor> lastB(N+2),
                   B(N+2),
                   BH(N+2);

    B.at(N) = psi(N)*dag(prime(psi(N),"Link"));
    BH.at(N) = psi(N)*H(N)*dag(prime(psi(N)));
    for(int n = N-1; n > 2; --n)
        {
        B.at(n) = B.at(n+1)*psi(n)*dag(prime(psi(n),"Link"));
        BH.at(n) = BH.at(n+1)*psi(n)*H(n)*dag(prime(psi(n)));
        }

    lastB = B;

    MPS last(psi);

    bool up = true;

    for(int ord = order, n = 0; ord >= 1; --ord, ++n)
        {
        const Real mpofac = -tau/(1.*ord);

        if(n > 0) lastB.swap(B);

        for(int sw = 1; sw <= nsweep; ++sw)
            {
            for(int b = 1, ha = 1; ha <= 2; sweepnext(b,ha,N))
                {
                ITensor lwf,rwf,
                       lwfH,rwfH;

                if(up)
                    {
                    lwf = (B.at(b-1) ? B.at(b-1)*psi(b) : psi(b) );
                    rwf = (B.at(b+2) ? B.at(b+2)*psi(b+1) : psi(b+1));

                    lwfH = (BH.at(b-1) ? BH.at(b-1)*last(b) : last(b));
                    lwfH *= H(b);
                    rwfH = (BH.at(b+2) ? BH.at(b+2)*last(b+1) : last(b+1));
                    rwfH *= H(b+1);
                    }
                else //dn
                    {
                    lwf = (B.at(b-1) ? B.at(b-1)*dag(prime(psi(b),"Link")) : dag(prime(psi(b),"Link")));
                    rwf = (B.at(b+2) ? B.at(b+2)*dag(prime(psi(b+1),"Link")) : dag(prime(psi(b+1),"Link")));

                    lwfH = (BH.at(b-1) ? BH.at(b-1)*dag(prime(last(b))) : dag(prime(last(b))));
                    lwfH *= H(b);
                    rwfH = (BH.at(b+2) ? BH.at(b+2)*dag(prime(last(b+1))) : dag(prime(last(b+1))));
                    rwfH *= H(b+1);
                    }

                auto wf = noPrime(lwf*rwf) + mpofac*noPrime(lwfH*rwfH);
                if(!up) wf.dag();

                res.svdBond(b,wf,(ha==1?Fromleft:Fromright),args+Args("UseSVD",true));

                if(up)
                    {
                    if(ha == 1)
                        {
                        B.at(b) = lwf * dag(prime(res(b),"Link"));
                        BH.at(b) = lwfH * dag(prime(res(b)));
                        }
                    else
                        {
                        B.at(b+1) = rwf * dag(prime(res(b+1),"Link"));
                        BH.at(b+1) = rwfH * dag(prime(res(b+1)));
                        }
                    }
                else //dn
                    {
                    if(ha == 1)
                        {
                        B.at(b) = lwf * res(b);
                        BH.at(b) = lwfH * res(b);
                        }
                    else
                        {
                        B.at(b+1) = rwf * res(b+1);
                        BH.at(b+1) = rwfH * res(b+1);
                        }
                    }
                }
            }

        last = res;

        up = !up;

        } // for ord

    }

//
// For now this is unsupported
//
void 
zipUpApplyMPO(MPS const& psi, 
              MPO const& K, 
              MPS& res, 
              Args const& args)
    {
    const
    bool allow_arb_position = args.getBool("AllowArbPosition",false);

    if(&psi == &res)
        Error("psi and res must be different MPS instances");

    auto N = length(psi);
    if(length(K) != N) 
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
	div(psi(i));
    for(int i = 1; i <= N; i++)
	div(K(i));
    cout << "Done Checking divergence in zip" << endl;
    */
#endif

    res = psi; 
    res.replaceTags("0","4","Link");
    res.replaceTags("0","1","Site");

    ITensor clust,nfork;
    vector<int> midsize(N);
    int maxdim = 1;
    for(int i = 1; i < N; i++)
        {
        if(i == 1) { clust = psi(i) * K(i); }
        else { clust = nfork * (psi(i) * K(i)); }
        if(i == N-1) break; //No need to SVD for i == N-1

        Index oldmid = rightLinkIndex(res,i); assert(oldmid.dir() == Out);
        nfork = ITensor(rightLinkIndex(psi,i),rightLinkIndex(K,i),oldmid);
        //if(clust.iten_size() == 0)	// this product gives 0 !!
	    //throw ResultIsZero("clust.iten size == 0");
        denmatDecomp(clust, res.ref(i), nfork,Fromleft,args);
        Index mid = commonIndex(res(i),nfork);
        //assert(mid.dir() == In);
        mid.dag();
        midsize[i] = dim(mid);
        maxdim = std::max(midsize[i],maxdim);
        assert(rightLinkIndex(res,i+1).dir() == Out);
        res.ref(i+1) = ITensor(mid,prime(res.sites()(i+1)),rightLinkIndex(res,i+1));
        }
    nfork = clust * psi(N) * K(N);
    //if(nfork.iten_size() == 0)	// this product gives 0 !!
    //throw ResultIsZero("nfork.iten size == 0");

    res.svdBond(N-1,nfork,Fromright,args);
    res.noPrime("Link");
    res.replaceTags("1","0","Site");
    res.position(1);
    } //void zipUpApplyMPO

} //namespace itensor
