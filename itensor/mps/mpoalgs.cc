//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
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
    if(!args.defined("RespectDegenerate")) args.add("RespectDegenerate",true);

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

    auto lA = linkInds(A);
    auto lB = linkInds(B);
    auto sA = uniqueSiteInds(A,B);
    auto sB = uniqueSiteInds(B,A);

    // Check that A and B have unique indices
    for(int i = 1; i <= N; ++i)
      {
      if(!sA(i)) throw ITError("Error in nmultMPO(A,B): MPO tensor A("+str(i)+") does not have a unique site index. You may have meant to call nmultMPO(A,prime(B)).");
      if(!sB(i)) throw ITError("Error in nmultMPO(A,B): MPO tensor B("+str(i)+") does not have a unique site index. You may have meant to call nmultMPO(A,prime(B)).");
      }

    res=A;
    res.ref(1) = ITensor(sA(1),sB(1),lA(1));

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

        nfork = ITensor(lA(i),lB(i),linkIndex(res,i));
        denmatDecomp(clust,res.ref(i),nfork,Fromleft,{args,"Tags=",tags(lA(i))});

        auto mid = commonIndex(res(i),nfork);
        mid.dag();
        res.ref(i+1) = ITensor(mid,sA(i+1),sB(i+1),rightLinkIndex(res,i+1));
        }

    nfork = clust * A(N) * B(N);

    res.svdBond(N-1,nfork,Fromright, args);
    res.orthogonalize();
    }

MPO
nmultMPO(MPO const& A,
         MPO const& B,
         Args args)
  {
  MPO res;
  nmultMPO(A,B,res,args);
  return res;
  }

//
// Define specific applyMPO methods
//

MPS
densityMatrixApplyMPOImpl(MPO const& K,
                          MPS const& x,
                          Args args = Args::global());

void
fitApplyMPOImpl(MPS const& psi,
                MPO const& K,
                MPS & res,
                Args const& args = Args::global());

MPS
applyMPO(MPO const& K,
         MPS const& x,
         Args args)
    {
    if( !x ) Error("Error in applyMPO, MPS is uninitialized.");
    if( !K ) Error("Error in applyMPO, MPO is uninitialized.");

    auto method = args.getString("Method","DensityMatrix");
    if(!args.defined("RespectDegenerate")) args.add("RespectDegenerate",true);

    MPS res;
    if(method == "DensityMatrix")
        {
        res = densityMatrixApplyMPOImpl(K,x,args);
        }
    else if(method == "Fit")
        {
        // Use the input MPS x to be applied as the
        // default starting state
        // TODO: consider using zipUpApplyMPOImpl as 
        // a way to get a better starting state
        auto sites = uniqueSiteInds(K,x);
        res = replaceSiteInds(x,sites);
        //res = x;
        fitApplyMPOImpl(x,K,res,args);
        }
    else
        {
        Error("applyMPO currently supports the following methods: 'DensityMatrix', 'Fit'");
        }

    return res;
    }


MPS
applyMPO(MPO const& K,
         MPS const& x,
         MPS const& x0,
         Args args)
    {
    if( !x ) Error("Error in applyMPO, MPS is uninitialized.");
    if( !K ) Error("Error in applyMPO, MPO is uninitialized.");
    if( !x0 ) Error("Error in applyMPO, guess MPS is uninitialized.");

    auto method = args.getString("Method","Fit");
    if(!args.defined("RespectDegenerate")) args.add("RespectDegenerate",true);

    MPS res = x0;
    if(method == "DensityMatrix")
        Error("applyMPO method 'DensityMatrix' does not accept an input MPS");
    else if(method == "Fit")
        fitApplyMPOImpl(x,K,res,args);
    else
        Error("applyMPO currently supports the following methods: 'DensityMatrix', 'Fit'");

    return res;
    }

//
// Implement specific applyMPO methods
//


MPS
densityMatrixApplyMPOImpl(MPO const& K,
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
    dargs.add("RespectDegenerate",args.getBool("RespectDegenerate",true));
    auto verbose = args.getBool("Verbose",false);
    auto normalize = args.getBool("Normalize",false);
    auto dowrite = args.defined("WriteDim") && (maxLinkDim(psi) >= args.getInt("WriteDim"));
    std::string writedir_ = "./";
    if(dowrite) writedir_ = mkTempDir("E",args.getString("WriteDir","./"));

    auto N = length(psi);

    for( auto n : range1(N) )
      {
      if( commonIndex(psi(n),K(n)) != siteIndex(psi,n) )
          Error("MPS and MPO have different site indices in applyMPO method 'DensityMatrix'");
      }

    auto rand_plev = 14741;

    auto res = psi;

    //Set up conjugate psi and K
    auto psic = psi;
    auto Kc = K;
    //TODO: use sim(linkInds), sim(siteInds)
    psic.prime(rand_plev);
    Kc.prime(rand_plev);

    // Make sure the original and conjugates match
    for(auto j : range1(N-1)) 
        Kc.ref(j).prime(-rand_plev,uniqueSiteIndex(Kc,psic,j));

    //Build environment tensors from the left
    if(verbose) print("Building environment tensors...");
    auto E = std::vector<ITensor>(N+1);
    E[1] = psi(1)*K(1)*dag(Kc(1))*dag(psic(1));
    for(int j = 2; j < N; ++j)
        {
        E[j] = E[j-1]*psi(j)*K(j)*dag(Kc(j))*dag(psic(j));
        if(dowrite)
            {
            writeToFile(format("%s/E_%03d",writedir_,j-1),E[j-1]);
            E[j-1] = ITensor();
            }
        }
    if(verbose) println("done");

    //O is the representation of the product of K*psi in the new MPS basis
    auto O = psi(N)*K(N);

    auto rho = E[N-1] * O * dag(prime(O,rand_plev));
    E[N-1] = ITensor();

    ITensor U,D;
    auto ts = tags(linkIndex(psi,N-1));
    auto spec = diagPosSemiDef(rho,U,D,{dargs,"Tags=",ts});
    if(verbose) printfln("  j=%02d truncerr=%.2E dim=%d",N-1,spec.truncerr(),dim(commonIndex(U,D)));

    res.ref(N) = dag(U);

    O = O*U*psi(N-1)*K(N-1);

    for(int j = N-1; j > 1; --j)
        {
        if(dowrite) readFromFile(format("%s/E_%03d",writedir_,j-1),E[j-1]);
        if(not maxdim_set)
            {
            //Infer maxdim from bond dim of original MPS
            //times bond dim of MPO
            //i.e. upper bound on order of rho
            auto cip = commonIndex(psi(j),E[j-1]);
            auto ciw = commonIndex(K(j),E[j-1]);
            auto maxdim = (cip) ? dim(cip) : 1l;
            maxdim *= (ciw) ? dim(ciw) : 1l;
            dargs.add("MaxDim",maxdim);
            }
        rho = E[j-1] * O * dag(prime(O,rand_plev));
        E[j-1] = ITensor();
        ts = tags(linkIndex(psi,j-1));
        auto spec = diagPosSemiDef(rho,U,D,{dargs,"Tags=",ts});
        O = O*U*psi(j-1)*K(j-1);
        res.ref(j) = dag(U);
        if(verbose) printfln("  j=%02d truncerr=%.2E dim=%d",j,spec.truncerr(),dim(commonIndex(U,D)));
        }

    if(normalize) O /= norm(O);
    res.ref(1) = O;
    res.leftLim(0);
    res.rightLim(2);

    if(dowrite)
        {
        const string cmdstr = "rm -fr " + writedir_;
        system(cmdstr.c_str());
        }
    return res;
    }

void
oneSiteFitApply(vector<ITensor> & E,
                Real fac,
                MPS const& x,
                MPO const& K,
                MPS & Kx,
                Args const& args)
    {
    auto N = length(x);
    auto verbose = args.getBool("Verbose",false);
    auto normalize = args.getBool("Normalize",false);
    auto sw = args.getInt("Sweep",1);

    for(int s = 1, ha = 1; ha <= 2; sweepnext1(s,ha,N))
        {
        if(verbose)
            {
            printfln("Sweep=%d, HS=%d, Site=%d",sw,ha,s);
            }
        auto ds = (ha==1 ? +1 : -1);

        auto nE = E[s-ds]*x(s)*K(s);
        auto P = nE*E[s+ds];

        if(order(P) > 3)
            {
            Print(P);
            Error("order > 3 of P");
            }

        P *= fac;
        if(normalize) P /= norm(P);

        if(s+ds >= 1 && s+ds <= N)
            {
            auto ci = commonIndex(Kx(s),Kx(s+ds));
            if(!hasIndex(P,ci))
                {
                Print(ci);
                Print(inds(P));
                Error("P does not have Index ci");
                }
            auto [U,S,V] = svd(P,{ci},args);
            (void)U;
            (void)S;
            Kx.ref(s) = dag(V);
            }
        else
            {
            Kx.ref(s) = dag(P);
            }

        // Update environment:
        E[s] = nE * Kx(s);
        }
    }

void
twoSiteFitApply(vector<ITensor> & E,
                Real fac,
                MPS const& x,
                MPO const& K,
                MPS & Kx,
                Args const& args)
    {
    auto N = length(x);
    auto verbose = args.getBool("Verbose",false);
    auto normalize = args.getBool("Normalize",false);
    auto sw = args.getInt("Sweep",1);

    for(int b = 1, ha = 1; ha <= 2; sweepnext(b,ha,N))
        {
        if(verbose)
            {
            printfln("Sweep=%d, HS=%d, Bond=(%d,%d)",sw,ha,b,b+1);
            }

        auto lwfK = E[b-1]*x(b);
        lwfK *= K(b);
        auto rwfK = E[b+2]*x(b+1);
        rwfK *= K(b+1);

        auto wfK = lwfK*rwfK;
        wfK *= fac;

        if(normalize) wfK /= norm(wfK);
        auto PH = LocalOp(K(b),K(b+1),E[b-1],E[b+2]);

        wfK.dag();
        auto spec = Kx.svdBond(b,wfK,(ha==1?Fromleft:Fromright),PH,args);

        if(verbose)
            {
            printfln("    Trunc. err=%.1E, States kept=%s",
                     spec.truncerr(),
                     showDim(linkIndex(Kx,b)) );
            }

        if(ha == 1)
            E[b] = lwfK * Kx(b);
        else
            E[b+1] = rwfK * Kx(b+1);
        }
    }

void
fitApplyMPOImpl(Real fac,
                MPS const& x,
                MPO const& K,
                MPS& Kx,
                Sweeps const& sweeps,
                Args args)
    {
    auto N = length(x);
    auto NCenterSites = args.getInt("NCenterSites",2);

    Kx.dag();
    // Make the indices of |Kx> and K|x> match
    Kx.replaceSiteInds(uniqueSiteInds(K,x));
    // Replace the link indices of |Kx> with similar ones
    // so they don't clash with the links of |x>
    Kx.replaceLinkInds(sim(linkInds(Kx)));
    Kx.position(1);

    auto E = vector<ITensor>(N+2,ITensor(1.));
    for(auto n = N; n > NCenterSites; --n)
        {
        E[n] = E[n+1]*x(n)*K(n)*Kx(n);
        }

    for(auto sw : range1(sweeps.nsweep()))
        {
        args.add("Sweep",sw);
        args.add("Cutoff",sweeps.cutoff(sw));
        args.add("MinDim",sweeps.mindim(sw));
        args.add("MaxDim",sweeps.maxdim(sw));
        args.add("Noise",sweeps.noise(sw));

        if(NCenterSites == 2)
            {
            twoSiteFitApply(E,fac,x,K,Kx,args);
            }
        else if(NCenterSites == 1)
            {
            oneSiteFitApply(E,fac,x,K,Kx,args);
            }
        else
            {
            Error(format("NCenterSites = %d not supported",NCenterSites));
            }
        }
    Kx.dag();
    }

void
fitApplyMPOImpl(Real fac,
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
    fitApplyMPOImpl(fac,psi,K,res,sweeps,args);
    }

void
fitApplyMPOImpl(MPS const& psi,
                MPO const& K,
                MPS& res,
                Args const& args)
    {
    fitApplyMPOImpl(1.,psi,K,res,args);
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
// Deprecated
//

//
// For now this is unsupported
// We can consider bringing it back for example
// as a way to get a default starting MPS
// for fitApplyMPOImpl
//

//void 
//zipUpApplyMPOImpl(MPS const& psi, 
//                  MPO const& K, 
//                  MPS& res, 
//                  Args const& args)
//    {
//    const
//    bool allow_arb_position = args.getBool("AllowArbPosition",false);
//
//    if(&psi == &res)
//        Error("psi and res must be different MPS instances");
//
//    auto N = length(psi);
//    if(length(K) != N) 
//        Error("Mismatched N in ApplyMPO() ZipUp method");
//
//    if(!itensor::isOrtho(psi) || itensor::orthoCenter(psi) != 1)
//        Error("Ortho center of psi must be site 1");
//
//    if(!allow_arb_position && (!itensor::isOrtho(K) || itensor::orthoCenter(K) != 1))
//        Error("Ortho center of K must be site 1");
//
//#ifdef DEBUG
//    checkQNs(psi);
//    checkQNs(K);
//    /*
//    cout << "Checking divergence in zip" << endl;
//    for(int i = 1; i <= N; i++)
//	div(psi(i));
//    for(int i = 1; i <= N; i++)
//	div(K(i));
//    cout << "Done Checking divergence in zip" << endl;
//    */
//#endif
//
//    res = psi; 
//    res.replaceTags("0","4","Link");
//    res.replaceTags("0","1","Site");
//
//    ITensor clust,nfork;
//    vector<int> midsize(N);
//    int maxdim = 1;
//    for(int i = 1; i < N; i++)
//        {
//        if(i == 1) { clust = psi(i) * K(i); }
//        else { clust = nfork * (psi(i) * K(i)); }
//        if(i == N-1) break; //No need to SVD for i == N-1
//
//        Index oldmid = linkIndex(res,i); assert(oldmid.dir() == Out);
//        nfork = ITensor(linkIndex(psi,i),linkIndex(K,i-1),oldmid);
//        //if(clust.iten_size() == 0)	// this product gives 0 !!
//	    //throw ResultIsZero("clust.iten size == 0");
//        denmatDecomp(clust, res.ref(i), nfork,Fromleft,args);
//        Index mid = commonIndex(res(i),nfork);
//        //assert(mid.dir() == In);
//        mid.dag();
//        midsize[i] = dim(mid);
//        maxdim = std::max(midsize[i],maxdim);
//        assert(linkIndex(res,i+1).dir() == Out);
//        res.ref(i+1) = ITensor(mid,prime(res.sites()(i+1)),linkIndex(res,i+1));
//        }
//    nfork = clust * psi(N) * K(N);
//    //if(nfork.iten_size() == 0)	// this product gives 0 !!
//    //throw ResultIsZero("nfork.iten size == 0");
//
//    res.svdBond(N-1,nfork,Fromright,args);
//    res.noPrime("Link");
//    res.replaceTags("1","0","Site");
//    res.position(1);
//    } //void zipUpApplyMPOImpl

//
// These versions calculate |res> = |psiA> + mpofac*H*|psiB>
// Currently they are unsupported

//Real
//fitApplyMPOImpl(MPS const& psiA, 
//                Real mpofac,
//                MPS const& psiB,
//                MPO const& K,
//                MPS& res,
//                Args const& args)
//    {
//    return fitApplyMPOImpl(1.,psiA,mpofac,psiB,K,res,args);
//    }
//
//
//Real
//fitApplyMPOImpl(Real mpsfac,
//                MPS const& psiA, 
//                Real mpofac,
//                MPS const& psiB,
//                MPO const& K,
//                MPS& res,
//                Args const& args)
//    {
//    if(&psiA == &res || &psiB == &res)
//        {
//        Error("fitApplyMPOImpl: Result MPS cannot be same as an input MPS");
//        }
//    auto N = length(psiA);
//    auto nsweep = args.getInt("Nsweep",1);
//
//    res.position(1);
//
//    vector<ITensor> B(N+2),
//                   E(N+2);
//
//    B.at(N) = psiA(N)*dag(prime(res(N),"Link"));
//    E.at(N) = psiB(N)*K(N)*dag(prime(res(N)));
//    for(int n = N-1; n > 2; --n)
//        {
//        B.at(n) = B.at(n+1)*psiA(n)*dag(prime(res(n),"Link"));
//        E.at(n) = E.at(n+1)*psiB(n)*K(n)*dag(prime(res(n)));
//        }
//
//
//    for(int sw = 1; sw <= nsweep; ++sw)
//        {
//        for(int b = 1, ha = 1; ha <= 2; sweepnext(b,ha,N))
//            {
//            ITensor lwf = (B.at(b-1) ? B.at(b-1)*psiA(b) : psiA(b));
//            ITensor rwf = (B.at(b+2) ? psiA(b+1)*B.at(b+2) : psiA(b+1));
//
//            ITensor lwfK = (E.at(b-1) ? E.at(b-1)*psiB(b) : psiB(b));
//            lwfK *= K(b);
//            ITensor rwfK = (E.at(b+2) ? E.at(b+2)*psiB(b+1) : psiB(b+1));
//            rwfK *= K(b+1);
//
//            ITensor wf = mpsfac*noPrime(lwf*rwf) + mpofac*noPrime(lwfK*rwfK);
//            wf.noPrime();
//
//            res.svdBond(b,wf,(ha==1?Fromleft:Fromright),args+Args("UseSVD",true));
//
//            if(ha == 1)
//                {
//                B.at(b) = lwf * dag(prime(res(b),"Link"));
//                E.at(b) = lwfK * dag(prime(res(b)));
//                }
//            else
//                {
//                B.at(b+1) = rwf * dag(prime(res(b+1),"Link"));
//                E.at(b+1) = rwfK * dag(prime(res(b+1)));
//                }
//            }
//        }
//
//    auto olp = B.at(3);
//    olp *= psiA(2);
//    olp *= dag(prime(res(2),"Link"));
//    olp *= psiA(1);
//    olp *= dag(prime(res(1),"Link"));
//
//    return olp.elt();
//    }

} //namespace itensor
