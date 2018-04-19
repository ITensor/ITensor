#include "itensor/all.h"
#include "TStateObserver.h"
#include "S2.h"

using namespace std;
using namespace itensor;

int 
main(int argc, char* argv[])
    {
    //Get parameter file
    if(argc != 2)
        {
        printfln("Usage: %s inputfile.",argv[0]);
        return 0;
        }
    auto input = InputGroup(argv[1],"input");

    auto N = input.getInt("N",10);

    auto beta = input.getReal("beta",1);
    auto tau = input.getReal("tau",0.01);

    auto maxm = input.getInt("maxm",1000);
    auto cutoff = input.getReal("cutoff",1E-11);

    auto verbose = input.getYesNo("verbose",false);

    Args args;
    args.add("Maxm",maxm);
    args.add("Cutoff",cutoff);
    args.add("Verbose",verbose);

    auto sites = SpinHalf(2*N);

    auto ampo = AutoMPO(sites);
    for(auto j : range1(N-1))
        {
        auto s1 = 2*j-1,
             s2 = 2*j+1;
        ampo += 0.5,"S+",s1,"S-",s2;
        ampo += 0.5,"S-",s1,"S+",s2;
        ampo +=     "Sz",s1,"Sz",s2;
        }

    auto H = MPO(ampo);

    auto expH = toExpH<ITensor>(ampo,tau);

    auto S2 = toMPO(makeS2(sites,{"SkipAncilla=",true}));

    //
    // Make initial 'wavefunction' which is a product
    // of perfect singlets between neighboring sites
    //
    auto psi = MPS(sites);
    for(int n = 1; n <= 2*N; n += 2)
        {
        auto s1 = sites(n);
        auto s2 = sites(n+1);
        auto wf = ITensor(s1,s2);
        wf.set(s1(1),s2(2), ISqrt2);
        wf.set(s1(2),s2(1), -ISqrt2);
        ITensor D;
        psi.Aref(n) = ITensor(s1);
        psi.Aref(n+1) = ITensor(s2);
        svd(wf,psi.Aref(n),D,psi.Aref(n+1));
        psi.Aref(n) *= D;
        }

    auto obs = TStateObserver<ITensor>(psi);

    auto ttotal = beta/2.;
    const int nt = int(ttotal/tau+(1e-9*(ttotal/tau)));
    if(fabs(nt*tau-ttotal) > 1E-9)
        {
        Error("Timestep not commensurate with total time");
        }
    printfln("Doing %d steps of tau=%f",nt,tau);

    auto targs = args;

    auto En = Vector(nt);
    auto Sus = Vector(nt);
    auto Betas = Vector(nt);

    Real tsofar = 0;
    for(int tt = 1; tt <= nt; ++tt)
        {
        psi = exactApplyMPO(expH,psi,args);

        //Normalize wavefunction
        psi.Aref(1) /= norm(psi.A(1));

        tsofar += tau;
        targs.add("TimeStepNum",tt);
        targs.add("Time",tsofar);
        targs.add("TotalTime",ttotal);
        obs.measure(targs);

        //Record beta value
        auto bb = (2*tsofar);
        Betas(tt-1) = bb;

        //
        // Measure Energy
        //
        auto en = overlap(psi,H,psi);
        printfln("\nEnergy/N %.4f %.20f",bb,en/N);
        En(tt-1) = en/N;

        //
        // Measure Susceptibility
        //
        auto s2val = overlap(psi,S2,psi);
        Sus(tt-1) = (s2val*bb/3.)/N;

        println();
        }

    std::ofstream enf("en.dat");
    std::ofstream susf("sus.dat");
    for(auto n : range(Betas))
        {
        enf << format("%.14f %.14f\n",Betas(n),En(n));
        susf << format("%.14f %.14f\n",Betas(n),Sus(n));
        }
    enf.close();
    susf.close();

    return 0;
    }
