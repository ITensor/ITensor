#include "itensor/all.h"
using namespace itensor;

using SHalfSOne = MixedSiteSet<SpinHalfSite,SpinOneSite>;

int 
main()
    {
    int N = 100;
    auto Jhh = 0.5;
    auto Joo = 0.5;
    auto Jho = 1.0;

    auto sites = SHalfSOne(N,{"ConserveQNs=",true});

    //
    // Use the AutoMPO feature to create the 
    // next-neighbor Heisenberg model
    //
    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N; j += 1)
        {
        ampo += Jho/2,"S+",j,"S-",j+1;
        ampo += Jho/2,"S-",j,"S+",j+1;
        ampo +=   Jho,"Sz",j,"Sz",j+1;
        }
    for(int j = 1; j < N-1; j += 2)
        {
        ampo += Jhh/2,"S+",j,"S-",j+2;
        ampo += Jhh/2,"S-",j,"S+",j+2;
        ampo +=   Jhh,"Sz",j,"Sz",j+2;
        }
    for(int j = 2; j < N-1; j += 2)
        {
        ampo += Joo/2,"S+",j,"S-",j+2;
        ampo += Joo/2,"S-",j,"S+",j+2;
        ampo +=   Joo,"Sz",j,"Sz",j+2;
        }
    auto H = toMPO(ampo);

    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
    auto state = InitState(sites);
    for(int i = 1, n = 1; i <= N; i += 2, n += 1)
        {
        if(n%2 == 1) state.set(i,"Up");
        else         state.set(i,"Dn");
        }
    for(int i = 2, n = 1; i <= N; i += 2, n += 1)
        {
        if(n%2 == 1) state.set(i,"Up");
        else         state.set(i,"Dn");
        }
    auto psi = MPS(state);

    //
    // inner calculates matrix elements of MPO's with respect to MPS's
    // inner(psi,H,psi) = <psi|H|psi>
    //
    printfln("Initial energy = %.5f", inner(psi,H,psi) );

    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. 
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).
    //
    auto sweeps = Sweeps(10);
    sweeps.maxdim() = 10,10,20,40,80,100,140,180,200;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;
    println(sweeps);

    //
    // Begin the DMRG calculation
    //
    auto energy = dmrg(psi,H,sweeps,"Quiet");

    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",energy);
    printfln("\nUsing inner = %.10f", inner(psi,H,psi) );

    return 0;
    }
