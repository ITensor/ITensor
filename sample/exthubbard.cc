#include "itensor/util/print_macro.h"
#include "itensor/mps/dmrg.h"
#include "itensor/mps/sites/hubbard.h"
#include "itensor/mps/autompo.h"

using namespace itensor;

int main(int argc, char* argv[])
    {
    //Parse the input file
    if(argc != 2) { printfln("Usage: %s inputfile_exthubbard",argv[0]); return 0; }
    InputGroup basic(argv[1],"basic");

    auto N = basic.getInt("N");
    auto Npart = basic.getInt("Npart",N); //number of particles, default is N (half filling)

    auto nsweeps = basic.getInt("nsweeps");
    auto t1 = basic.getReal("t1",1);
    auto t2 = basic.getReal("t2",0);
    auto U = basic.getReal("U",0);
    auto V1 = basic.getReal("V1",0);
    auto quiet = basic.getYesNo("quiet",false);

    InputGroup table(basic,"sweeps");
    Sweeps sweeps(nsweeps,table);
    println(sweeps);

    //
    // Initialize the site degrees of freedom.
    //
    Hubbard sites(N);

    //
    // Create the Hamiltonian using AutoMPO
    //
    AutoMPO ampo(sites);
    for(int i = 1; i <= N; ++i) 
        {
        ampo += U,"Nupdn",i;
        }
    for(int b = 1; b < N; ++b)
        {
        //Note the + signs for the Hermitian conjugate terms
        ampo += -t1,"Cdagup",b,"Cup",b+1;
        ampo += +t1,"Cup",b,"Cdagup",b+1;
        ampo += -t1,"Cdagdn",b,"Cdn",b+1;
        ampo += +t1,"Cdn",b,"Cdagdn",b+1;

        ampo += V1,"Ntot",b,"Ntot",b+1;
        }
    for(int b = 1; b < N-1; ++b)
        {
        ampo += -t2,"Cdagup",b,"Cup",b+2;
        ampo += +t2,"Cup",b,"Cdagup",b+2;
        ampo += -t2,"Cdagdn",b,"Cdn",b+2;
        ampo += +t2,"Cdn",b,"Cdagdn",b+2;
        }
    auto H = IQMPO(ampo);

    //
    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
    InitState initState(sites);
    int p = Npart;
    for(int i = N; i >= 1; --i) 
        {
        if(p > i)
            {
            println("Doubly occupying site ",i);
            initState.set(i,"UpDn");
            p -= 2;
            }
        else
        if(p > 0)
            {
            println("Singly occupying site ",i);
            initState.set(i,(i%2==1 ? "Up" : "Dn"));
            p -= 1;
            }
        else
            {
            initState.set(i,"Emp");
            }
        }

    IQMPS psi(initState);

    Print(totalQN(psi));

    //
    // Begin the DMRG calculation
    //
    auto energy = dmrg(psi,H,sweeps,{"Quiet",quiet});

    //
    // Measure spin densities
    //
    Vector upd(N),dnd(N);
    for(int j = 1; j <= N; ++j)
        {
        psi.position(j);
        upd(j-1) = (dag(prime(psi.A(j),Site))*sites.op("Nup",j)*psi.A(j)).real();
        dnd(j-1) = (dag(prime(psi.A(j),Site))*sites.op("Ndn",j)*psi.A(j)).real();
        }

    println("Up Density:");
    for(int j = 1; j <= N; ++j)
        printfln("%d %.10f",j,upd(j-1));
    println();

    println("Dn Density:");
    for(int j = 1; j <= N; ++j)
        printfln("%d %.10f",j,dnd(j-1));
    println();

    println("Total Density:");
    for(int j = 1; j <= N; ++j)
        printfln("%d %.10f",j,(upd(j-1)+dnd(j-1)));
    println();

    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",energy);

    return 0;
    }
