#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

int main(int argc, char* argv[])
    {
    //Parse the input file
    if(argc != 2) { printfln("Usage: %s inputfile_exthubbard",argv[0]); return 0; }
    auto input = InputGroup(argv[1],"input");

    auto N = input.getInt("N");
    auto Npart = input.getInt("Npart",N); //number of particles, default is N (half filling)

    auto nsweeps = input.getInt("nsweeps");
    auto t1 = input.getReal("t1",1);
    auto t2 = input.getReal("t2",0);
    auto U = input.getReal("U",0);
    auto V1 = input.getReal("V1",0);
    auto quiet = input.getYesNo("quiet",false);

    auto table = InputGroup(input,"sweeps");
    auto sweeps = Sweeps(nsweeps,table);
    println(sweeps);

    //
    // Initialize the site degrees of freedom.
    //
    auto sites = Electron(N);

    //
    // Create the Hamiltonian using AutoMPO
    //
    auto ampo = AutoMPO(sites);
    for(int i = 1; i <= N; ++i) 
        {
        ampo += U,"Nupdn",i;
        }
    for(int b = 1; b < N; ++b)
        {
        ampo += -t1,"Cdagup",b,"Cup",b+1;
        ampo += -t1,"Cdagup",b+1,"Cup",b;
        ampo += -t1,"Cdagdn",b,"Cdn",b+1;
        ampo += -t1,"Cdagdn",b+1,"Cdn",b;
        ampo += V1,"Ntot",b,"Ntot",b+1;
        }
    for(int b = 1; b < N-1; ++b)
        {
        ampo += -t2,"Cdagup",b,"Cup",b+2;
        ampo += -t2,"Cdagup",b+2,"Cup",b;
        ampo += -t2,"Cdagdn",b,"Cdn",b+2;
        ampo += -t2,"Cdagdn",b+2,"Cdn",b;
        }
    auto H = toMPO(ampo);

    //
    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
    auto state = InitState(sites);
    int p = Npart;
    for(int i = N; i >= 1; --i) 
        {
        if(p > i)
            {
            println("Doubly occupying site ",i);
            state.set(i,"UpDn");
            p -= 2;
            }
        else
        if(p > 0)
            {
            println("Singly occupying site ",i);
            state.set(i,(i%2==1 ? "Up" : "Dn"));
            p -= 1;
            }
        else
            {
            state.set(i,"Emp");
            }
        }

    auto psi0 = MPS(state);

    Print(totalQN(psi0));

    //
    // Begin the DMRG calculation
    //
    auto [energy,psi] = dmrg(H,psi0,sweeps,{"Quiet",quiet});

    //
    // Measure spin densities
    //
    Vector upd(N),dnd(N);
    for(int j = 1; j <= N; ++j)
        {
        psi.position(j);
        upd(j-1) = elt(dag(prime(psi(j),"Site"))*op(sites,"Nup",j)*psi(j));
        dnd(j-1) = elt(dag(prime(psi(j),"Site"))*op(sites,"Ndn",j)*psi(j));
        }

    println("Up Density:");
    for(int j = 0; j < N; ++j)
        printfln("%d %.10f",1+j,upd(j));
    println();

    println("Dn Density:");
    for(int j = 0; j < N; ++j)
        printfln("%d %.10f",1+j,dnd(j));
    println();

    println("Total Density:");
    for(int j = 0; j < N; ++j)
        printfln("%d %.10f",1+j,(upd(j)+dnd(j)));
    println();

    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",energy);

    return 0;
    }
