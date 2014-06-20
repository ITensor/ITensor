#include "core.h"
#include "sites/hubbard.h"
#include "hams/ExtendedHubbard.h"

using namespace std;
using namespace itensor;

int main(int argc, char* argv[])
    {
    //Parse the input file
    if(argc != 2) { printfln("Usage: %s inputfile",argv[0]); return 0; }
    InputGroup basic(argv[1],"basic");

    int N = 0;
    basic.GetIntM("N",N); //the 'M' stands for mandatory
    int Npart = N; //number of particles, default is N (half filling)
    basic.GetInt("Npart",Npart);

    int nsweeps = 0;
    basic.GetIntM("nsweeps",nsweeps);
    Real t1 = 0;
    basic.GetRealM("t1",t1);
    Real t2 = 0;
    basic.GetRealM("t2",t2);
    Real U = 0;
    basic.GetRealM("U",U);
    Real V1 = 0;
    basic.GetRealM("V1",V1);

    InputGroup table(basic,"sweeps");
    Sweeps sweeps(nsweeps,table);
    println(sweeps);

    //
    // Initialize the site degrees of freedom.
    //
    Hubbard sites(N);

    //
    // Create the Hamiltonian matrix product operator.
    // Here we use the IQMPO class which is an MPO of 
    // IQTensors, tensors whose indices are sorted
    // with respect to quantum numbers
    //
    IQMPO H = ExtendedHubbard(sites,
                              Opt("U",U)
                              & Opt("t1",t1)
                              & Opt("t2",t2)
                              & Opt("V1",V1));

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

    println(totalQN(psi));

    //
    // Begin the DMRG calculation
    //
    Real En = dmrg(psi,H,sweeps,"Quiet");

    //
    // Measure spin densities
    //
    Vector upd(N),dnd(N);
    for(int j = 1; j <= N; ++j)
        {
        psi.position(j);
        upd(j) = Dot(conj(primed(psi.A(j),Site)),sites.op("Nup",j)*psi.A(j));
        dnd(j) = Dot(conj(primed(psi.A(j),Site)),sites.op("Ndn",j)*psi.A(j));
        }

    println("Up Density:");
    for(int j = 1; j <= N; ++j)
        printfln("%d %.10f",j,upd(j));
    println();

    println("Dn Density:");
    for(int j = 1; j <= N; ++j)
        printfln("%d %.10f",j,dnd(j));
    println();

    println("Total Density:");
    for(int j = 1; j <= N; ++j)
        printfln("%d %.10f",j,(upd(j)+dnd(j)));
    println();

    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",En);

    return 0;
    }
