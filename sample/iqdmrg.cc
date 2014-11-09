#include "core.h"
#include "sites/spinhalf.h"
#include "sites/spinone.h"
#include "autompo.h"

using namespace itensor;

int 
main(int argc, char* argv[])
    {
    int N = 100;

    //
    // Initialize the site degrees of freedom.
    //
    //SpinHalf sites(N); //make a chain of N spin 1/2's
    SpinOne sites(N); //make a chain of N spin 1's

    //
    // Use the AutoMPO feature to create the 
    // next-neighbor Heisenberg model.
    //
    // Here we convert the AutoMPO information
    // into an IQMPO, a matrix-product operator
    // which automatically tracks quantum
    // number information.
    //
    AutoMPO ampo(sites);
    for(int j = 1; j < N; ++j)
        {
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        ampo +=     "Sz",j,"Sz",j+1;
        }
    auto H = IQMPO(ampo);

    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
    InitState initState(sites);
    for(int i = 1; i <= N; ++i) 
        {
        if(i%2 == 1)
            initState.set(i,"Up");
        else
            initState.set(i,"Dn");
        }

    IQMPS psi(initState);

    //
    // psiHphi calculates matrix elements of MPO's with respect to MPS's
    // psiHphi(psi,H,psi) = <psi|H|psi>
    //
    printfln("Initial energy = %.5f", psiHphi(psi,H,psi) );

    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. 
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).
    //
    Sweeps sweeps(5);
    sweeps.maxm() = 10,20,100,100,200;
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
    printfln("\nUsing psiHphi = %.10f", psiHphi(psi,H,psi) );

    println("\nTotal QN of Ground State = ",totalQN(psi));

    return 0;
    }
