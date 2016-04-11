#include "itensor/all.h"

using namespace itensor;

int main(int argc, char* argv[])
    {
    int N = 100;

    println("Input J2 value:");
    Real J2 = 0;
    std::cin >> J2;

    //
    // Initialize the site degrees of freedom.
    //
    auto sites = SpinHalf(N);

    //
    // Create the Hamiltonian matrix product operator.
    // Here we use the IQMPO class which is an MPO of 
    // IQTensors, tensors whose indices are sorted
    // with respect to quantum numbers
    //
    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N; ++j)
        {
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        ampo +=     "Sz",j,"Sz",j+1;
        }
    for(int j = 1; j < N-1; ++j)
        {
        ampo += 0.5*J2,"S+",j,"S-",j+2;
        ampo += 0.5*J2,"S-",j,"S+",j+2;
        ampo +=     J2,"Sz",j,"Sz",j+2;
        }
    auto H = IQMPO(ampo);

    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
    auto state = InitState(sites);
    for(int i = 1; i <= N; ++i) 
        state.set(i,(i%2==1 ? "Up" : "Dn"));

    auto psi = IQMPS(state);

    //
    // overlap calculates matrix elements of MPO's with respect to MPS's
    // overlap(psi,H,psi) = <psi|H|psi>
    //
    printfln("Initial energy = %.5f",overlap(psi,H,psi));

    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. Here less than 10 maxm
    // values are provided, so all remaining sweeps will use the
    // last maxm (= 200).
    //
    auto sweeps = Sweeps(5);
    sweeps.maxm() = 50,50,100,100,200;
    sweeps.cutoff() = 1E-8;
    println(sweeps);

    //
    // Begin the DMRG calculation
    //
    auto energy = dmrg(psi,H,sweeps,"Quiet");

    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",energy);
    printfln("\nUsing overlap = %.10f\n", overlap(psi,H,psi) );

    println("\nTotal QN of Ground State = ",totalQN(psi));

    //
    // Measure S.S on every bond
    //
    for(int b = 1; b < N; ++b)
        {
        psi.position(b);
        auto ketzz = psi.A(b)*psi.A(b+1)*sites.op("Sz",b)*sites.op("Sz",b+1);
        auto ketpm = psi.A(b)*psi.A(b+1)*sites.op("Sp",b)*sites.op("Sm",b+1)*0.5;
        auto ketmp = psi.A(b)*psi.A(b+1)*sites.op("Sm",b)*sites.op("Sp",b+1)*0.5;
        auto bra = dag(psi.A(b)*psi.A(b+1));
        bra.prime(Site);
        auto SdS = (bra*ketzz).real() + (bra*ketpm).real() + (bra*ketmp).real();
        printfln("S.S b %d = %.10f",b,SdS);
        }

    return 0;
    }
