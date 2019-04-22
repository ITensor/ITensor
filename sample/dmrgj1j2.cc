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
    auto sites = SpinHalf(N,{"ConserveQNs=",true});

    //
    // Create the Hamiltonian matrix product operator.
    // Here we use the MPO class which is an MPO of 
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
    auto H = toMPO(ampo);

    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
    auto state = InitState(sites);
    for(int i = 1; i <= N; ++i) 
        state.set(i,(i%2==1 ? "Up" : "Dn"));

    auto psi0 = MPS(state);

    //
    // inner calculates matrix elements of MPO's with respect to MPS's
    // inner(psi0,H,psi0) = <psi0|H|psi0>
    //
    printfln("Initial energy = %.5f",inner(psi0,H,psi0));

    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. Here less than 10 maxdim
    // values are provided, so all remaining sweeps will use the
    // last maxdim (= 200).
    //
    auto sweeps = Sweeps(5);
    sweeps.maxdim() = 50,50,100,100,200;
    sweeps.cutoff() = 1E-8;
    println(sweeps);

    //
    // Begin the DMRG calculation
    //
    auto [energy,psi] = dmrg(H,psi0,sweeps,"Quiet");

    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",energy);
    printfln("\nUsing inner = %.10f\n", inner(psi,H,psi) );

    println("\nTotal QN of Ground State = ",totalQN(psi));

    //
    // Measure S.S on every bond
    //
    for(int b = 1; b < N; ++b)
        {
        psi.position(b);
        auto ketzz = psi(b)*psi(b+1)*op(sites,"Sz",b)*op(sites,"Sz",b+1);
        auto ketpm = psi(b)*psi(b+1)*op(sites,"Sp",b)*op(sites,"Sm",b+1)*0.5;
        auto ketmp = psi(b)*psi(b+1)*op(sites,"Sm",b)*op(sites,"Sp",b+1)*0.5;
        auto bra = dag(psi(b)*psi(b+1));
        bra.prime("Site");
        auto SdS = elt(bra*ketzz) + elt(bra*ketpm) + elt(bra*ketmp);
        printfln("S.S b %d = %.10f",b,SdS);
        }

    return 0;
    }
