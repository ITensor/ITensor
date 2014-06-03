#include "core.h"
#include "sites/spinhalf.h"
#include "hams/J1J2Chain.h"

using namespace std;
using namespace itensor;

int main(int argc, char* argv[])
    {
    int N = 100;

    cout << "Input J2 value:" << endl;
    Real J2 = 0;
    cin >> J2;

    //
    // Initialize the site degrees of freedom.
    //
    SpinHalf sites(N);

    //
    // Create the Hamiltonian matrix product operator.
    // Here we use the IQMPO class which is an MPO of 
    // IQTensors, tensors whose indices are sorted
    // with respect to quantum numbers
    //
    IQMPO H = J1J2Chain(sites,Opt("J2",J2));

    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
    InitState initState(sites);
    for(int i = 1; i <= N; ++i) 
        initState.set(i,(i%2==1 ? "Up" : "Dn"));

    IQMPS psi(initState);

    Global::opts(Opt("WriteM", 10));
    psi.doWrite(true);
    //
    // psiHphi calculates matrix elements of MPO's with respect to MPS's
    // psiHphi(psi,H,psi) = <psi|H|psi>
    //
    printfln("Initial energy = %.5f",psiHphi(psi,H,psi));

    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. Here less than 10 maxm
    // values are provided, so all remaining sweeps will use the
    // last maxm (= 200).
    //
    Sweeps sweeps(5);
    sweeps.maxm() = 50,50,100,100,200;
    sweeps.cutoff() = 1E-8;
    cout << sweeps;

    //
    // Begin the DMRG calculation
    //
    Real En = dmrg(psi,H,sweeps,"Quiet");

    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",En);
    printfln("\nUsing psiHphi = %.10f\n", psiHphi(psi,H,psi) );

    println("\nTotal QN of Ground State = ",totalQN(psi));

    //
    // Measure S.S on every bond
    //
    psi.doWrite(false);
    for(int b = 1; b < N; ++b)
        {
        psi.position(b);
        IQTensor ketzz = psi.A(b)*psi.A(b+1)*sites.op("Sz",b)*sites.op("Sz",b+1);
        IQTensor ketpm = psi.A(b)*psi.A(b+1)*sites.op("Sp",b)*sites.op("Sm",b+1)*0.5;
        IQTensor ketmp = psi.A(b)*psi.A(b+1)*sites.op("Sm",b)*sites.op("Sp",b+1)*0.5;
        IQTensor bra = conj(psi.A(b)*psi.A(b+1));
        bra.prime(Site);
        Real SdS = Dot(bra,ketzz) + Dot(bra,ketpm) + Dot(bra,ketmp);
        printfln("S.S b %d = %.10f",b,SdS);
        }

    return 0;
    }
