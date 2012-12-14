#include "core.h"
#include "model/spinhalf.h"
#include "model/spinone.h"
#include "hams/Heisenberg.h"
using boost::format;
using namespace std;

int main(int argc, char* argv[])
    {
    int N = 100;

    //
    // Initialize the site degrees of freedom.
    //
    SpinOne model(N);
    //SpinHalf model(N);

    //
    // Create the Hamiltonian matrix product operator.
    // Here we use the IQMPO class which is an MPO of 
    // IQTensors, tensors whose indices are sorted
    // with respect to quantum numbers
    //
    IQMPO H = Heisenberg(model);

    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
    InitState initState(N);
    for(int i = 1; i <= N; ++i) 
        initState(i) = (i%2==1 ? model.Up(i) : model.Dn(i));

    IQMPS psi(model,initState);

    //
    // psiHphi calculates matrix elements of MPO's with respect to MPS's
    // psiHphi(psi,H,psi) = <psi|H|psi>
    //
    cout << format("Initial energy = %.5f\n")%psiHphi(psi,H,psi);

    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. Here less than 10 maxm
    // values are provided, so all remaining sweeps will use the
    // last maxm (= 200).
    //
    Sweeps sweeps(10);
    sweeps.maxm() = 50,50,100,100,200;
    sweeps.cutoff() = 1E-8;
    cout << sweeps;

    //
    // Begin the DMRG calculation
    //
    Real En = dmrg(psi,H,sweeps,Quiet());

    //
    // Print the final energy reported by DMRG
    //
    cout << format("\nGround State Energy = %.10f\n")%En;
    cout << format("\nUsing psiHphi = %.10f\n") % psiHphi(psi,H,psi);

    cout << "\nTotal QN of Ground State = " << totalQN(psi) << "\n";

    return 0;
    }
