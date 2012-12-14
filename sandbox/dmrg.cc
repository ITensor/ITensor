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
    SpinOne model(N);    // make a chain of N spin 1's
    //SpinHalf model(N); // make a chain of N spin 1/2's

    //
    // Create the Hamiltonian matrix product operator (MPO)
    //
    MPO H = Heisenberg(model,Opt("J",0.1));

    //
    // Set the initial wavefunction matrix product state (MPS)
    // to be a Neel state.
    //
    InitState initState(N);
    for(int i = 1; i <= N; ++i) 
        initState(i) = (i%2==1 ? model.Up(i) : model.Dn(i));

    MPS psi(model,initState);

    //
    // psiHphi calculates matrix elements of MPO's with respect to MPS's
    // psiHphi(psi,H,psi) = <psi|H|psi>
    //
    cout << format("Initial energy = %.5f") % psiHphi(psi,H,psi) << endl;

    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. Here less than 10 maxm
    // values are provided, so all remaining sweeps will use the
    // last maxm (= 200).
    //
    Sweeps sweeps(10);
    sweeps.maxm() = 50,50,100,100,200;
    sweeps.cutoff() = 1E-10;
    cout << sweeps;

    //
    // Begin the DMRG calculation
    //
    Real En = dmrg(psi,H,sweeps,Quiet());

    //
    // Print the final energy reported by DMRG
    //
    cout << format("\nGround State Energy = %.10f")%En << endl;

    return 0;
    }
