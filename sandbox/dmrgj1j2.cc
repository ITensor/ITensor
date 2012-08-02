#include "core.h"
#include "model/spinhalf.h"
#include "hams/J1J2Chain.h"
using boost::format;
using namespace std;

int main(int argc, char* argv[])
    {
    int N = 100;

    cout << "Input J2 value:" << endl;
    Real J2 = 0;
    cin >> J2;

    //
    // Initialize the site degrees of freedom.
    //
    SpinHalf model(N);

    //
    // Create the Hamiltonian matrix product operator.
    // Here we use the IQMPO class which is an MPO of 
    // IQTensors, tensors whose indices are sorted
    // with respect to quantum numbers
    //
    IQMPO H = J1J2Chain(model,J2,1);

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
    Sweeps sweeps(5);
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

    //
    // Measure S.S on every bond
    //
    for(int b = 1; b < N; ++b)
        {
        psi.position(b);
        IQTensor ketzz = psi.AA(b)*psi.AA(b+1)*model.sz(b)*model.sz(b+1);
        IQTensor ketpm = psi.AA(b)*psi.AA(b+1)*model.sp(b)*model.sm(b+1)*0.5;
        IQTensor ketmp = psi.AA(b)*psi.AA(b+1)*model.sm(b)*model.sp(b+1)*0.5;
        IQTensor bra = conj(psi.AA(b)*psi.AA(b+1));
        bra.primesite();
        Real SdS = Dot(bra,ketzz) + Dot(bra,ketpm) + Dot(bra,ketmp);
        cout << format("S.S b %d = %.10f") % b % SdS << endl;
        }

    return 0;
    }
