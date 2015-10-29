#include "itensor/eigensolver.h"
#include "itensor/mps/localmpo.h"
#include "itensor/mps/sweeps.h"
#include "itensor/mps/sites/spinone.h"
#include "itensor/mps/hams/Heisenberg.h"

using namespace itensor;

int
main(int argc, char* argv[])
    {
    int N = 100;

    //Model objects represent a collection of 
    //lattice degrees of freedom of a certain type
    SpinOne sites(N);

    //Get Hamiltonian MPO
    auto H = MPO(Heisenberg(sites));

    //Create random initial MPS 
    MPS psi(sites);

    //Define DMRG sweeps
    Sweeps sweeps(5);
    sweeps.maxm() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;

    //Some stuff needed to solve
    //projected eigenvalue problem
    LocalMPO<ITensor> Heff(H);

    Real energy = NAN;

    Args args;

    //Loop over sweeps
    for(int sw = 1; sw <= sweeps.nsweep(); ++sw)
        {
        //Loop over bonds
        for(int b = 1, ha = 1; ha != 3; sweepnext(b,ha,N))
            {
            printfln("Sweep=%d, HS=%d, Bond=(%d,%d)",sw,ha,b,b+1);

            //Grow effective H
            psi.position(b);
            Heff.position(b,psi);

            //Compute ITensor phi which is
            //two-site "superblock" wavefunction
            auto phi = psi.A(b)*psi.A(b+1);
            //Solve effective eigenvalue problem
            energy = davidson(Heff,phi);

            //Update accuracy parameters
            //to pass to svd
            args.add("Cutoff",sweeps.cutoff(sw));
            args.add("Maxm",sweeps.maxm(sw));
            args.add("Minm",sweeps.minm(sw));

            //Define aliases (references) A and B 
            //to MPS tensors at sites b and b+1
            auto& A = psi.Anc(b);   //nc means 'non-const'
            auto& B = psi.Anc(b+1); //nc means 'non-const'
            //Singular values will be put on 
            //the diagonal of ITensor D
            ITensor D;

            //Add code:
            //
            // SVD phi and store results
            // into A, D, and B
            //

            if(ha == 1) //sweeping right
                {
                //Add code:
                //
                // Shift orthogonality center
                // to the right
                //
                }
            else
            if(ha == 2) //sweeping left
                {
                //Add code:
                //
                // Shift orthogonality center
                // to the left
                //
                }


            } // for loop over b

        } // for loop over sw

    return 0;
    }
