#include "eigensolver.h"
#include "localmpo.h"
#include "sweeps.h"
#include "sites/spinone.h"
#include "hams/Heisenberg.h"

using namespace itensor;

int
main(int argc, char* argv[])
    {
    const int N = 100;

    //Model objects represent a collection of 
    //lattice degrees of freedom of a certain type
    SpinOne model(N);

    //Get Hamiltonian
    MPO H = Heisenberg(model);

    //Create MPS
    MPS psi(model); //random starting state
    psi.position(1);

    //Define DMRG sweeps
    Sweeps sweeps(5);
    sweeps.maxm() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;

    //Some stuff we'll need
    //to solve local
    //eigenvalue problem
    LocalMPO<ITensor> Heff(H);

    Real energy = NAN;

    OptSet opts;

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

            //Solve effective eigenvalue problem
            ITensor phi = psi.A(b)*psi.A(b+1);
            energy = davidson(Heff,phi);

            //Update accuracy parameters for svd
            opts.add("Cutoff",sweeps.cutoff(sw));
            opts.add("Maxm",sweeps.maxm(sw));
            opts.add("Minm",sweeps.minm(sw));

            //Define tensor (references/aliases)
            //to hold SVD results
            ITensor& A = psi.Anc(b);   //nc means 'non-const'
            ITensor& B = psi.Anc(b+1); //nc means 'non-const'
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
