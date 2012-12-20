#include "eigensolver.h"
#include "localmpo.h"
#include "Sweeps.h"
#include "model/spinone.h"
#include "hams/Heisenberg.h"

using namespace std;
using boost::format;

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
    Eigensolver solver;
    solver.debugLevel(1);

    Real energy = NAN;

    //Loop over sweeps
    for(int sw = 1; sw <= sweeps.nsweep(); ++sw)
        {
        //Loop over bonds
        for(int b = 1, ha = 1; ha != 3; sweepnext(b,ha,N))
            {
            cout << format("Sweep=%d, HS=%d, Bond=(%d,%d)")
                    % sw
                    % ha
                    % b
                    % (b+1)
                    << endl;

            //Grow effective H
            psi.position(b);
            Heff.position(b,psi);

            //Solve effective eigenvalue problem
            ITensor phi = psi.A(b)*psi.A(b+1);
            energy = solver.davidson(Heff,phi);

            //Construct SVDWorker and set
            //accuracy parameters
            SVDWorker W;
            W.cutoff(sweeps.cutoff(sw)); 
            W.minm(sweeps.minm(sw)); 
            W.maxm(sweeps.maxm(sw));

            //Define tensor (references)
            //to hold SVD results
            ITensor& A = psi.Anc(b);   //nc means 'non-const'
            ITensor& B = psi.Anc(b+1); //nc means 'non-const'
            ITSparse D;

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
