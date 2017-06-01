#include "itensor/all.h"

using namespace itensor;

int main()
    {
    int N = 100;

    //sites objects represent the Hilbert space, 
    //a collection of "physical" indices
    auto sites = SpinOne(N);

    //Use AutoMPO to make Hamiltonian MPO
    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N; ++j)
        {
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        ampo +=     "Sz",j,"Sz",j+1;
        }
    auto H = MPO(ampo);

    //Create MPS
    auto psi = MPS(sites); //random starting state

    //Define DMRG sweeps
    auto sweeps = Sweeps(5);
    sweeps.maxm() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;

    //Some stuff needed to solve
    //projected eigenvalue problem
    auto Heff = LocalMPO<ITensor>(H);

    Real energy = NAN;

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

            //Update accuracy parameters
            //to pass to svd
            auto args = Args("Cutoff",sweeps.cutoff(sw),
                             "Maxm",sweeps.maxm(sw),
                             "Minm",sweeps.minm(sw));

            //Define tensor (references/aliases)
            //to hold SVD results
            auto& A = psi.Aref(b);  //ref means reference
            auto& B = psi.Aref(b+1); //ref means reference
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
                // to the right by multiplying
                // singular values into B
                //
                }
            else
            if(ha == 2) //sweeping left
                {
                //Add code:
                //
                // Shift orthogonality center
                // to the right by multiplying
                // singular values into A
                //
                }

            } // for loop over b

        } // for loop over sw

    return 0;
    }
