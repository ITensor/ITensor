//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_LATTICE_SQUARE_H_
#define __ITENSOR_LATTICE_SQUARE_H_

#include "itensor/mps/lattice/latticebond.h"

namespace itensor {

LatticeGraph inline
squareLattice(int Nx, 
              int Ny,
              Args const& args = Args::global())
    {
    auto yperiodic = args.getBool("YPeriodic",false);
    // Periodicity on y is meaningless for one dimensional chain or a ladder
    yperiodic = yperiodic && (Ny > 2);
    auto N = Nx*Ny;
    auto Nbond = 2*N-Ny + (yperiodic ? 0 : -Nx);
    LatticeGraph latt; 
    latt.reserve(Nbond);
    for(int n = 1; n <= N; ++n)
        {
        int x = (n-1)/Ny+1; 
        int y = (n-1)%Ny+1;

        //X-direction bond
        if(x < Nx) latt.emplace_back(n,n+Ny,x,y,x+1,y);

        if(Ny > 1)
            {
            //Y-direction bond
            if(y < Ny) latt.emplace_back(n,n+1,x,y,x,y+1);
            //Periodic bond
            if(yperiodic && y == 1) latt.emplace_back(n,n+Ny-1,x,y,x,y+Ny);
            }
        }
    if(int(latt.size()) != Nbond) Error("Square latt wrong number of bonds");
    return latt;
    }

} //namespace itensor

#endif
