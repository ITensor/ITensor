//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_LATTICE_SQUARE_H_
#define __ITENSOR_LATTICE_SQUARE_H_

#include "latticebond.h"


namespace itensor {

Lattice inline
squareLattice(int Nx, 
              int Ny,
              const Args& args = Global::args())
    {
    auto yperiodic = args.getBool("YPeriodic",false);
    auto N = Nx*Ny;
    auto Nbond = 2*N-Ny + (yperiodic ? 0 : -Nx);
    Lattice latt; 
    latt.reserve(Nbond);
    for(int n = 1; n <= N; ++n)
        {
        int x = (n-1)/Ny+1, 
            y = (n-1)%Ny+1;

        //X-direction bond
        if(x < Nx) latt.emplace_back(n,n+Ny);

        if(Ny > 1)
            {
            //Y-direction bond
            if(y < Ny) latt.emplace_back(n,n+1);
            //Periodic bond
            if(yperiodic && y == 1) latt.emplace_back(n,n+Ny-1);
            }
        }
    if(int(latt.size()) != Nbond) Error("Square latt wrong number of bonds");
    return latt;
    }

} //namespace itensor

#endif
