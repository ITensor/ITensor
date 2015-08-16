//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_LATTICE_TRIANGULAR_H_
#define __ITENSOR_LATTICE_TRIANGULAR_H_

#include "latticebond.h"

namespace itensor {

Lattice inline
triangularLattice(int Nx, 
                  int Ny,
                  const Args& args = Global::args())
    {
    auto yperiodic = args.getBool("YPeriodic",true);
    auto N = Nx*Ny;
    auto Nbond = 3*N-2*Ny + (yperiodic ? 0 : -2*Nx+1);
    Lattice latt; 
    latt.reserve(Nbond);

    for(int n = 1; n <= N; ++n)
        {
        const int x = (n-1)/Ny+1, 
                  y = (n-1)%Ny+1;

        //X-direction bonds
        if(x < Nx) latt.emplace_back(n,n+Ny);

        if(Ny > 1) //2d bonds
            {
            //vertical bond / Y-periodic diagonal bond
            if((n+1 <= N) && ((y < Ny) || yperiodic)) 
                latt.emplace_back(n,n+1);

            //Periodic vertical bond
            if(yperiodic && y == 1) latt.emplace_back(n,n+Ny-1);

            //Diagonal bonds
            if(x < Nx && y < Ny) latt.emplace_back(n,n+Ny+1);
            }
        }

    if(int(latt.size()) != Nbond) Error("Wrong number of bonds");

    return latt;
    }

} //namespace itensor

#endif
