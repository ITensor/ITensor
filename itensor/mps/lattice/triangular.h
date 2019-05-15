//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef __ITENSOR_LATTICE_TRIANGULAR_H_
#define __ITENSOR_LATTICE_TRIANGULAR_H_

#include "itensor/mps/lattice/latticebond.h"

namespace itensor {

LatticeGraph inline
triangularLattice(int Nx, 
                  int Ny,
                  Args const& args = Args::global())
    {
    auto yperiodic = args.getBool("YPeriodic",true);
    // Periodicity on y is meaningless for one dimensional chain or a ladder
    yperiodic = yperiodic && (Ny > 2);
    auto N = Nx*Ny;
    auto Nbond = 3*N-2*Ny + (yperiodic ? 0 : -2*Nx+1);
    LatticeGraph latt; 
    latt.reserve(Nbond);

    for(int n = 1; n <= N; ++n)
        {
        int x = (n-1)/Ny+1; 
        int y = (n-1)%Ny+1;

        //X-direction bonds
        if(x < Nx) latt.emplace_back(n,n+Ny);

        if(Ny > 1) //2d bonds
            {
            //vertical bond / Y-periodic diagonal bond
            if((n+1 <= N) && ((y < Ny) || yperiodic)) 
                {
                latt.emplace_back(n,n+1);
                }

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
