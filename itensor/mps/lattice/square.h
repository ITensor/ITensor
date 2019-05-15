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

LatticeGraph inline
squareNextNeighbor(int Nx, 
                   int Ny,
                   Args const& args = Args::global())
    {
    auto yperiodic = args.getBool("YPeriodic",false);
    // Periodicity on y is meaningless for one dimensional chain or a ladder
    yperiodic = yperiodic && (Ny > 2);
    auto N = Nx*Ny;
    LatticeGraph latt; 
    for(int n = 1; n <= N; ++n)
        {
        int x = (n-1)/Ny+1; 
        int y = (n-1)%Ny+1;

        //First-neighbor bonds
        if(x < Nx) 
            {
            //X-direction bond
            latt.emplace_back(n,n+Ny,x,y,x+1,y,"1");
            }
        if(Ny > 1)
            {
            //Y-direction bond
            if(y < Ny) latt.emplace_back(n,n+1,x,y,x,y+1,"1");
            //Periodic Y bond
            if(yperiodic && y == 1) latt.emplace_back(n,n+Ny-1,x,y,x,y+Ny,"1");
            }

        //Second-neighbor bonds
        if(x < Nx && Ny > 1)
            {
            //Next-Neighbor X +Y
            if(y < Ny) latt.emplace_back(n,n+Ny+1,x,y,x+1,y+1,"2");
            //Next-Neighbor X -Y
            if(y > 1) latt.emplace_back(n,n+Ny-1,x,y,x+1,y-1,"2");
            //Periodic Next-Neighbor bonds
            if(yperiodic && y == Ny) 
                {
                //Periodic Next-Neighbor X +Y
                latt.emplace_back(n,n+1,x,Ny,x+1,1,"2");
                }
            if(yperiodic && y == 1) 
                {
                //Periodic Next-Neighbor X -Y
                latt.emplace_back(n,n+2*Ny-1,x,1,x+1,Ny,"2");
                }
            }
        }
    return latt;
    }

} //namespace itensor

#endif
