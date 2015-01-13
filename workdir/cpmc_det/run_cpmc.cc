#include "core.h"
#include "cpmc.h"

using std::cout;
using std::endl;
using std::string;
using namespace itensor;

//
//  A script to set the input parameters and run a CPMC calculation
//
// Huy Nguyen, Hao Shi, Jie Xu and Shiwei Zhang
// ©2014 v1.0
// Package homepage: http://cpmc-lab.wm.edu
// Distributed under the 
// <a href="matlab: web('http://cpc.cs.qub.ac.uk/licence/licence.html')">
// Computer Physics Communications Non-Profit Use License</a>
// Any publications resulting from either applying or building on the present package 
// should cite the following journal article (in addition to the relevant literature 
// on the method):
// "CPMC-Lab: A Matlab Package for Constrained Path Monte Carlo Calculations" 
// Comput. Phys. Commun. (2014)
//

int 
main(int argc, char* argv[])
    {
    //
    // system parameters:
    //
    int Lx = 4;     // The number of lattice sites in the x direction
    int Ly = 1;     // The number of lattice sites in the y direction

    int N_up = 2;   // The number of spin-up electrons
    int N_dn = 2;   // The number of spin-down electrons

    Real U = 1.0;   // The on-site repulsion strength in the Hubbard Hamiltonian
    Real tx = 1.0;  // The hopping amplitude between nearest-neighbor sites in the 
                    // x direction
    Real ty = 1.0;  // The hopping amplitude between nearest neighbor sites in the 
                    // y direction

    //
    // run parameters:
    //
    Real deltau = 0.01;     // The imaginary time step
    int N_wlk = 10;          // The number of random walkers
    int N_blksteps = 10;     // The number of random walk steps in each block
    int N_eqblk = 10;        // The number of blocks used to equilibrate the random 
                            // walk before energy measurement takes place
    int N_blk = 10;          // The number of blocks used in the measurement phase
    int itv_modsvd = 1;     // The interval between two adjacent modified Gram-Schmidt 
                            // re-orthonormalization of the random walkers. No 
                            // re-orthonormalization if itv_modsvd > N_blksteps
    int itv_pc = 1;         // The interval between two adjacent population controls. 
                            // No population control if itv_pc > N_blksteps
    int itv_Em = 1;         // The interval between two adjacent energy measurements
  
    //
    // initialize output values:
    //
    Real E_ave = 0.0,
         E_err = 0.0;
    
    //
    // invoke the main function
    //
    CPMC_Lab(E_ave, E_err, Lx, Ly, N_up, N_dn, U, tx, ty, deltau, N_wlk, N_blksteps,
             N_eqblk, N_blk, itv_modsvd, itv_pc, itv_Em);
    
    //
    // output
    //
    cout << "E_ave = " << E_ave << endl;
    cout << "E_err = " << E_err << endl;
    cout << endl;
    
    return 0;
    }


