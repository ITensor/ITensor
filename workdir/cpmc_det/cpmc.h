#ifndef __CPMC_H_
#define __CPMC_H_

#include "real.h"

namespace itensor {

//
// Perform a constrained path Monte Carlo calculatiion.
// Input
//  Lx: The number of lattice sites in the x direction.
//  Ly: The number of lattice sites in the y direction.
//  N_up: The number of spin-up electrons
//  N_dn: The number of spin-down electrons
//  U: The on-site repulsion strength in the Hubbard Hamiltonian
//  tx: The hopping amplitude between nearest-neighbor sites in the x direction
//  ty: The hopping amplitude between nearest neighbor sites in the y direction
//  deltau: The imaginary time step
//  N_wlk: The number of random walkers
//  N_blksteps: The number of random walk steps in each block
//  N_eqblk: The number of blocks used to equilibrate the random walk before energy 
//           measurement takes place
//  N_blk: The number of blocks used in the measurement phase
//  itv_modsvd: The interval between two adjacent modified Gram-Schmidt 
//              re-orthonormalization of the random walkers.
//  itv_pc: The interval between two adjacent population controls
//  itv_Em: The interval between two adjacent energy measurements
// Output:
//  E_ave: the ground state energy
//  E_err: the standard error in the ground state energy
//

void
CPMC_Lab(Real& E_ave, Real& E_err, int Lx, int Ly, int N_up, int N_dn, 
         Real U, Real tx, Real ty, Real deltau, int N_wlk, int N_blksteps, 
         int N_eqblk, int N_blk, int itv_modsvd, int itv_pc, int itv_Em);

};

#endif
