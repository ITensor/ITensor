#include "cpmc.h"
#include "matrix.h"
#include "math.h"

using namespace std;

#define rNum()    ( (double)rand() / (1.0+(double)RAND_MAX) )

namespace itensor {

//
// Propagate a walker by the kinetic energy propagator exp(-deltau*K/2)
// Inputs:
//  phi: the matrix of a single walker
//  w: the weight of that walker
//  O: the overlap of that walker
//  Proj_k_half: the matrix of the operator exp(-deltau*K/2)
//  Phi_T: the matrix of the trial wave function
//  N_up: the number of spin up electrons
//  N_par: the total number of electrons
// Outputs:
//  phi: the matrix of the propagated single walker
//  w: the weight of the propagated walker
//  O: the overlap of the propagated walker   
//  invO_matrix_up: the inverse of the spin up sector of the walker's overlap matrix 
//  invO_matrix_dn: the inverse of the spin down sector of the walker's overlap matrix
// 

void
halfK(Matrix& phi, Real& w, Real& O, Matrix& invO_matrix_up, 
      Matrix& invO_matrix_dn, Matrix Proj_k_half, Matrix Phi_T, 
      int N_up, int N_par)
    {
    //
    // propagate the walker by exp(-deltau*K/2)
    //
    phi = Proj_k_half*phi;
    
    int N_site = phi.Nrows();
   
    //
    // update the inverse of the overlap
    //
    Real detinvO_matrix_up = 1.0;
    Real detinvO_matrix_dn = 1.0;
    if(N_up > 0)
        {
        Matrix Phi_T_up = Phi_T.SubMatrix(1,N_site,1,N_up);
        Matrix phi_up = phi.SubMatrix(1,N_site,1,N_up);
        invO_matrix_up = Inverse(Phi_T_up.t()*phi_up);
        detinvO_matrix_up = Determinant(invO_matrix_up);
        }
    
    if(N_par-N_up > 0)
        {
        Matrix Phi_T_dn = Phi_T.SubMatrix(1,N_site,N_up+1,N_par);
        Matrix phi_dn = phi.SubMatrix(1,N_site,N_up+1,N_par);
        invO_matrix_dn = Inverse(Phi_T_dn.t()*phi_dn);
        detinvO_matrix_dn = Determinant(invO_matrix_dn);
        }
    // calculate the new overlap 
    Real O_new = 1.0/(detinvO_matrix_up*detinvO_matrix_dn);
    Real O_ratio = O_new/O;

    //
    // enforce the constrained path condition
    //
    // If the new weight is negative (O_raio<0), kill the walker by setting 
    // its weight to zero
    // real(O_ratio) enforces the phase-free approximation in case of complex 
    // phases (because the condition O_ratio>0 only checks the real part of O_ratio)
    if(O_ratio > 0)
        {
        O=O_new;
        w=w*O_ratio;
        }
    else
        {
        w = 0;
        }

    } // End halfK()

//
// Sample the auxiliary field over a single lattice site for a single walker and 
// propagate that walker by the potential energy propagator exp(-deltau*V)
// Inputs:
//  phi: a single row in the matrix of a walker (corresponding to the amplitude 
//  of all electrons over a single lattice site in that walker)
//  phi_T: the matrix of the trial wave function
//  N_up: the number of spin up electrons
//  N_par: the total number of electrons
//  O: the overlap of the aforementioned walker
//  w: the weight of the aforementioned walker
//  invO_matrix_up: the inverse of the spin up sector of the walker's overlap matrix 
//  invO_matrix_dn: the inverse of the spin down sector of the walker's overlap matrix
//  aux_fld: the 2x2 matrix containing all the possible values of the quantity 
//  exp(gamma*s(sigma)*x_i) (used in V.m only)
// Outputs:
//  phi: the propagated row of the aforementioned walker
//  O: the overlap after propagation of the aforementioned walker
//  w: the weight after propagation of the aforementioned walker
//  invO_matrix_up: the updated inverse of the spin up sector of the walker's 
//  overlap matrix 
//  invO_matrix_dn: the updated inverse of the spin down sector of the walker's 
//  overlap matrix 
//

void
V(Vector& phi, Vector phi_T, int N_up, int N_par, Real& O, 
  Real& w, Matrix& invO_matrix_up, Matrix& invO_matrix_dn, 
  Matrix aux_fld)
    {
    //
    // Pre-allocate matrices:
    //
    Vector Gii(2);
    Matrix RR(2,2);
    Gii = 0.0;
    RR = 0.0;

    //
    // Calculate the Green's function
    //
    Vector temp1_up(N_up);
    Vector temp2_up(N_up);
    Vector temp1_dn(N_par-N_up);
    Vector temp2_dn(N_par-N_up);
    temp1_up = 0.0;
    temp2_up = 0.0;
    temp1_dn = 0.0;
    temp2_dn = 0.0;
    if(N_up > 0)
        {
        temp1_up = phi.SubVector(1,N_up)*invO_matrix_up;
        temp2_up = invO_matrix_up*phi_T.SubVector(1,N_up);
        Gii(1) = temp1_up*phi_T.SubVector(1,N_up);
        }

    if(N_par-N_up > 0)
        {
        temp1_dn = phi.SubVector(N_up+1,N_par)*invO_matrix_dn;
        temp2_dn = invO_matrix_dn*phi_T.SubVector(N_up+1,N_par);
        Gii(2) = temp1_dn*phi_T.SubVector(N_up+1,N_par);
        }
    RR(1,1) = (aux_fld(1,1)-1.0)*Gii(1) + 1.0;
    RR(1,2) = (aux_fld(1,2)-1.0)*Gii(1) + 1.0;
    RR(2,1) = (aux_fld(2,1)-1.0)*Gii(2) + 1.0;
    RR(2,2) = (aux_fld(2,2)-1.0)*Gii(2) + 1.0;

    //
    // Perform the importance sampling and propagate the walker
    //
    // compute overlaps
    Vector O_ratio_temp(2);
    O_ratio_temp(1) = RR(1,1)*RR(2,1);
    O_ratio_temp(2) = RR(1,2)*RR(2,2);
    Vector O_ratio_temp_real(2);
    O_ratio_temp_real = 0.0;
    if(O_ratio_temp(1) > 0.0)
        O_ratio_temp_real(1) = O_ratio_temp(1);
    if(O_ratio_temp(2) > 0.0)
        O_ratio_temp_real(2) = O_ratio_temp(2);
    // the normalization for the importance-sampled pdf
    Real sum_O_ratio_temp_real = O_ratio_temp_real(1)+O_ratio_temp_real(2);
    // if both auxiliary fields lead to negative overlap then kill the walker
    if(sum_O_ratio_temp_real <= 0)
        w=0;
    if(w > 0)
        {
        // Otherwise update the weight
        w=w*0.5*sum_O_ratio_temp_real;
        
        int x_spin;
        if(O_ratio_temp_real(1)/sum_O_ratio_temp_real >= rNum())
            x_spin=1;
        else
            x_spin=2;
        // propagates the walker with the chosen auxiliary field
        if(N_up > 0)
            phi.SubVector(1,N_up) = phi.SubVector(1,N_up)*aux_fld(1,x_spin);
        if(N_par-N_up > 0)
            phi.SubVector(N_up+1,N_par) = phi.SubVector(N_up+1,N_par)*
                                                    aux_fld(2,x_spin);
        
        // Update the overlap using Sherman-Morrison
        O = O*O_ratio_temp(x_spin);
        invO_matrix_up = invO_matrix_up + (1-aux_fld(1,x_spin))/
                                 RR(1,x_spin)*temp2_up*temp1_up;
        invO_matrix_dn = invO_matrix_dn + (1-aux_fld(2,x_spin))/
                                 RR(2,x_spin)*temp2_dn*temp1_dn;
        }
    } // End V()

//
// Calculate the mixed estimator for the ground state energy of a walker
// Inputs:
//  H_k: the one-body kinetic Hamiltonian
//  phi: the matrix of a single walker
//  Phi_T: the matrix of the trial wave function
//  invO_matrix_up: the inverse of the spin up sector of the walker's overlap matrix 
//  invO_matrix_dn: the inverse of the spin down sector of the walker's overlap matrix
//  N_up: the number of spin up electrons
//  N_par: the total number of electrons of both spins
//  U: the on-site repulsion strength in the Hubbard model
// Outputs:
//  e: the mixed estimator for the ground state energy of the input walker
//

Real
measure(Matrix H_k, Matrix phi, Matrix Phi_T, Matrix invO_matrix_up, 
        Matrix invO_matrix_dn, int N_up, int N_par, Real U)
    {
    Real e = 0.0;
    int N_sites = phi.Nrows();

    //
    // calculate the single-particle Green's function matrix for each spin:
    //
    Vector diag_G_up(N_sites);
    diag_G_up = 0.0;
    Vector diag_G_dn(N_sites);
    diag_G_dn = 0.0;
    Matrix G_up(N_sites,N_sites);
    Matrix G_dn(N_sites,N_sites);
    G_up = 0.0;
    G_dn = 0.0;
    
    //
    // calculate the potential energy:
    //
    if(N_up > 0)
        {
        Matrix temp_up = phi.SubMatrix(1,N_sites,1,N_up)*invO_matrix_up;
        G_up = temp_up*(Phi_T.SubMatrix(1,N_sites,1,N_up)).t();
        diag_G_up = G_up.Diagonal();
        }
    if(N_par - N_up > 0)
        {
        Matrix temp_dn = phi.SubMatrix(1,N_sites,N_up+1,N_par)*invO_matrix_dn;
        G_dn = temp_dn*(Phi_T.SubMatrix(1,N_sites,N_up+1,N_par)).t();
        diag_G_dn = G_dn.Diagonal();
        }
    Real n_int = diag_G_up*diag_G_dn;
    Real potentialEnergy = n_int*U;
   
    //
    // calculate the kinetic energy:
    //
    Real kineticEnergy = 0.0;
    for(int i = 1; i <= N_sites; i++)
        {
        for(int j = 1; j <= N_sites; j++)
            kineticEnergy += H_k(i,j)*(G_up(i,j)+G_dn(i,j));
        }
   
    //
    // calculate the total energy:
    //
    e = potentialEnergy + kineticEnergy;

    return e;
    } // End measure()

//
// Perform one step of the random walk
// Inputs:
//  phi: the whole ensemble of walkers
//  N_wlk: the number of walkers
//  N_sites: the total number of lattice sites
//  w: the array of weights of all the walkers
//  O: the array of overlaps of all the walkers
//  E: the total energy of all walkers
//  W: the total weight of all walkers
//  H_k: the one-body kinetic Hamiltonian
//  Proj_k_half: the matrix of the operator exp(-deltau*K/2)
//  flag_mea: the flag (1 or 0) that specifies whether the energy should the 
//  measured in this step
//  Phi_T: the matrix of the trial wave function
//  N_up: the number of spin up electrons
//  N_par: the total number of electrons
//  U: the on-site repulsion strength in the Hubbard model
//  fac_norm: the exponent of the pre-factor exp(-deltau*(H-E_T))
//  aux_fld: the 2x2 matrix containing all the possible values of the quantity 
//  exp(gamma*s(sigma)*x_i) (used in V.m only)
// Outputs:
//  phi: the ensemble of walkers after propagation
//  w: the new array of weights of all walkers
//  O: the new array of overlaps of all walkers
//  E: the new total energy of all walkers
//  W: the new total weight of all walkers
//

void
stepwlk(std::vector<Matrix>& phi, int N_wlk, int N_sites, Vector& w, 
        Vector& O, Real& E, Real& W, Matrix H_k, Matrix Proj_k_half, 
        int flag_mea, Matrix Phi_T, int N_up, int N_par, Real U, 
        Real fac_norm, Matrix aux_fld)
    {
    //
    // Propagate each walker:
    //
    Vector e(N_wlk); // Array containing the energy
    e = 0.0;
    Matrix invO_matrix_up;
    Matrix invO_matrix_dn;

    for(int i_wlk = 1; i_wlk <= N_wlk; i_wlk++)
        {
        Matrix Phi = phi[i_wlk];
        if(w(i_wlk) > 0)
            {
            // multiply by the pre-factor exp(-deltau*(E_T)) in the ground-state 
            // projector and by the prefactor exp(-0.5*U*(N_up+N_dn)) in the 
            // Hirsch transformation
            w(i_wlk) = w(i_wlk)*exp(fac_norm);
            // propagate by the kinetic term exp(-1/2*deltau*K)
            halfK(Phi, w(i_wlk), O(i_wlk), invO_matrix_up, invO_matrix_dn, 
                  Proj_k_half, Phi_T, N_up, N_par);
            if(w(i_wlk) > 0)
                {
                // propagate each lattice site of a walker by the potential term:
                for(int j_site = 1; j_site <= N_sites; j_site++)
                    {
                    if(w(i_wlk) > 0)
                        {
                        Vector Phi_site = Phi.Row(j_site);
                        V(Phi_site, Phi_T.Row(j_site), N_up, N_par, O(i_wlk), 
                          w(i_wlk), invO_matrix_up, invO_matrix_dn, aux_fld);
                        Phi.Row(j_site) = Phi_site;
                        }
                    }
                }
            if(w(i_wlk) > 0)
                {
                // propagate by the kinetic term exp(-1/2*deltau*K)
                halfK(Phi, w(i_wlk), O(i_wlk), invO_matrix_up, invO_matrix_dn, 
                      Proj_k_half, Phi_T, N_up, N_par);
                if(w(i_wlk) > 0)
                    {
                    // measure the energy if needed:
                    if(flag_mea == 1)
                        e(i_wlk) = measure(H_k, Phi, Phi_T,  invO_matrix_up, 
                                           invO_matrix_dn, N_up, N_par, U); 
                    }
                }
            }
        phi[i_wlk] = Phi;
        }
    
    //
    // Compute the ensemble's total energy and weight if measurement took place
    //
    if(flag_mea == 1)
        {
        for(int i_wlk = 1; i_wlk <= N_wlk; i_wlk++)
            {
            if(w(i_wlk) > 0)
                {
                E = E + e(i_wlk)*w(i_wlk);
                W = W + w(i_wlk);
                }
            }
        }
    } // End stepwlk()

//
// Perform the modified Gram-Schmidt orthogonalization to stabilize the walkers
// Inputs:
//  Phi: the whole ensemble of walkers
//  N_wlk: the number of walkkers
//  O: the array of overlaps of all walkers
//  N_up: the number of spin up electrons
//  N_par: the total number of electrons
// Outputs:
//  Phi: the stabilized ensemble of walkers
//  O: the updated array of overlaps
//

void
stblz(std::vector<Matrix>& Phi, int N_wlk, Vector& O, int N_up, int N_par)
    {
    int N_sites = Phi[1].Nrows();
    int N_dn = N_par - N_up;
    
    //
    // Perform the QR decomposition on each walker
    //
    // Keep only the Q matrices and discard the R matrices
    for(int i_wlk = 1; i_wlk <= N_wlk; i_wlk++)
        {
        Matrix phi = Phi[i_wlk];
        Real det_R_up = 1.0;
        Real det_R_dn = 1.0;
        // for the spin up sector: 
        if(N_up > 0)
            {
            Matrix Q_up;
            Matrix R_up;
            Matrix phi_up(N_sites,N_sites);
            phi_up = 0.0;
            phi_up.SubMatrix(1,N_sites,1,N_up) = phi.SubMatrix(1,N_sites,1,N_up);
            QRDecomp(phi_up,Q_up,R_up);
            phi.SubMatrix(1,N_sites,1,N_up) = Q_up.SubMatrix(1,N_sites,1,N_up);
            det_R_up = Determinant(R_up.SubMatrix(1,N_up,1,N_up));
            }
        // for the spin down sector:
        if(N_dn > 0)
            {
            Matrix Q_dn;
            Matrix R_dn;
            Matrix phi_dn(N_sites,N_sites);
            phi_dn = 0.0;
            phi_dn.SubMatrix(1,N_sites,1,N_dn) = phi.SubMatrix(1,N_sites,N_up+1,N_par);
            QRDecomp(phi_dn,Q_dn,R_dn);
            phi.SubMatrix(1,N_sites,N_up+1,N_par) = Q_dn.SubMatrix(1,N_sites,1,N_dn);
            det_R_dn = Determinant(R_dn.SubMatrix(1,N_dn,1,N_dn));
            }
        Phi[i_wlk] = phi;
        // Update the weight of each walker
        O(i_wlk) = O(i_wlk)/det_R_up/det_R_dn;
        }
    } // End stblz()

//
// Perform population control with a simple "combing" method
// Inputs:
//  Phi: the whole ensemble of walkers
//  w: array containing the weights of all the walkers
//  O: array containing the overlaps of all the walkers
//  N_wlk: the number of walkers
//  N_sites: the total number of lattice sites
//  N_par: the total number of electrons
// Outputs:
//  Phi: the new ensemble of walkers after population control
//  w: the new array of weights
//  O: the new array of overlaps
// 

void
pop_cntrl(std::vector<Matrix>& Phi, Vector& w, Vector& O, int N_wlk, int N_sites, int N_par)
    {
    //
    // Preparation
    //
    // Create empty matrices that will hold the outputs
    // in the end the number of walkers will still be N_wlk
    std::vector<Matrix> new_Phi(N_wlk+1);
    Matrix zeros(N_sites, N_par);
    zeros = 0.0;
    for(int i = 1; i <= N_wlk; i++)
        new_Phi[i] = zeros;
    Vector new_O(N_wlk);
    new_O = 0.0;
    
    Real sum_w = 0.0;
    for(int i = 1; i <= N_wlk; i++)
        sum_w += w(i);
    // scaling factor to bring the current total weight back to the original level (=N_wlk)
    Real d = 1.0*N_wlk/sum_w;
    // start the "comb" at a random position to avoid bias against the first walker 
    sum_w = -rNum();
    int n_wlk = 0;

    //
    // Apply the comb
    //
    for(int i_wlk = 1; i_wlk <= N_wlk; i_wlk++)
        {
        sum_w += w(i_wlk)*d;
        int n = ceil(sum_w);
        for(int j = n_wlk+1; j <= n; j++)
            {
            new_Phi[j] = Phi[i_wlk];
            new_O(j) = O(i_wlk);
            }
        n_wlk = n;
        }

    //
    // Return the new population, weights and overlaps:
    //
    Phi = new_Phi;
    O = new_O;
    // All new walkers have weights to 1 and the total weight = N_wlk
    w = Vector(N_wlk,1.0);
    } // End pop_cntrl()

//
// Generate the one-body kinetic term of the Hubbard Hamiltonian with 
// the given parameters
// Input:
//  Lx: The number of lattice sites in the x direction.
//  Ly: The number of lattice sites in the y direction.
//  tx: The hopping amplitude between nearest-neighbor sites in the x direction
//  ty: The hopping amplitude between nearest neighbor sites in the y direction
// Output
//  H_k: The one-body kinetic Hamiltonian in the form of a square matrix of 
//  size (Lx*Ly) 
//

Matrix
H_K(int Lx, int Ly, Real tx, Real ty)
    {
    int N_sites = Lx*Ly;
    Matrix H_k(N_sites,N_sites);
    H_k = 0.0;

    int r = 0;
    for(int iy = 1; iy <= Ly; iy++)
        {
        for(int jx = 1; jx <= Lx; jx++)
            {
            r++;            // r=(iy-1)*Lx+jx;
            if(Lx != 1)
                {
                if(jx == 1)
                    {
                    H_k(r,r+1) = H_k(r,r+1) - tx;
                    }
                else if(jx == Lx)
                    {
                    H_k(r,r-1) = H_k(r,r-1) - tx;
                    }
                else
                    {
                    H_k(r,r-1) = -tx;
                    H_k(r,r+1) = -tx;
                    }
                }

            if(Ly != 1)
                {
                if(iy == 1)
                    {
                    H_k(r,r+Lx) = H_k(r,r-Lx) - ty;
                    }
                else if(iy == Ly)
                    {
                    H_k(r,r-Lx) = H_k(r,r-Lx) - ty;
                    }
                else
                    {
                    H_k(r,r-Lx) = -ty;
                    H_k(r,r-Lx) = -ty;
                    }
                }
            }
        }
   
    return H_k;
    } // End H_K()

void
initialization(int Lx, int Ly, Real tx, Real ty, int N_up, int N_dn,
               Real deltau, Real U, int N_wlk, int N_blk,
               int& N_sites, int& N_par, Matrix& H_k, Matrix& Proj_k_half, 
               Matrix& Phi_T, std::vector<Matrix>& Phi, Vector& w, Vector& O, 
               Vector& E_blk, Vector& W_blk, Real& fac_norm, Real& gamma, 
               Matrix& aux_fld)
    {
    //
    //  Initialize internal quantities
    //
    N_sites = Lx*Ly;
    N_par = N_up + N_dn;
    //  form the one-body kinetic Hamiltonian
    H_k = H_K(Lx,Ly,tx,ty);
    // the matrix of the operator exp(-deltau*K/2)
    Proj_k_half = Exp(-0.5*deltau*H_k);

    //
    //  Initialize the trial wave function and calculate the ensemble's initial 
    //  energy 
    //
    
    // Diagonalize the one-body kinetic Hamiltonian to get the non-interacting 
    // single-particle orbitals:
    Vector E_nonint;
    Matrix psi_nonint;
    EigenValues(H_k, E_nonint, psi_nonint);

    // assemble the non-interacting single-particle orbitals into a Slater 
    // determinant:
    Phi_T = Matrix(N_sites,N_par);
    if(N_up > 0)
        Phi_T.SubMatrix(1,N_sites,1,N_up) = psi_nonint.SubMatrix(1,N_sites,
                                            1,N_up);
    if(N_dn > 0)
        Phi_T.SubMatrix(1,N_sites,N_up+1,N_par) = psi_nonint.SubMatrix(1,
                                                  N_sites,1,N_dn);
    // the kinetic energy of the trial wave function
    Real E_K = 0.0;
    for(int i = 1; i <= N_up; i++)
        E_K += E_nonint(i);
    for(int i = 1; i <= N_dn; i++)
        E_K += E_nonint(i);

    // the potential energy of the trial wave function
    Vector n_r_up(N_sites),
           n_r_dn(N_sites);
    n_r_up = 0.0, n_r_dn = 0.0;

    if(N_up > 0)
        {
        Matrix Phi_T_up = Phi_T.SubMatrix(1,N_sites,1,N_up);
        Matrix Lambda_up = Phi_T_up*Phi_T_up.t();
        n_r_up = Lambda_up.Diagonal();
        }
    if(N_dn > 0)
        {
        Matrix Phi_T_dn = Phi_T.SubMatrix(1,N_sites,N_up+1,N_par);
        Matrix Lambda_dn = Phi_T_dn*Phi_T_dn.t();
        n_r_dn = Lambda_dn.Diagonal();
        }

    Real E_V = U*n_r_up*n_r_dn;
    // the total energy of the trial wave function = the initial trial energy
    Real E_T = E_K + E_V;

    //
    // Assemble the initial population of walkers
    //

    // initiate each walker to be the trial wave function
    for(int i = 1; i <= N_wlk; i++)
        {
        // Phi[i] is the ith walker. Each is a matrix of size N_sites by N_par
        // The left N_sites by N_up block is the spin up sector
        // The rest is the spin down sector
        // They are propagated independently and only share the auxiliary field
        Phi[i] = Phi_T;
        }

    // initiate the weight and overlap of each walker to 1
    w = Vector(N_wlk);
    O = Vector(N_wlk);
    w = 1.0, O = 1.0;

    // the arrays that store the energy and weight at each block
    E_blk = Vector(N_blk);
    W_blk = Vector(N_blk);
    E_blk = 0.0, W_blk = 0.0;

    //
    // initialize auxiliary field constants
    //

    // exponent of the prefactor exp(-deltau*(-E_T)) in the ground state 
    // projector 
    // fac_norm also include -0.5*U*(N_up+N_dn), the exponent of the 
    // prefactor in the Hirsch transformation
    fac_norm = (E_T-0.5*U*N_par)*deltau;
    // gamma in Hirsch's transformation
    gamma = acosh(exp(0.5*deltau*U));
    // aux_fld is the 2x2 matrix containing all the possible values of 
    // the quantity exp(-gamma*s(sigma)*x_i)
    aux_fld = Matrix(2,2);
    aux_fld = 0.0;
    
    // The first index corresponds to spin up or down
    // The second index corresponds to the auxiliary field x_i=1 or x_i=-1
    for(int i = 1; i <= 2; i++)
        {
        for(int j = 1; j <= 2; j++)
            {
            if((i+j)%2 == 0)
                aux_fld(i,j)=exp(gamma);
            else
                aux_fld(i,j)=exp(-gamma);
            }
        }
  
    } // End initialization()

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
         int N_eqblk, int N_blk, int itv_modsvd, int itv_pc, int itv_Em)
    {
    //
    //  Initialization
    //
    int N_sites, N_par;
    Real fac_norm, gamma;
    Vector w, O, E_blk, W_blk;
    Matrix H_k, Proj_k_half, Phi_T, aux_fld;
    std::vector<Matrix> Phi(N_wlk+1);
    // initialize internal constants, form the trial wave function and assemble the initial 
    // population of walkers
    initialization(Lx, Ly, tx, ty, N_up, N_dn, deltau, U, N_wlk, N_blk, N_sites, N_par, 
                   H_k, Proj_k_half, Phi_T, Phi, w, O, E_blk, W_blk, fac_norm, 
                   gamma, aux_fld);

    // randomize the random number generator seed based on the current time
    srand(time(NULL));

    int flag_mea = 0; // determine when a measurement should take place
    Real E = 0.0;
    Real W = 0.0;

    //
    // Equilibration phase
    //
    for(int i_blk = 1; i_blk <= N_eqblk; i_blk++)
        {
        for(int j_step = 1; j_step <= N_blksteps; j_step++)
            {
            stepwlk(Phi, N_wlk, N_sites, w, O, E, W, H_k, Proj_k_half, flag_mea, Phi_T, 
                    N_up, N_par, U, fac_norm, aux_fld);
            if(j_step % itv_modsvd == 0)
                stblz(Phi, N_wlk, O, N_up, N_par); // re-orthonormalize the walkers
            if(j_step % itv_pc == 0)
                pop_cntrl(Phi, w, O, N_wlk, N_sites, N_par); // population control
            }
        }
    
    // 
    // Measurement phase
    //
    for(int i_blk = 1; i_blk <= N_blk; i_blk++)
        {
        for(int j_step = 1; j_step <= N_blksteps; j_step++)
            {
            if(j_step % itv_Em == 0)
                flag_mea = 1;
            else
                flag_mea = 0;
            // propagate the walkers:
            stepwlk(Phi, N_wlk, N_sites, w, O, E_blk(i_blk), W_blk(i_blk), H_k, 
                    Proj_k_half, flag_mea, Phi_T, N_up, N_par, U, fac_norm, aux_fld);
            if (j_step % itv_modsvd == 0)
                stblz(Phi, N_wlk, O, N_up, N_par); // re-orthonormalize the walkers
            if (j_step % itv_pc == 0)
                pop_cntrl(Phi, w, O, N_wlk, N_sites, N_par); // population control
            if (j_step % itv_Em ==0)
                {
                // update the exponent of the pre-factor exp(-deltau*(H-E_T))
                fac_norm = (E_blk(i_blk)/W_blk(i_blk)-0.5*U*N_par)*deltau;
                }
            }
        E_blk(i_blk)=E_blk(i_blk)/W_blk(i_blk);
        }
   
    //
    // Results
    //
    for(int i = 1; i <= N_blk; i++)
        E_ave += E_blk(i);
    E_ave = E_ave/N_blk;
    if(N_blk > 1)
        {
        for(int i = 1; i <= N_blk; i++)
            E_err += (E_blk(i) - E_ave)*(E_blk(i) - E_ave);
        E_err = std::sqrt(E_err) / std::sqrt(N_blk - 1);
        E_err = E_err / std::sqrt(N_blk);
        }
    } // End CPMC_Lab

};

