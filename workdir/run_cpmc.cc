#include "core.h"
#include "cpmc.h"

using std::cout;
using std::endl;
using std::string;
using namespace itensor;

#define rNum()    ( (double)rand() / (1.0+(double)RAND_MAX) )

void
halfK(Matrix& phi, Real& w, Real& O, Matrix& invO_matrix_up, Matrix& invO_matrix_dn, Matrix Proj_k_half, Matrix Phi_T, int N_up, int N_par)
    {
    phi = Proj_k_half*phi;
    
    int N_site = phi.Nrows();
   
    Matrix Phi_T_up;
    Matrix phi_up;
    Matrix Phi_T_dn;
    Matrix phi_dn;
    Real detinvO_matrix_up = 1.0;
    Real detinvO_matrix_dn = 1.0;

    if(N_up > 0)
        {
        Phi_T_up = Phi_T.SubMatrix(1,N_site,1,N_up);
        phi_up = phi.SubMatrix(1,N_site,1,N_up);
        invO_matrix_up = Inverse(Phi_T_up.t()*phi_up);
        cout << invO_matrix_up;
        detinvO_matrix_up = Determinant(invO_matrix_up);
        }
    
    if(N_par-N_up > 0)
        {
        Phi_T_dn = Phi_T.SubMatrix(1,N_site,N_up+1,N_par);
        phi_dn = phi.SubMatrix(1,N_site,N_up+1,N_par);
        invO_matrix_dn = Inverse(Phi_T_dn.t()*phi_dn);
        detinvO_matrix_dn = Determinant(invO_matrix_dn);
        }
    
    Real O_new = 1.0/(detinvO_matrix_up*detinvO_matrix_dn);
    Real O_ratio=O_new/O;

    if(O_ratio > 0)
        {
        O=O_new;
        w=w*O_ratio;
        }
    else
        {
        w = 0;
        }

    }

void
V(Vector& phi, Vector phi_T, int N_up, int N_par, Real& O, Real& w, Matrix& invO_matrix_up, Matrix& invO_matrix_dn, Matrix aux_fld)
    {
    Vector Gii(2);
    Matrix RR(2,2);
    Gii = 0.0;
    RR = 0.0;

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

    Vector O_ratio_temp(2);
    O_ratio_temp(1) = RR(1,1)*RR(2,1);
    O_ratio_temp(2) = RR(1,2)*RR(2,2);

    Vector O_ratio_temp_real(2);
    O_ratio_temp_real = 0.0;
    if(O_ratio_temp(1) > 0.0)
        O_ratio_temp_real(1) = O_ratio_temp(1);
    if(O_ratio_temp(2) > 0.0)
        O_ratio_temp_real(2) = O_ratio_temp(2);
    Real sum_O_ratio_temp_real = O_ratio_temp_real(1)+O_ratio_temp_real(2);

    if(sum_O_ratio_temp_real <= 0)
        w=0;

    if(w > 0)
        {
        w=w*0.5*sum_O_ratio_temp_real;
        
        int x_spin;
        if(O_ratio_temp_real(1)/sum_O_ratio_temp_real >= rNum())
            x_spin=1;
        else
            x_spin=2;
        if(N_up > 0)
            phi.SubVector(1,N_up)=phi.SubVector(1,N_up)*aux_fld(1,x_spin);
        if(N_par-N_up > 0)
            phi.SubVector(N_up+1,N_par)=phi.SubVector(N_up+1,N_par)*aux_fld(2,x_spin);
    
        O=O*O_ratio_temp(x_spin);
        invO_matrix_up=invO_matrix_up+(1-aux_fld(1,x_spin))/RR(1,x_spin)*temp2_up*temp1_up;
        invO_matrix_dn=invO_matrix_dn+(1-aux_fld(2,x_spin))/RR(2,x_spin)*temp2_dn*temp1_dn;
        }
    }

Real
measure(Matrix H_k, Matrix phi, Matrix Phi_T, Matrix invO_matrix_up, Matrix invO_matrix_dn, int N_up, int N_par, Real U)
    {
    Real e = 0.0;
    int N_sites = phi.Nrows();

    Vector diag_G_up(N_sites);
    diag_G_up = 0.0;
    Vector diag_G_dn(N_sites);
    diag_G_dn = 0.0;
    Matrix G_up(N_sites,N_sites);
    Matrix G_dn(N_sites,N_sites);
    G_up = 0.0;
    G_dn = 0.0;

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
    
    Real kineticEnergy = 0.0;
    for(int i = 1; i <= N_sites; i++)
        {
        for(int j = 1; j <= N_sites; j++)
            {
            kineticEnergy += H_k(i,j)*(G_up(i,j)+G_dn(i,j));
            }
        }
    
    e = potentialEnergy + kineticEnergy;

    return e;
    }

void
stepwlk(std::vector<Matrix>& phi, int N_wlk, int N_sites, Vector& w, Vector& O, Real& E, Real& W, Matrix H_k, Matrix Proj_k_half, int flag_mea, Matrix Phi_T, int N_up, int N_par, Real U, Real fac_norm, Matrix aux_fld)
    {
    Vector e(N_wlk);
    e = 0.0;
    Matrix invO_matrix_up;
    Matrix invO_matrix_dn;

    for(int i_wlk = 1; i_wlk <= N_wlk; i_wlk++)
        {
        Matrix Phi = phi[i_wlk];
        if(w(i_wlk) > 0)
            {
            w(i_wlk) = w(i_wlk)*exp(fac_norm);
            halfK(Phi, w(i_wlk), O(i_wlk), invO_matrix_up, invO_matrix_dn, Proj_k_half, Phi_T, N_up, N_par);
            if(w(i_wlk) > 0)
                {
                for(int j_site = 1; j_site <= N_sites; j_site++)
                    {
                    if(w(i_wlk) > 0)
                        {
                        Vector Phi_site = Phi.Row(j_site);
                        V(Phi_site, Phi_T.Row(j_site), N_up, N_par, O(i_wlk), w(i_wlk), invO_matrix_up, invO_matrix_dn, aux_fld);
                        Phi.Row(j_site) = Phi_site;
                        }
                    }
                }
            if(w(i_wlk) > 0)
                {
                halfK(Phi, w(i_wlk), O(i_wlk), invO_matrix_up, invO_matrix_dn, Proj_k_half, Phi_T, N_up, N_par);
                if(w(i_wlk) > 0)
                    {
                    if(flag_mea == 1)
                        {
                        e(i_wlk) = measure(H_k, Phi, Phi_T,  invO_matrix_up, invO_matrix_dn, N_up, N_par, U); 
                        }
                    }
                }
            }
        phi[i_wlk] = Phi;
        }

    // Compute the ensemble's total energy and weight if measurement took place
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
    }

void
stblz(std::vector<Matrix>& Phi, int N_wlk, Vector& O, int N_up, int N_par)
    {
    int N_sites = Phi[1].Nrows();
    int N_dn = N_par - N_up;

    for(int i_wlk = 1; i_wlk <= N_wlk; i_wlk++)
        {
        Matrix phi = Phi[i_wlk];
        Real det_R_up = 1.0;
        Real det_R_dn = 1.0;

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
        O(i_wlk) = O(i_wlk)/det_R_up/det_R_dn;
        }
    }

void
pop_cntrl(std::vector<Matrix>& Phi, Vector& w, Vector& O, int N_wlk, int N_sites, int N_par)
    {
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
    Real d = 1.0*N_wlk/sum_w;
    
    sum_w = -rNum();
    int n_wlk = 0;

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

    Phi = new_Phi;
    O = new_O;

    w = Vector(N_wlk,1.0);
    }

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

int 
main(int argc, char* argv[])
    {
    //
    // system parameters:
    //
    int Lx = 4;     // The number of lattice sites in the x direction
    int Ly = 1;     // The number of lattice sites in the y direction

    int N_up = 2;   // The number of spin-up electrons
    int N_dn = 0;   // The number of spin-down electrons

    Real U = 0.0;   // The on-site repulsion strength in the Hubbard Hamiltonian
    Real tx = 1.0;  // The hopping amplitude between nearest-neighbor sites in the 
                    // x direction
    Real ty = 1.0;  // The hopping amplitude between nearest neighbor sites in the 
                    // y direction

    //
    // run parameters:
    //
    Real deltau = 0.01;     // The imaginary time step
    int N_wlk = 1;          // The number of random walkers
    int N_blksteps = 1;     // The number of random walk steps in each block
    int N_eqblk = 1;        // The number of blocks used to equilibrate the random 
                            // walk before energy measurement takes place
    int N_blk = 1;          // The number of blocks used in the measurement phase
    int itv_modsvd = 1;     // The interval between two adjacent modified Gram-Schmidt 
                            // re-orthonormalization of the random walkers. No 
                            // re-orthonormalization if itv_modsvd > N_blksteps
    int itv_pc = 1;         // The interval between two adjacent population controls. 
                            // No population control if itv_pc > N_blksteps
    int itv_Em = 1;         // The interval between two adjacent energy measurements
   
    //
    // invoke the main function
    //
    int N_sites = Lx*Ly;
    int N_par = N_up + N_dn;

    Matrix H_k(N_sites,N_sites);
    H_k = 0.0;

    int r = 0;
    for(int iy = 1; iy <= Ly; iy++)
        {
        for(int jx = 1; jx <= Lx; jx++)
            {
            r++;
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
   
    Matrix Proj_k_half = Exp(-0.5*deltau*H_k);

    int n = N_sites;

    Vector E_nonint;
    Matrix psi_nonint;

    EigenValues(H_k, E_nonint, psi_nonint);

    Matrix Phi_T(N_sites,N_par);
    if(N_up > 0)
        {
        Phi_T.SubMatrix(1,N_sites,1,N_up) = psi_nonint.SubMatrix(1,N_sites,1,N_up);
        }
    if(N_dn > 0)
        {
        Phi_T.SubMatrix(1,N_sites,N_up+1,N_par) = psi_nonint.SubMatrix(1,N_sites,1,N_dn);;
        }

    Real E_K = 0.0;

    for(int i = 1; i <= N_up; i++)
        E_K += E_nonint(i);
    for(int i = 1; i <= N_dn; i++)
        E_K += E_nonint(i);

    Vector n_r_up(N_sites);
    Vector n_r_dn(N_sites);
    n_r_up = 0.0;
    n_r_dn = 0.0;

    if(N_up > 0)
        {
        Matrix Lambda_up = Phi_T.SubMatrix(1,N_sites,1,N_up)*(Phi_T.SubMatrix(1,N_sites,1,N_up)).t();
        n_r_dn = Lambda_up.Diagonal();
        }
    if(N_dn > 0)
        {
        Matrix Lambda_dn = Phi_T.SubMatrix(1,N_sites,N_up+1,N_par)*(Phi_T.SubMatrix(1,N_sites,N_up+1,N_par)).t();
        n_r_dn = Lambda_dn.Diagonal();
        }

    Real E_V = U*n_r_up*n_r_dn;
    Real E_T = E_K + E_V;

    // Assemble the initial population of walkers
    std::vector<Matrix> Phi(N_wlk+1);
    for(int i = 1; i <= N_wlk; i++)
        {
        Phi[i] = Phi_T;
        }

    Vector w(N_wlk);
    Vector O(N_wlk);

    w = 1.0, O = 1.0;

    Vector E_blk(N_blk);
    Vector W_blk(N_blk);

    E_blk = 0.0, W_blk = 0.0;

    // initialize auxiliary field constants
    Real fac_norm = (E_T-0.5*U*N_par)*deltau;
    Real gamma = acosh(exp(0.5*deltau*U));
    Matrix aux_fld(2,2);
    aux_fld = 0.0;

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
   
    srand(time(NULL));

    int flag_mea = 0;
    Real E = 0.0;
    Real W = 0.0;

    // Equilibration phase
    for(int i_blk = 1; i_blk <= N_eqblk; i_blk++)
        {
        for(int j_step = 1; j_step <= N_blksteps; j_step++)
            {
            stepwlk(Phi, N_wlk, N_sites, w, O, E, W, H_k, Proj_k_half, flag_mea, Phi_T, N_up, N_par, U, fac_norm, aux_fld);
            if(j_step % itv_modsvd == 0)
                {
                stblz(Phi, N_wlk, O, N_up, N_par);
                }
            if(j_step % itv_pc == 0)
                {
                pop_cntrl(Phi, w, O, N_wlk, N_sites, N_par);
                }
            }
        }
    
    // Measurement phase
    for(int i_blk = 1; i_blk <= N_blk; i_blk++)
        {
        for(int j_step = 1; j_step <= N_blksteps; j_step++)
            {
            if(j_step % itv_Em == 0)
                flag_mea = 1;
            else
                flag_mea = 0;
            stepwlk(Phi, N_wlk, N_sites, w, O, E_blk(i_blk), W_blk(i_blk), H_k, Proj_k_half, flag_mea, Phi_T, N_up, N_par, U, fac_norm, aux_fld);
            if (j_step % itv_modsvd == 0)
                {
                stblz(Phi, N_wlk, O, N_up, N_par);
                }
            if (j_step % itv_pc == 0)
                pop_cntrl(Phi, w, O, N_wlk, N_sites, N_par);
            if (j_step % itv_Em ==0)
                fac_norm = (E_blk(i_blk)/W_blk(i_blk)-0.5*U*N_par)*deltau;
            }
        E_blk(i_blk)=E_blk(i_blk)/W_blk(i_blk);
        }
    
    Real E_ave = 0.0;
    for(int i = 1; i <= N_blk; i++)
        E_ave += E_blk(i);
    E_ave = E_ave/N_blk;
    Real E_err = 0.0;
    if(N_blk > 1)
        {
        for(int i = 1; i <= N_blk; i++)
            E_err += (E_blk(i) - E_ave)*(E_blk(i) - E_ave);
        E_err = sqrt(E_err) / sqrt(N_blk - 1);
        E_err = E_err / sqrt(N_blk);
        }

    cout << "E_ave = " << E_ave << endl;
    cout << "E_err = " << E_err << endl;
    cout << endl;
    
    return 0;
    }


