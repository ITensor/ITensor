#ifndef __CPMC_H_
#define __CPMC_H_

#include <vector>
#include "storelink.h"

namespace itensor {

class Matrix;
class Vector;

void
halfK(Matrix& phi, Real& w, Real& O, Matrix& invO_matrix_up, 
      Matrix& invO_matrix_dn, Matrix Proj_k_half, Matrix Phi_T, 
      int N_up, int N_par);

void
V(Vector& phi, Vector phi_T, int N_up, int N_par, Real& O, Real& w, 
  Matrix& invO_matrix_up, Matrix& invO_matrix_dn, Matrix aux_fld);

Real
measure(Matrix H_k, Matrix phi, Matrix Phi_T, Matrix invO_matrix_up, 
        Matrix invO_matrix_dn, int N_up, int N_par, Real U);

void
stepwlk(std::vector<Matrix>& phi, int N_wlk, int N_sites, Vector& w, 
        Vector& O, Real& E, Real& W, Matrix H_k, Matrix Proj_k_half, 
        int flag_mea, Matrix Phi_T, int N_up, int N_par, Real U, 
        Real fac_norm, Matrix aux_fld);

void
stblz(std::vector<Matrix>& Phi, int N_wlk, Vector& O, int N_up, int N_par);

void
pop_cntrl(std::vector<Matrix>& Phi, Vector& w, Vector& O, int N_wlk, 
          int N_sites, int N_par);

};

#endif
