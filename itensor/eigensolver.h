//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_EIGENSOLVER_H
#define __ITENSOR_EIGENSOLVER_H
#include "iqcombiner.h"


/* Notes on optimization
 * 
 * - May be able to avoid additional
 *   work on the last iteration.
 *   See line 140 of david.cc
 *
 * - Is the MatrixRef Davidon
 *   actually doing fewer iterations?
 *
 */

class Eigensolver
    {
    public:

    Eigensolver(int maxiter = 2, Real errgoal = 1E-4, int numget = 1);

    template <class SparseT, class Tensor> Real 
    davidson(const SparseT& A, Tensor& phi) const;

    //Accessor methods ------------

    Real 
    errgoal() const { return errgoal_; }
    void 
    errgoal(Real val) { errgoal_ = val; }

    int 
    numGet() const { return numget_; }
    void 
    numGet(int val) { numget_ = val; }

    int 
    maxIter() const { return maxiter_; }
    void 
    maxIter(int val) { maxiter_ = val; }

    int 
    debugLevel() const { return debug_level_; }
    void 
    debugLevel(int val) { debug_level_ = val; }

    //Other methods ------------

    private:

    //Function object which applies the mapping
    // f(x,theta) = 1/(theta - x)
    class DavidsonPrecond
        {
        public:
            DavidsonPrecond(Real theta)
                : theta_(theta)
                { }
            Real
            operator()(Real val) const
                {
                if(theta_ == val)
                    return 0;
                else
                    return 1.0/(theta_-val);
                }
        private:
            Real theta_;
        };

    class LanczosPrecond
        {
        public:
            LanczosPrecond(Real theta)
                : theta_(theta)
                { }
            Real
            operator()(Real val) const
                {
                return 1.0/(theta_-1+1E-33);
                }
        private:
            Real theta_;
        };

    int maxiter_;
    Real errgoal_;
    int numget_;
    int debug_level_;

    }; //class Eigensolver


inline Eigensolver::
Eigensolver(int maxiter, Real errgoal, int numget)
    : maxiter_(maxiter),
      errgoal_(errgoal),
      numget_(numget),
      debug_level_(-1)
    { }

template <class SparseT, class Tensor> 
inline Real Eigensolver::
davidson(const SparseT& A, Tensor& phi) const
    {
    phi *= 1.0/phi.norm();

    const int maxsize = A.size();
    const int actual_maxiter = min(maxiter_,maxsize);
    Real lambda = 1E30, 
         last_lambda = lambda,
         qnorm = 1E30;

    std::vector<Tensor> B(actual_maxiter+2),
                       AB(actual_maxiter+2);

    //Storage for Matrix that gets diagonalized 
    Matrix M(actual_maxiter+2,actual_maxiter+2);

    MatrixRef Mref(M.SubMatrix(1, 1, 1, 1));

    //Get diagonal of A to use later
    Tensor Adiag(phi);
    A.diag(Adiag);

    for(int ii = 1; ii <= actual_maxiter; ++ii)
        {
        //Diagonalize conj(B)*A*B
        //and compute the residual q
        Tensor q;
        if(ii == 1)
            {
            B[1] = phi;
            A.product(B[1],AB[1]);

            //No need to diagonalize
            lambda = Dot(B[1],AB[1]);
            Mref = lambda;

            //Calculate residual q
            q = B[1];
            q *= -lambda;
            q += AB[1]; 
            }
        else // ii != 1
            {
            //Diagonalize M
            Vector D;
            Matrix U;
            EigenValues(Mref,D,U);

            //lambda is the minimum eigenvalue of M
            lambda = D(1);

            //Compute corresponding eigenvector
            //phi of A from the min evec of M
            //(and start calculating residual q)
            phi = U(1,1)*B[1];
            q   = U(1,1)*AB[1];
            for(int k = 2; k <= ii; ++k)
                {
                phi += U(k,1)*B[k];
                q   += U(k,1)*AB[k];
                }

            //Calculate residual q
            q += (-lambda)*phi;
            }

        //Check convergence
        qnorm = q.norm();
        if( (qnorm < errgoal_ && fabs(lambda-last_lambda) < errgoal_) 
            || qnorm < 1E-12 )
            {
            if(debug_level_ > 0)
                {
                std::cout << boost::format("I %d q %.0E E %.10f")
                             % ii
                             % qnorm
                             % lambda 
                             << std::endl;
                }
            return lambda;
            }

        /*
        if(debug_level_ >= 1 && qnorm > 3)
            {
            std::cerr << "Large qnorm = " << qnorm << "\n";
            Print(Mref);
            Vector D;
            Matrix U;
            EigenValues(Mref,D,U);
            Print(D);

            std::cerr << "ii = " << ii  << "\n";
            Matrix Borth(ii,ii);
            for(int i = 1; i <= ii; ++i)
            for(int j = i; j <= ii; ++j)
                {
                Borth(i,j) = Dot(B[i],B[j]);
                Borth(j,i) = Borth(i,j);
                }
            Print(Borth);
            exit(0);
            }
        */

        if(debug_level_ > 1 || (ii == 1 && debug_level_ > 0))
            {
            std::cout << boost::format("I %d q %.0E E %.10f")
                         % ii
                         % qnorm
                         % lambda 
                         << std::endl;
            }

        //Apply Davidson preconditioner
        {
        DavidsonPrecond dp(lambda);
        Tensor cond(Adiag);
        cond.mapElems(dp);
        q /= cond;
        }

        //Do Gram-Schmidt on xi
        //to include it in the subbasis
        Tensor& d = B[ii+1];
        d = q;
        Vector Bd(ii);
        for(int k = 1; k <= ii; ++k)
            {
            Bd(k) = Dot(B[k],d);
            }
        d = Bd(1)*B[1];
        for(int k = 2; k <= ii; ++k)
            {
            d += Bd(k)*B[k];
            }
        d *= -1;
        d += q;
        d *= 1.0/(d.norm()+1E-33);

        last_lambda = lambda;

        //Expand AB and M
        //for next step
        if(ii < actual_maxiter)
            {
            A.product(d,AB[ii+1]);

            //Add new row and column to M
            Mref << M.SubMatrix(1,ii+1,1,ii+1);
            Vector newCol(ii+1);
            for(int k = 1; k <= ii+1; ++k)
                {
                newCol(k) = Dot(B[k],AB[ii+1]);
                }
            Mref.Column(ii+1) = newCol;
            Mref.Row(ii+1) = newCol;
            }

        } //for(ii)

    if(debug_level_ > 0)
        {
        std::cout << boost::format("I %d q %.0E E %.10f")
                     % actual_maxiter
                     % qnorm
                     % lambda 
                     << std::endl;
        }

    return lambda;

    } //Eigensolver::davidson


#endif
