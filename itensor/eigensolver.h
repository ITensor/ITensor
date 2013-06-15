//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_EIGENSOLVER_H
#define __ITENSOR_EIGENSOLVER_H
#include "iqcombiner.h"

#define Cout std::cout
#define Endl std::endl
#define Format boost::format

/*
template<class Tensor>
void
orthog(std::vector<Tensor>& T, int num, int numpass, int start = 1);
*/


/* Notes on optimization
 * 
 * - May be faster to do a direct solution
 *   for small sizes.
 *
 */

class Eigensolver
    {
    public:

    Eigensolver(const OptSet& opts = Global::opts());

    //
    // Uses the Davidson algorithm to find the 
    // minimal eigenvector of the sparse matrix A.
    // (LocalT objects must implement the methods product, size and diag.)
    // Returns the minimal eigenvalue lambda such that
    // A phi = lambda phi.
    //
    template <class LocalT, class Tensor> 
    Real 
    davidson(const LocalT& A, Tensor& phi) const;

    //
    // Uses the Davidson algorithm to find the minimal
    // eigenvector of the generalized eigenvalue problem
    // A phi = lambda B phi.
    // (B should have positive definite eigenvalues.)
    //
    template <class LocalTA, class LocalTB, class Tensor> 
    Real
    genDavidson(const LocalTA& A, const LocalTB& B, Tensor& phi) const;

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
    minIter() const { return miniter_; }
    void 
    minIter(int val) { miniter_ = val; }

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

    //Function object which applies the mapping
    // f(x,theta) = 1/(theta - 1)
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

    //Function object which applies the mapping
    // f(x) = (x < cut ? 0 : 1/x);
    class PseudoInverter
        {
        public:
            PseudoInverter(Real cut = MIN_CUT)
                :
                cut_(cut)
                { }

            Real
            operator()(Real val) const
                {
                if(fabs(val) < cut_)
                    return 0;
                else
                    return 1./val;
                }
        private:
            Real cut_;
        };

    int maxiter_;
    int miniter_;
    Real errgoal_;
    int numget_;
    int debug_level_;

    }; //class Eigensolver


inline Eigensolver::
Eigensolver(const OptSet& opts)
    : 
    miniter_(1)
    { 
    maxiter_ = opts.getInt("MaxIter",2);
    errgoal_ = opts.getReal("ErrGoal",1E-4);
    numget_ = opts.getInt("NumGet",1);
    debug_level_ = opts.getInt("DebugLevel",-1);
    }

template <class LocalT, class Tensor> 
Real inline Eigensolver::
davidson(const LocalT& A, Tensor& phi) const
    {
    typedef typename Tensor::SparseT
    SparseT;

    phi *= 1.0/phi.norm();

    bool complex_diag = false;

    const int maxsize = A.size();
    const int actual_maxiter = min(maxiter_,maxsize-1);
    if(debug_level_ >= 2)
        {
        Cout << Format("maxsize-1 = %d, maxiter = %d, actual_maxiter = %d") 
                % (maxsize-1) % maxiter_ % actual_maxiter << Endl;
        }

    std::vector<Tensor> V(actual_maxiter+2),
                       AV(actual_maxiter+2);

    //Storage for Matrix that gets diagonalized 
    Matrix MR(actual_maxiter+2,actual_maxiter+2),
           MI(actual_maxiter+2,actual_maxiter+2);
    MR = NAN; //set to NAN to ensure failure if we use uninitialized elements
    MI = NAN;

    //Mref holds current projection of A into V's
    MatrixRef MrefR(MR.SubMatrix(1,1,1,1)),
              MrefI(MI.SubMatrix(1,1,1,1));

    //Get diagonal of A to use later
    const Tensor Adiag = A.diag();

    Real lambda = NAN, //current lowest eigenvalue
         last_lambda = NAN,
         qnorm = NAN; //norm of residual, primary convergence criterion

    V[0] = phi;
    A.product(V[0],AV[0]);

    Complex z = BraKet(V[0],AV[0]);
    const Real initEn = z.real();

    if(debug_level_ > 2)
        Cout << Format("Initial Davidson energy = %.10f") % initEn << Endl;

    const Real enshift = 0;
    //const Real enshift = initEn;

    int iter = 0;
    for(int ii = 0; ii <= actual_maxiter; ++ii)
        {
        //Diagonalize conj(V)*A*V
        //and compute the residual q

        const int ni = ii+1; 
        Tensor& q = V.at(ni);

        if(ii == 0)
            {
            lambda = initEn;
            MrefR = lambda-enshift;
            MrefI = 0;
            //Calculate residual q
            q = V[0];
            q *= -lambda;
            q += AV[0]; 
            }
        else // ii != 0
            {
            //Diagonalize M
            Vector D;
            Matrix UR;
            if(complex_diag)
                {
                Matrix UI;
                HermitianEigenvalues(MrefR,MrefI,D,UR,UI);

                //Compute corresponding eigenvector
                //phi of A from the min evec of M
                //(and start calculating residual q)

                phi = (UR(1,1)*Complex_1+UI(1,1)*Complex_i)*V[0];
                q   = (UR(1,1)*Complex_1+UI(1,1)*Complex_i)*AV[0];
                for(int k = 1; k <= ii; ++k)
                    {
                    const Complex cfac = (UR(k+1,1)*Complex_1+UI(k+1,1)*Complex_i);
                    phi += cfac*V[k];
                    q   += cfac*AV[k];
                    }
                }
            else
                {
                EigenValues(MrefR,D,UR);

                //Compute corresponding eigenvector
                //phi of A from the min evec of M
                //(and start calculating residual q)
                phi = UR(1,1)*V[0];
                q   = UR(1,1)*AV[0];
                for(int k = 1; k <= ii; ++k)
                    {
                    phi += UR(k+1,1)*V[k];
                    q   += UR(k+1,1)*AV[k];
                    }
                }

            //lambda is the minimum eigenvalue of M
            lambda = D(1)+enshift;

            //Calculate residual q
            q += (-lambda)*phi;

            //Fix sign
            if(UR(1,1) < 0)
                {
                phi *= -1;
                q *= -1;
                }

            if(debug_level_ >= 3)
                {
                Cout << "complex_diag = " 
                     << (complex_diag ? "true" : "false") << Endl;
                Cout << "D = " << D;
                Cout << Format("lambda = %.10f") % (D(1)+enshift) << Endl;
                }

            }

        //Check convergence
        qnorm = q.norm();

        if(qnorm < 1E-20)
            {
            if(debug_level_ >= 3) //Explain why breaking out of Davidson loop early
                {
                Cout << Format("Breaking out of Davidson because qnorm = %.2E < 1E-20") 
                        % qnorm 
                        << Endl;
                }
            goto done;
            }

        const bool converged = (qnorm < errgoal_ && fabs(lambda-last_lambda) < errgoal_) 
                               || qnorm < max(1E-12,errgoal_ * 1.0e-3);

        if((converged && ii >= miniter_) || (ii == actual_maxiter))
            {
            if(debug_level_ >= 3) //Explain why breaking out of Davidson loop early
                {
                if((qnorm < errgoal_ && fabs(lambda-last_lambda) < errgoal_))
                    Cout << "Breaking out of Davidson because errgoal reached" << Endl;
                else
                if(qnorm < max(1E-12,errgoal_ * 1.0e-3) && ii >= miniter_)
                    Cout << "Breaking out of Davidson because small errgoal obtained" << Endl;
                else
                if(ii == actual_maxiter)
                    Cout << "Breaking out of Davidson because ii == actual_maxiter" << Endl;
                }
            goto done;
            }
        
        if(debug_level_ >= 2 || (ii == 0 && debug_level_ >= 1))
            {
            Cout << Format("I %d q %.0E E %.10f")
                           % ii
                           % qnorm
                           % lambda 
                           << Endl;
            }

        if(ii < actual_maxiter)
            {
            //On all but last step,
            //compute next Krylov/trial vector by
            //first applying Davidson preconditioner
            //formula then orthogonalizing against
            //other vectors

            //Apply Davidson preconditioner
            {
            DavidsonPrecond dp(lambda);
            Tensor cond(Adiag);
            cond.mapElems(dp);
            q /= cond;
            }

            //Do Gram-Schmidt on d (Npass times)
            //to include it in the subbasis
            const int Npass = 2;
            std::vector<Complex> Vq(ni);

            int count = 0;
            for(int pass = 1; pass <= Npass; ++pass)
                {
                ++count;
                for(int k = 0; k < ni; ++k)
                    {
                    Vq[k] = BraKet(V[k],q);
                    }

                for(int k = 0; k < ni; ++k)
                    {
                    q += (-Vq[k].real())*V[k];
                    if(Vq[k].imag() != 0)
                        {
                        q += (-Vq[k].imag()*Complex_i)*V[k];
                        }
                    }

                Real qn = q.norm();

                if(qn < 1E-10)
                    {
                    //Orthogonalization failure,
                    //try randomizing
                    if(debug_level_ >= 2)
                        Cout << "Vector not independent, randomizing" << Endl;
                    q = V.at(ni-1);
                    q.randomize();

                    if(q.indices().dim() <= ni)
                        {
                        //Not be possible to orthogonalize if
                        //max size of q (vecSize after randomize)
                        //is size of current basis
                        if(debug_level_ >= 3)
                            Cout << "Breaking out of Davidson: max Hilbert space size reached" << Endl;
                        goto done;
                        }

                    qn = q.norm();
                    --pass;
                    }

                q *= 1./qn;
                }


            /*
            if(debug_level_ >= 4)
                {
                //Check V's are orthonormal
                Matrix Vo(ni+1,ni+1); 
                Vo = NAN;
                for(int r = 1; r <= ni+1; ++r)
                for(int c = r; c <= ni+1; ++c)
                    {
                    Vo(r,c) = Dot(conj(V[r-1]),V[c-1]);
                    Vo(c,r) = Vo(r,c);
                    }
                Print(Vo);
                }
            */

            last_lambda = lambda;

            //Expand AV and M
            //for next step
            A.product(V[ni],AV[ni]);

            //Add new row and column to M
            MrefR << MR.SubMatrix(1,ni+1,1,ni+1);
            MrefI << MI.SubMatrix(1,ni+1,1,ni+1);
            Vector newColR(ni+1),
                   newColI(ni+1);
            for(int k = 0; k <= ni; ++k)
                {
                z = BraKet(V.at(k),AV.at(ni));
                newColR(k+1) = z.real();
                newColI(k+1) = z.imag();
                }
            newColR(ni+1) -= enshift;

            MrefR.Column(ni+1) = newColR;
            MrefR.Row(ni+1) = newColR;

            if(!complex_diag && Norm(newColI) > errgoal_)
                { complex_diag = true; }

            MrefI.Column(ni+1) = newColI;
            MrefI.Row(ni+1) = -newColI;

            } //if ii < actual_maxiter

        ++iter;

        } //for(ii)

    done:

    if(debug_level_ >= 3)
        {
        //Check V's are orthonormal
        Matrix Vo_final(iter+1,iter+1); 
        Vo_final = NAN;
        for(int r = 1; r <= iter+1; ++r)
        for(int c = r; c <= iter+1; ++c)
            {
            z = BraKet(V[r-1],V[c-1]);
            Vo_final(r,c) = z.real();
            Vo_final(c,r) = Vo_final(r,c);
            }
        Print(Vo_final);
        }


    if(debug_level_ > 0)
        {
        Cout << Format("I %d q %.0E E %.10f")
                       % iter
                       % qnorm
                       % lambda 
                       << Endl;
        }

    return lambda;

    } //Eigensolver::davidson

template <class LocalTA, class LocalTB, class Tensor> 
inline Real Eigensolver::
genDavidson(const LocalTA& A, const LocalTB& B, Tensor& phi) const
    {
    typedef typename Tensor::SparseT
    SparseT;

    //B-normalize phi
    {
    Tensor Bphi;
    B.product(phi,Bphi);
    Real phiBphi = Dot(conj(phi),Bphi);
    phi *= 1.0/sqrt(phiBphi);
    }


    const int maxsize = A.size();
    const int actual_maxiter = min(maxiter_,maxsize);
    Real lambda = 1E30, 
         last_lambda = lambda,
         qnorm = 1E30;

    std::vector<Tensor> V(actual_maxiter+2),
                       AV(actual_maxiter+2),
                       BV(actual_maxiter+2);

    //Storage for Matrix that gets diagonalized 
    //M is the projected form of A
    //N is the projected form of B
    Matrix M(actual_maxiter+2,actual_maxiter+2),
           N(actual_maxiter+2,actual_maxiter+2);

    MatrixRef Mref(M.SubMatrix(1,1,1,1)),
              Nref(N.SubMatrix(1,1,1,1));

    Vector D;
    Matrix U;

    //Get diagonal of A,B to use later
    //Tensor Adiag(phi);
    //A.diag(Adiag);
    //Tensor Bdiag(phi);
    //B.diag(Bdiag);

    int iter = 0;
    for(int ii = 1; ii <= actual_maxiter; ++ii)
        {
        ++iter;
        //Diagonalize conj(V)*A*V
        //and compute the residual q
        Tensor q;
        if(ii == 1)
            {
            V[1] = phi;
            A.product(V[1],AV[1]);
            B.product(V[1],BV[1]);

            //No need to diagonalize
            Mref = Dot(conj(V[1]),AV[1]);
            Nref = Dot(conj(V[1]),BV[1]);
            lambda = Mref(1,1)/(Nref(1,1)+1E-33);

            //Calculate residual q
            q = BV[1];
            q *= -lambda;
            q += AV[1]; 
            }
        else // ii != 1
            {
            //Diagonalize M
            //Print(Mref);
            //std::cerr << "\n";
            //Print(Nref);
            //std::cerr << "\n";

            GeneralizedEV(Mref,Nref,D,U);

            //lambda is the minimum eigenvalue of M
            lambda = D(1);

            //Calculate residual q
            q   = U(1,1)*AV[1];
            q   = U(1,1)*(AV[1]-lambda*BV[1]);
            for(int k = 2; k <= ii; ++k)
                {
                q = U(k,1)*(AV[k]-lambda*BV[k]);
                }
            }

        //Check convergence
        qnorm = q.norm();
        if( (qnorm < errgoal_ && fabs(lambda-last_lambda) < errgoal_) 
            || qnorm < 1E-12 )
            {
            break; //Out of ii loop to return
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
            Matrix Vorth(ii,ii);
            for(int i = 1; i <= ii; ++i)
            for(int j = i; j <= ii; ++j)
                {
                Vorth(i,j) = Dot(conj(V[i]),V[j]);
                Vorth(j,i) = Vorth(i,j);
                }
            Print(Vorth);
            exit(0);
            }
        */

        if(debug_level_ > 1 || (ii == 1 && debug_level_ > 0))
            {
            Cout << Format("I %d q %.0E E %.10f")
                           % ii
                           % qnorm
                           % lambda 
                           << Endl;
            }

        //Apply generalized Davidson preconditioner
        //{
        //Tensor cond = lambda*Bdiag - Adiag;
        //PseudoInverter inv;
        //cond.mapElems(inv);
        //q /= cond;
        //}

        /*
         * According to Kalamboukis Gram-Schmidt not needed,
         * presumably because the new entries of N will
         * contain the B-overlap of any new vectors and thus
         * account for their non-orthogonality.
         * Since B-orthogonalizing new vectors would be quite 
         * expensive, what about doing regular Gram-Schmidt?
         * The idea is it won't hurt if it's a little wrong
         * since N will account for it, but if B is close
         * to the identity then it should help a lot.
         *
         */

//#define DO_GRAM_SCHMIDT

        //Do Gram-Schmidt on xi
        //to include it in the subbasis
        Tensor& d = V[ii+1];
        d = q;
#ifdef  DO_GRAM_SCHMIDT
        Vector Vd(ii);
        for(int k = 1; k <= ii; ++k)
            {
            Vd(k) = Dot(conj(V[k]),d);
            }
        d = Vd(1)*V[1];
        for(int k = 2; k <= ii; ++k)
            {
            d += Vd(k)*V[k];
            }
        d *= -1;
        d += q;
#endif//DO_GRAM_SCHMIDT
        d *= 1.0/(d.norm()+1E-33);

        last_lambda = lambda;

        //Expand AV, M and BV, N
        //for next step
        if(ii < actual_maxiter)
            {
            A.product(d,AV[ii+1]);
            B.product(d,BV[ii+1]);


            //Add new row and column to N
            Vector newCol(ii+1);
            Nref << N.SubMatrix(1,ii+1,1,ii+1);
            for(int k = 1; k <= ii+1; ++k)
                {
                newCol(k) = Dot(conj(V[k]),BV[ii+1]);

                if(newCol(k) < 0)
                    {
                    //if(k > 1)
                        //Error("Can't fix sign of new basis vector");
                    newCol(k) *= -1;
                    BV[ii+1] *= -1;
                    d *= -1;
                    }
                }
            Nref.Column(ii+1) = newCol;
            Nref.Row(ii+1) = newCol;

            //Add new row and column to M
            Mref << M.SubMatrix(1,ii+1,1,ii+1);
            for(int k = 1; k <= ii+1; ++k)
                {
                newCol(k) = Dot(conj(V[k]),AV[ii+1]);
                }
            Mref.Column(ii+1) = newCol;
            Mref.Row(ii+1) = newCol;

            }

        } //for(ii)

    if(debug_level_ > 0)
        {
        Cout << Format("I %d q %.0E E %.10f")
                       % iter
                       % qnorm
                       % lambda 
                       << Endl;
        }

    //Compute eigenvector phi before returning
#ifdef DEBUG
    if(U.Nrows() != iter)
        {
        Print(U.Nrows());
        Print(iter);
        Error("Wrong size: U.Nrows() != iter");
        }
#endif
    phi = U(1,1)*V[1];
    for(int k = 2; k <= iter; ++k)
        {
        phi += U(k,1)*V[k];
        }

    return lambda;

    } //Eigensolver::genDavidson

/*
template<class Tensor>
void
orthog(std::vector<Tensor>& T, int num, int numpass, int start)
    {
    const int size = T[start].maxSize();
    if(num > size)
        {
        Print(num);
        Print(size);
        Error("num > size");
        }

    for(int n = start; n <= num+(start-1); ++n)
        { 
        Tensor& col = T.at(n);
        Real norm = col.norm();
        if(norm == 0)
            {
            col.randomize();
            norm = col.norm();
            //If norm still zero, may be
            //an IQTensor with no blocks
            if(norm == 0)
                {
                PrintDat(col);
                Error("Couldn't randomize column");
                }
            }
        col /= norm;
        col.scaleTo(1);
        
        if(n == start) continue;

        for(int pass = 1; pass <= numpass; ++pass)
            {
            Vector dps(n-start);
            for(int m = start; m < n; ++m)
                {
                dps(m-start+1) = Dot(conj(col),T.at(m));
                }
            Tensor ovrlp = dps(1)*T.at(start);
            for(int m = start+1; m < n; ++m)
                {
                ovrlp += dps(m-start+1)*T.at(m);
                }
            ovrlp *= -1;
            col += ovrlp;

            norm = col.norm();
            if(norm == 0)
                {
                Error("Couldn't normalize column");
                }
            col /= norm;
            col.scaleTo(1);

            }
        }
    } // orthog(vector<Tensor> ... )
    */

#undef Cout
#undef Endl
#undef Format

#endif
