//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_EIGENSOLVER_H
#define __ITENSOR_EIGENSOLVER_H
#include "iqtensor.h"
#include "matrix/algs.h"


namespace itensor {

//
// Use the Davidson algorithm to find the 
// eigenvector of the Hermitian matrix A with minimal eigenvalue.
// (BigMatrixT objects must implement the methods product, size and diag.)
// Returns the minimal eigenvalue lambda such that
// A phi = lambda phi.
//
template <class BigMatrixT, class Tensor> 
Real 
davidson(const BigMatrixT& A, Tensor& phi,
         const Args& args = Global::args());

//
// Use Davidson to find the N eigenvectors with smallest 
// eigenvalues of the Hermitian matrix A, given a vector of N 
// initial guesses (zero indexed).
// (BigMatrixT objects must implement the methods product, size and diag.)
// Returns a vector of the N smallest eigenvalues corresponding
// to the set of eigenvectors phi.
//
template <class BigMatrixT, class Tensor> 
std::vector<Real>
davidson(const BigMatrixT& A, 
         std::vector<Tensor>& phi,
         const Args& args = Global::args());

template <class BigMatrixT, class Tensor> 
std::vector<Complex>
complexDavidson(const BigMatrixT& A, 
                std::vector<Tensor>& phi,
                const Args& args = Global::args());




//
//
// Implementations
//
//


int inline
findEig(int which,        //zero-indexed; so is return value
        const Vec& DR, //real part of eigenvalues
        const Vec& DI) //imag part of eigenvalues
    {
    const auto L = DR.size();
#ifdef DEBUG
    if(DI.size() != L) Error("Vectors must have same length in findEig");
#endif
    Vec A2(L);
    int n = -1;

    //Find maximum norm2 of all the eigs
    Real maxj = -1;
    for(int ii = 1; ii <= L; ++ii) 
        {
        A2(ii) = sqr(DR(ii))+sqr(DI(ii));
        //A2(ii) = fabs(DR(ii));
        if(A2(ii) > maxj) 
            {
            maxj = A2(ii);
            n = ii;
            }
        }

    //if which > 0, find next largest norm2, etc.
    for(int j = 1; j <= which; ++j)
        {
        Real nextmax = -1;
        for(int ii = 1; ii <= L; ++ii)
            {
            if(maxj > A2(ii) && A2(ii) > nextmax)
                {
                nextmax = A2(ii);
                n = ii;
                }
            }
        maxj = nextmax;
        }
    return (n-1);
    }

template <class BigMatrixT, class Tensor> 
Real
davidson(const BigMatrixT& A, Tensor& phi,
         const Args& args)
    {
    std::vector<Tensor> v(1);
    v.front() = phi;
    std::vector<Real> eigs = davidson(A,v,args);
    phi = v.front();
    return eigs.front();
    }

template <class BigMatrixT, class Tensor> 
std::vector<Real>
davidson(const BigMatrixT& A, 
         std::vector<Tensor>& phi,
         const Args& args)
    {
    auto debug_level_ = args.getInt("DebugLevel",-1);
    Real Approx0 = 1E-12;
    std::vector<Complex> ceigs = complexDavidson(A,phi,args);
    std::vector<Real> eigs(ceigs.size());
    for(size_t j = 0; j < ceigs.size(); ++j)
        {
        eigs.at(j) = ceigs.at(j).real();
        if(debug_level_ > 2 && ceigs.at(j).imag() > Approx0*ceigs.at(j).real())
            {
            printfln("Warning: dropping imaginary part of eigs[%d] = (%.4E,%.4E).", 
                     j , ceigs.at(j).real(), ceigs.at(j).imag());
            }
        }
    return eigs;
    }



template <class BigMatrixT, class Tensor> 
std::vector<Complex>
complexDavidson(const BigMatrixT& A, 
                std::vector<Tensor>& phi,
                const Args& args)
    {
    auto maxiter_ = args.getInt("MaxIter",2);
    auto errgoal_ = args.getReal("ErrGoal",1E-4);
    auto debug_level_ = args.getInt("DebugLevel",-1);
    auto miniter_ = args.getInt("MinIter",1);

    Real Approx0 = 1E-12;

    size_t nget = phi.size();
    if(nget == 0)
        {
        Error("No initial vectors passed to davidson.");
        }

    for(size_t j = 0; j < nget; ++j)
        {
        auto nrm = norm(phi[j]);
        if(nrm == 0.0)
            Error("norm of 0 in davidson");
        phi[j] *= 1.0/nrm;
        }

    bool complex_diag = false;

    auto maxsize = A.size();
    auto actual_maxiter = std::min(maxiter_,maxsize-1);
    if(debug_level_ >= 2)
        {
        printfln("maxsize-1 = %d, maxiter = %d, actual_maxiter = %d",
                 (maxsize-1), maxiter_, actual_maxiter);
        }

    if(area(phi.front().inds()) != maxsize)
        {
        Print(area(phi.front().inds()));
        Print(A.size());
        Error("davidson: size of initial vector should match linear matrix size");
        }

    std::vector<Tensor> V(actual_maxiter+2),
                       AV(actual_maxiter+2);

    //Storage for Matrix that gets diagonalized 
    //set to NAN to ensure failure if we use uninitialized elements
    Mat MR(actual_maxiter+2,actual_maxiter+2,NAN),
        MI(actual_maxiter+2,actual_maxiter+2,NAN);

    //Mref holds current projection of A into V's
    auto MrefR = subMatrix(MR,1,1,1,1);
    auto MrefI = subMatrix(MI,1,1,1,1);

    //Get diagonal of A to use later
    //auto Adiag = A.diag();

    auto last_lambda = Cplx(1000,0);

    Real qnorm = NAN;

    V[0] = phi.front();
    A.product(V[0],AV[0]);

    auto z = (dag(V[0])*AV[0]).cplx();
    auto initEn = z.real();

    if(debug_level_ > 2)
        printfln("Initial Davidson energy = %.10f",initEn);

    size_t t = 0; //which eigenvector we are currently targeting
    Vec D,DI;
    Mat UR,UI;

    std::vector<Complex> eigs(nget,Complex(NAN,NAN));

    int iter = 0;
    for(int ii = 0; ii <= actual_maxiter; ++ii)
        {
        //Diagonalize dag(V)*A*V
        //and compute the residual q

        auto ni = ii+1; 
        auto& q = V.at(ni);
        auto& phi_t = phi.at(t);
        auto& lambda = eigs.at(t);

        //Step A (or I) of Davidson (1975)
        if(ii == 0)
            {
            lambda = Cplx(initEn,0.);
            std::fill(MrefR.begin(),MrefR.end(),lambda.real());
            std::fill(MrefI.begin(),MrefI.end(),0);
            //Calculate residual q
            q = V[0];
            q *= -lambda.real();
            q += AV[0]; 
            }
        else // ii != 0
            {
            //Diagonalize M
            int w = t; //w 'which' variable needed because
                       //eigs returned by ComplexEigenValues and
                       //GenEigenvalues aren't sorted
            if(complex_diag)
                {
                Error("Complex Davidson not currently supported");
                //diagHermitian(MrefR,MrefI,UR,UI,D);
                //DI.resize(D.size());
                //std::fill(DI.begin(),DI.end(),0);

                ////Compute corresponding eigenvector
                ////phi_t of A from the min evec of M
                ////(and start calculating residual q)

                //phi_t = (UR(1,1+w)*Complex_1+UI(1,1+w)*Complex_i)*V[0];
                //q   = (UR(1,1+w)*Complex_1+UI(1,1+w)*Complex_i)*AV[0];
                //for(int k = 1; k <= ii; ++k)
                //    {
                //    auto cfac = (UR(k+1,1+w)*Complex_1+UI(k+1,1+w)*Complex_i);
                //    phi_t += cfac*V[k];
                //    q   += cfac*AV[k];
                //    }
                }
            else
                {
                diagSymmetric(MrefR,UR,D);
                DI.resize(D.size());
                std::fill(DI.begin(),DI.end(),0);

                //Compute corresponding eigenvector
                //phi_t of A from the min evec of M
                //(and start calculating residual q)
                phi_t = UR(1,1+w)*V[0];
                q   = UR(1,1+w)*AV[0];
                for(int k = 1; k <= ii; ++k)
                    {
                    phi_t += UR(k+1,1+w)*V[k];
                    q   += UR(k+1,1+w)*AV[k];
                    }
                }

            //lambda is the w^th eigenvalue of M
            lambda = Complex(D(1+w),DI(1+w));

            //Step B of Davidson (1975)
            //Calculate residual q
            if(lambda.imag() <= Approx0)
                q += (-lambda.real())*phi_t;
            else
                q += (-lambda)*phi_t;

            //Fix sign
            if(UR(1,1+w) < 0)
                {
                phi_t *= -1;
                q *= -1;
                }

            if(debug_level_ >= 3)
                {
                println("complex_diag = ", complex_diag ? "true" : "false");
                print("D = ",D);
                printfln("lambda = %.10f",D(1));
                }

            }

        //Step C of Davidson (1975)
        //Check convergence
        qnorm = norm(q);

        bool converged = (qnorm < errgoal_ && abs(lambda-last_lambda) < errgoal_) 
                         || qnorm < std::max(Approx0,errgoal_ * 1E-3);

        last_lambda = lambda;

        if((qnorm < 1E-20) || (converged && ii >= miniter_) || (ii == actual_maxiter))
            {
            if(t < (nget-1) && ii < actual_maxiter) 
                {
                ++t;
                last_lambda = Complex(1000,0);
                }
            else
                {
                if(debug_level_ >= 3) //Explain why breaking out of Davidson loop early
                    {
                    if((qnorm < errgoal_ && abs(lambda-last_lambda) < errgoal_))
                        {
                        printfln("Exiting Davidson because errgoal=%.0E reached",errgoal_);
                        }
                    else
                    if(ii < miniter_ || qnorm < std::max(Approx0,errgoal_ * 1.0e-3))
                        {
                        printfln("Exiting Davidson because small residual=%.0E obtained",qnorm);
                        }
                    else
                    if(ii == actual_maxiter)
                        {
                        println("Exiting Davidson because ii == actual_maxiter");
                        }
                    }

                goto done;
                }
            }
        
        if(debug_level_ >= 2 || (ii == 0 && debug_level_ >= 1))
            {
            printf("I %d q %.0E E",iter,qnorm);
            for(size_t j = 0; j < eigs.size(); ++j)
                {
                if(std::isnan(eigs[j].real())) break;
                if(fabs(eigs[j].imag()) > Approx0)
                    printf(" (%.10f,%.10f)",eigs[j].real(),eigs[j].imag());
                else
                    printf(" %.10f",eigs[j].real());
                }
            println();
            }

        //Compute next trial vector by
        //first applying Davidson preconditioner
        //formula then orthogonalizing against
        //other vectors

        //Step D of Davidson (1975)
        //Apply Davidson preconditioner
        //TODO
        //if(Adiag)
        //    {
        //    //Function which applies the mapping
        //    // f(x,theta) = 1/(theta - x)
        //    auto precond = [theta=lambda.real()](Real val)
        //        {
        //        return (theta==val) ? 0 : 1./(theta-val);
        //        };
        //    auto cond(Adiag);
        //    cond.apply(precond);
        //    q /= cond;
        //    }

        //Step E and F of Davidson (1975)
        //Do Gram-Schmidt on d (Npass times)
        //to include it in the subbasis
        int Npass = 1;
        std::vector<Complex> Vq(ni);

        int count = 0;
        for(int pass = 1; pass <= Npass; ++pass)
            {
            ++count;
            for(int k = 0; k < ni; ++k)
                {
                Vq[k] = (dag(V[k])*q).cplx();
                }

            for(int k = 0; k < ni; ++k)
                {
                q += (-Vq[k].real())*V[k];
                if(Vq[k].imag() != 0)
                    {
                    q += (-Vq[k].imag()*Complex_i)*V[k];
                    }
                }

            auto qnrm = norm(q);

            if(qnrm < 1E-10)
                {
                //Orthogonalization failure,
                //try randomizing
                if(debug_level_ >= 2)
                    println("Vector not independent, randomizing");
                q = V.at(ni-1);
                q = randomize(q);

                if(ni >= maxsize)
                    {
                    //Not be possible to orthogonalize if
                    //max size of q (vecSize after randomize)
                    //is size of current basis
                    if(debug_level_ >= 3)
                        println("Breaking out of Davidson: max Hilbert space size reached");
                    goto done;
                    }

                if(count > Npass * 3)
                    {
                    // Maybe the size of the matrix is only 1?
                    if(debug_level_ >= 3)
                        println("Breaking out of Davidson: count too big");
                    goto done;
                    }

                qnrm = norm(q);
                --pass;
                }

            q *= 1./qnrm;
            }

        if(debug_level_ >= 3)
            {
            if(fabs(norm(q)-1.0) > 1E-10)
                {
                Print(norm(q));
                Error("q not normalized after Gram Schmidt.");
                }
            }


        //Step G of Davidson (1975)
        //Expand AV and M
        //for next step
        A.product(V[ni],AV[ni]);

        //Step H of Davidson (1975)
        //Add new row and column to M
        MrefR = subMatrix(MR,1,ni+1,1,ni+1);
        MrefI = subMatrix(MI,1,ni+1,1,ni+1);
        Vec newColR(ni+1),
            newColI(ni+1);
        for(int k = 0; k <= ni; ++k)
            {
            z = (dag(V.at(k))*AV.at(ni)).cplx();
            newColR(k+1) = z.real();
            newColI(k+1) = z.imag();
            }
        column(MrefR,ni+1) &= newColR;
        column(MrefI,ni+1) &= newColI;
        row(MrefR,ni+1) &= newColR;
        row(MrefI,ni+1) &= newColI;
        row(MrefI,ni+1) *= -1;

        if(!complex_diag && norm(newColI) > errgoal_)
            {
            complex_diag = true;
            }

        ++iter;

        } //for(ii)

    done:

    //Compute any remaining eigenvalues and eigenvectors requested
    //(zero indexed) value of t indicates how many have been "targeted" so far
    for(size_t j = t+1; j < nget; ++j)
        {
        eigs.at(j) = Cplx(D(1+j),DI(1+j));

        auto& phi_j = phi.at(j);
        bool complex_evec = (norm(column(UI,1+t)) > Approx0);

        auto Nr = UR.Nrows();

        phi_j = UR(1,1+t)*V[0];
        for(int k = 1; k < Nr; ++k)
            {
            phi_j += UR(1+k,1+t)*V[k];
            }
        if(complex_evec)
            {
            phi_j += Complex_i*UI(1,1+t)*V[0];
            for(int k = 1; k < Nr; ++k)
                {
                phi_j += Complex_i*UI(1+k,1+t)*V[k];
                }
            }
        }

    if(debug_level_ >= 3)
        {
        //Check V's are orthonormal
        Mat Vo_final(iter+1,iter+1,NAN); 
        for(int r = 1; r <= iter+1; ++r)
        for(int c = r; c <= iter+1; ++c)
            {
            z = (dag(V[r-1])*V[c-1]).cplx();
            Vo_final(r,c) = abs(z);
            Vo_final(c,r) = Vo_final(r,c);
            }
        println("Vo_final = \n",Vo_final);
        }

    if(debug_level_ > 0)
        {
        printf("I %d q %.0E E",iter,qnorm);
        for(size_t j = 0; j < eigs.size(); ++j)
            {
            if(std::isnan(eigs[j].real())) break;
            if(fabs(eigs[j].imag()) > Approx0)
                printf(" (%.10f,%.10f)",eigs[j].real(),eigs[j].imag());
            else
                printf(" %.10f",eigs[j].real());
            }
        println();
        }

    return eigs;

    } //complexDavidson

} //namespace itensor


#endif
