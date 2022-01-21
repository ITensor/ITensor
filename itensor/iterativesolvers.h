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
#ifndef __ITENSOR_ITERATIVESOLVERS_H
#define __ITENSOR_ITERATIVESOLVERS_H
#include "itensor/util/iterate.h"
#include "itensor/itensor.h"
#include "itensor/tensor/algs.h"


namespace itensor {

//
// Use the Davidson algorithm to find the 
// eigenvector of the Hermitian matrix A with minimal eigenvalue.
// (BigMatrixT objects must implement the methods product, size and diag.)
// Returns the minimal eigenvalue lambda such that
// A phi = lambda phi.
//
template <class BigMatrixT>
Real 
davidson(BigMatrixT const& A, 
         ITensor& phi,
         Args const& args = Args::global());

//
// Use Davidson to find the N eigenvectors with smallest 
// eigenvalues of the Hermitian matrix A, given a vector of N 
// initial guesses (zero indexed).
// (BigMatrixT objects must implement the methods product, size and diag.)
// Returns a vector of the N smallest eigenvalues corresponding
// to the set of eigenvectors phi.
//
template <class BigMatrixT>
std::vector<Real>
davidson(BigMatrixT const& A, 
         std::vector<ITensor>& phi,
         Args const& args = Args::global());

//
// Use GMRES to iteratively solve A x = b for x.
// (BigMatrixT objects must implement the methods product and size.)
// Initial guess x is overwritten with the output.
//
template<typename BigMatrixT, typename BigVectorT>
void
gmres(BigMatrixT const& A,
      BigVectorT const& b,
      BigVectorT& x,
      Args const& args = Args::global());

//
// Use the Krylov subspace method to compute
// phi' = exp(t*A)*phi
// for the BigMatrixT(e.g. localmpo or localmposet) A,
// and the cplx constant t.
// It does not compute the matrix exponential explicitly
// but instead compute the action of the exponential matrix on the vector.
// After finished, phi -> phi'.
template <typename BigMatrixT, typename ElT>
void
applyExp(BigMatrixT const& A,
         ITensor& phi,
         ElT t,
         Args const& args = Args::global());

//
//
// Implementations
//
//


template <class BigMatrixT>
Real
davidson(BigMatrixT const& A, 
         ITensor& phi,
         Args const& args)
    {
    auto v = std::vector<ITensor>(1);
    v.front() = phi;
    auto eigs = davidson(A,v,args);
    phi = v.front();
    return eigs.front();
    }

template <class BigMatrixT>
std::vector<Real>
davidson(BigMatrixT const& A, 
         std::vector<ITensor>& phi,
         Args const& args)
    {
    auto maxiter_ = args.getSizeT("MaxIter",2);
    auto errgoal_ = args.getReal("ErrGoal",1E-14);
    auto debug_level_ = args.getInt("DebugLevel",-1);
    auto miniter_ = args.getSizeT("MinIter",1);

    Real Approx0 = 1E-12;

    auto nget = phi.size();
    if(nget == 0) Error("No initial vectors passed to davidson.");
    for(auto j : range(nget))
        {
        auto nrm = norm(phi[j]);
        while(nrm == 0.0) 
            {
            phi[j].randomize();
            nrm = norm(phi[j]);
            }
        phi[j] *= 1./nrm;
        }

    size_t maxsize = A.size();
    size_t actual_maxiter = std::min(maxiter_,size_t(maxsize-1));
    if(debug_level_ >= 2)
        {
        printfln("maxsize-1 = %d, maxiter = %d, actual_maxiter = %d",
                 (maxsize-1), maxiter_, actual_maxiter);
        }

    if(dim(inds(phi.front())) != maxsize)
        {
        println("dim(inds(phi.front())) = ",dim(inds(phi.front())));
        println("A.size() = ",A.size());
        Error("davidson: size of initial vector should match linear matrix size");
        }

    auto V = std::vector<ITensor>(actual_maxiter+2);
    auto AV = std::vector<ITensor>(actual_maxiter+2);

    //Storage for Matrix that gets diagonalized 
    //set to NAN to ensure failure if we use uninitialized elements
    auto M = CMatrix(actual_maxiter+2,actual_maxiter+2);
    for(auto& el : M) el = Cplx(NAN,NAN);

    auto NC = CVector(actual_maxiter+2);

    //Mref holds current projection of A into V's
    auto Mref = subMatrix(M,0,1,0,1);

    //Get diagonal of A to use later
    //auto Adiag = A.diag();

    Real qnorm = NAN;

    Vector D;
    CMatrix U;

    Real last_lambda = 1000.;
    auto eigs = std::vector<Real>(nget,NAN);

    V[0] = phi.front();
TIMER_START(31);
    A.product(V[0],AV[0]);
TIMER_STOP(31);

    auto initEn = real(eltC((dag(V[0])*AV[0])));

    if(debug_level_ > 2)
        printfln("Initial Davidson energy = %.10f",initEn);

    auto t = size_t(0); //which eigenvector we are currently targeting

    auto iter = size_t(0);
    for(auto ii : range(actual_maxiter+1))
        {
        //Diagonalize dag(V)*A*V
        //and compute the residual q

        auto ni = ii+1; 
        auto& q = V[ni];
        auto& phi_t = phi.at(t);
        auto& lambda = eigs.at(t);

        //Step A (or I) of Davidson (1975)
        if(ii == 0)
            {
            lambda = initEn;
            stdx::fill(Mref,lambda);
            //Calculate residual q

            q = AV[0] - lambda*V[0];
            }
        else // ii != 0
            {
            Mref *= -1;
            if(debug_level_ > 3)
                {
                println("Mref = \n",Mref);
                }
            diagHermitian(Mref,U,D);
            Mref *= -1;
            D *= -1;
            lambda = D(t);
            phi_t = U(0,t)*V[0];
            q     = U(0,t)*AV[0];

            for(auto k : range1(ii))
                {
                phi_t += U(k,t)*V[k];
                q     += U(k,t)*AV[k];
                }

            //Step B of Davidson (1975)
            //Calculate residual q
            q += (-lambda)*phi_t;

            //Fix sign
            if(U(0,t).real() < 0)
                {
                phi_t *= -1;
                q *= -1;
                }
            if(debug_level_ >= 3)
                {
                println("D = ",D);
                printfln("lambda = %.10f",lambda);
                }
            //printfln("ii=%d, full q = \n%f",ii,q);
            }

        //Step C of Davidson (1975)
        //Check convergence
        qnorm = norm(q);

        bool converged = (qnorm < errgoal_ && std::abs(lambda-last_lambda) < errgoal_) 
                         || qnorm < std::max(Approx0,errgoal_ * 1E-3);

        last_lambda = lambda;

        if((qnorm < 1E-20) || (converged && ii >= miniter_) || (ii == actual_maxiter))
            {
            if(t < (nget-1) && ii < actual_maxiter) 
                {
                ++t;
                last_lambda = 1000.;
                }
            else
                {
                if(debug_level_ >= 3) //Explain why breaking out of Davidson loop early
                    {
                    if((qnorm < errgoal_ && std::fabs(lambda-last_lambda) < errgoal_))
                        printfln("Exiting Davidson because errgoal=%.0E reached",errgoal_);
                    else if(ii < miniter_ || qnorm < std::max(Approx0,errgoal_ * 1.0e-3))
                        printfln("Exiting Davidson because small residual=%.0E obtained",qnorm);
                    else if(ii == actual_maxiter)
                        println("Exiting Davidson because ii == actual_maxiter");
                    }

                goto done;
                }
            }
        
        if(debug_level_ >= 2 || (ii == 0 && debug_level_ >= 1))
            {
            printf("I %d q %.0E E",iter,qnorm);
            for(auto eig : eigs)
                {
                if(std::isnan(eig)) break;
                printf(" %.10f",eig);
                }
            println();
            }

        //Compute next trial vector by
        //first applying Davidson preconditioner
        //formula then orthogonalizing against
        //other vectors

        //Step D of Davidson (1975)
        //Apply Davidson preconditioner

        //
        //TODO add preconditioner step (may require
        //non-contracting product to do efficiently)
        //
        //if(Adiag)
        //    {
        //    //Function which applies the mapping
        //    // f(x,theta) = 1/(theta - x)
        //    auto precond = [theta=lambda.real()](Real val)
        //        {
        //        return (theta==val) ? 0 : 1./(theta-val);
        //        };
        //    auto cond= Adiag;
        //    cond.apply(precond);
        //    q /= cond;
        //    }

        //Step E and F of Davidson (1975)
        //Do Gram-Schmidt on d (Npass times)
        //to include it in the subbasis
        int Npass = 1;
        auto Vq = std::vector<Cplx>(ni);
        int pass = 1;
        int tot_pass = 0;
        while(pass <= Npass)
            {
            if(debug_level_ >= 3) println("Doing orthog pass");
            ++tot_pass;
            for(auto k : range(ni))
                {
                Vq[k] = eltC(dag(V[k])*q);
                //printfln("pass=%d Vq[%d] = %s",pass,k,Vq[k]);
                }
            for(auto k : range(ni))
                {
                q += (-Vq[k])*V[k];
                }
            auto qnrm = norm(q);
            //printfln("pass=%d qnrm=%s",pass,qnrm);
            if(qnrm < 1E-10)
                {
                //Orthogonalization failure,
                //try randomizing
                if(debug_level_ >= 2) println("Vector not independent, randomizing");
                q = V[ni-1];
                q.randomize();
                qnrm = norm(q);
                //Do another orthog pass
                --pass;
                if(debug_level_ >= 3) printfln("Now pass = %d",pass);

                if(ni >= maxsize)
                    {
                    //Not be possible to orthogonalize if
                    //max size of q (vecSize after randomize)
                    //is size of current basis
                    if(debug_level_ >= 3)
                        println("Breaking out of Davidson: max Hilbert space size reached");
                    goto done;
                    }

                if(tot_pass > Npass * 3)
                    {
                    // Maybe the size of the matrix is only 1?
                    if(debug_level_ >= 3)
                        println("Breaking out of Davidson: orthog step too big");
                    goto done;
                    }
                }
            q *= 1./qnrm;
            //q.scaleTo(1.);
            ++pass;
            }
        if(debug_level_ >= 3) println("Done with orthog step, tot_pass=",tot_pass);

        //Check V's are orthonormal
        //Mat Vo(ni+1,ni+1,NAN); 
        //for(int r = 1; r <= ni+1; ++r)
        //for(int c = r; c <= ni+1; ++c)
        //    {
        //    z = eltC(dag(V[r-1])*V[c-1]);
        //    Vo(r,c) = abs(z);
        //    Vo(c,r) = Vo(r,c);
        //    }
        //println("Vo = \n",Vo);

        if(debug_level_ >= 3)
            {
            if(std::fabs(norm(q)-1.0) > 1E-10)
                {
                println("norm(q) = ",norm(q));
                Error("q not normalized after Gram Schmidt.");
                }
            }

        //Step G of Davidson (1975)
        //Expand AV and M
        //for next step
TIMER_START(31);
        A.product(V[ni],AV[ni]);
TIMER_STOP(31);

        //Step H of Davidson (1975)
        //Add new row and column to M
        Mref = subMatrix(M,0,ni+1,0,ni+1);
        auto newCol = subVector(NC,0,1+ni);
        for(auto k : range(ni+1))
            {
            newCol(k) = eltC(dag(V.at(k))*AV.at(ni));
            }
        column(Mref,ni) &= newCol;
        row(Mref,ni) &= conj(newCol);

        ++iter;

        } //for(ii)

    done:

    //TODO: put this back?
    //for(auto& T : phi)
    //    {
    //    if(T.scale().logNum() > 2) T.scaleTo(1.);
    //    }

    //Compute any remaining eigenvalues and eigenvectors requested
    //(zero indexed) value of t indicates how many have been "targeted" so far
    if(debug_level_ >= 2 && t+1 < nget) printfln("Max iter. reached, computing remaining %d evecs",nget-t-1);
    for(auto j : range(t+1,nget))
        {
        eigs.at(j) = D(j);
        auto& phi_j = phi.at(j);
        auto Nr = size_t(nrows(U));
        phi_j = U(0,j)*V[0];
        for(auto k : range1(std::min(V.size(),Nr)-1))
            {
            phi_j += U(k,j)*V[k];
            }
        }

    if(debug_level_ >= 4)
        {
        //Check V's are orthonormal
        auto Vo_final = CMatrix(iter+1,iter+1);
        for(auto r : range(iter+1))
        for(auto c : range(r,iter+1))
            {
            auto z = eltC(dag(V[r])*V[c]);
            Vo_final(r,c) = std::abs(z);
            Vo_final(c,r) = Vo_final(r,c);
            }
        println("Vo_final = \n",Vo_final);
        }

    if(debug_level_ > 0)
        {
        printf("I %d q %.0E E",iter,qnorm);
        for(auto eig : eigs)
            {
            if(std::isnan(eig)) break;
            printf(" %.10f",eig);
            }
        println();
        }

    return eigs;
    }

namespace gmres_details {

template<class Matrix, class T, class BigVectorT>
void
update(BigVectorT &x, int const k, Matrix const& h, std::vector<T>& s, std::vector<BigVectorT> const& v)
    {
    std::vector<T> y(s);

    // Backsolve:
    for (int i = k; i >= 0; i--)
        {
        y[i] /= h(i,i);
        for (int j = i - 1; j >= 0; j--)
            y[j] -= h(j,i) * y[i];
        }

    for (int j = 0; j <= k; j++)
        x += y[j] * v[j];
    }

template<typename T>
void
generatePlaneRotation(T const& dx, T const& dy, T& cs, T& sn)
    {
    if(dy == 0.0)
        {
        cs = 1.0;
        sn = 0.0;
        }
    else if(std::abs(dy) > std::abs(dx))
        {
        auto temp = dx / dy;
        sn = 1.0 / std::sqrt( 1.0 + temp*temp );
        cs = temp * sn;
        }
    else
        {
        auto temp = dy / dx;
        cs = 1.0 / std::sqrt( 1.0 + temp*temp );
        sn = temp * cs;
        }
    }

void inline
applyPlaneRotation(Real& dx, Real& dy, Real const& cs, Real const& sn)
    {
    auto temp =  cs * dx + sn * dy;
    dy = -sn * dx + cs * dy;
    dx = temp;
    }

void inline
applyPlaneRotation(Cplx& dx, Cplx& dy, Cplx const& cs, Cplx const& sn)
    {
    auto temp =  std::conj(cs) * dx + std::conj(sn) * dy;
    dy = -sn * dx + cs * dy;
    dx = temp;
    }

template<typename BigVectorT>
void
dot(BigVectorT const& A, BigVectorT const& B, Real& res)
    {
    res = elt(dag(A)*B);
    }

template<typename BigVectorT>
void
dot(BigVectorT const& A, BigVectorT const& B, Cplx& res)
    {
    res = eltC(dag(A)*B);
    }

}//namespace gmres_details

template<typename T, typename BigMatrixT, typename BigVectorT>
void
gmresImpl(BigMatrixT const& A,
          BigVectorT const& b,
          BigVectorT& x,
          BigVectorT& Ax,
          Args const& args)
    {
    auto n = A.size();
    auto max_iter = args.getInt("MaxIter",n);
    auto m = args.getInt("RestartIter",max_iter);
    auto tol = args.getReal("ErrGoal",1E-14);
    auto debug_level_ = args.getInt("DebugLevel",-1);

    auto H = Mat<T>(m+1,m+1);

    int i;
    int j = 1;
    int k;

    std::vector<T> s(m+1);
    std::vector<T> cs(m+1);
    std::vector<T> sn(m+1);
    BigVectorT w = x;

    auto normb = norm(b);

    auto r = b - Ax;
    auto beta = norm(r);

    if(normb == 0.0)
        normb = 1.0;

    auto resid = norm(r)/normb;
    if(resid <= tol)
        {
        tol = resid;
        max_iter = 0;
        }

    std::vector<BigVectorT> v(m+1);

    while(j <= max_iter)
        {

        v[0] = r/beta;
        //v[0].scaleTo(1.0);

        std::fill(s.begin(), s.end(), 0.0);
        s[0] = beta; 

        for(i = 0; i < m && j <= max_iter; i++, j++)
            {
            BigVectorT w = x;
            A.product(v[i],w);

            // Begin Arnoldi iteration
            // TODO: turn into a function?
            for(k = 0; k<=i; ++k)
                {
                gmres_details::dot(w, v[k], H(k,i));
                w -= H(k,i)*v[k];
                }
            auto normw = norm(w);

            if(debug_level_ > 0)
                println("norm(w) = ", normw);

            H(i+1,i) = normw;
            if(normw != 0)
                {
                v[i+1] = w/H(i+1,i);
                //v[i+1].scaleTo(1.0);
                }
            //else
            //    {
            //    // Maybe this should be a warning?
            //    // Also, maybe check if it is very close to zero?
            //    // GMRES generally is converged at this point anyway
            //    println("Warning: norm of new Krylov vector is zero.");
            //    }

            for(k = 0; k<i; ++k)
                gmres_details::applyPlaneRotation(H(k,i), H(k+1,i), cs[k], sn[k]);

            gmres_details::generatePlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
            gmres_details::applyPlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
            gmres_details::applyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);

            resid = std::abs(s[i+1])/normb;

            if(resid < tol)
                {
                gmres_details::update(x, i, H, s, v);
                return;
                }

            } // end for loop

            gmres_details::update(x, i-1, H, s, v);
            A.product(x, Ax);
            r = b - Ax;
            beta = norm(r);
            resid = beta/normb;
            if(resid < tol)
                return;

        } // end while loop

    }


template<typename BigMatrixT, typename BigVectorT>
void
gmres(BigMatrixT const& A,
      BigVectorT const& b,
      BigVectorT& x,
      Args const& args)
    {
    auto debug_level_ = args.getInt("DebugLevel",-1);

    // Precompute Ax to figure out whether A or x is
    // complex, maybe there is a cleaner code design
    // to avoid this?
    // Otherwise we would need to require that BigMatrixT
    // has a function isComplex()
    BigVectorT Ax = x;
    A.product(x, Ax); 
    if(isComplex(b) || isComplex(Ax))
        {
        if(debug_level_ > 0)
            println("Calling complex version of gmresImpl()");
        gmresImpl<Cplx>(A,b,x,Ax,args);
        }
    else
        {
        if(debug_level_ > 0)
            println("Calling real version of gmresImpl()");
        gmresImpl<Real>(A,b,x,Ax,args);
        }
    }

int inline
findEig(Vector const& vr, Vector const& vi, std::string whichEig)
    {
    int n = -1;
    Real foundval = NAN;
    for(size_t i = 0; i < vr.size(); i++)
      {
      if(whichEig == "LargestMagnitude")
        {
        auto ival = abs(Complex(vr(i),vi(i)));
        if(i == 0)
          {
          foundval = ival;
          n = 0;
          }
        else if(ival > foundval)
          {
          foundval = ival;
          n = i;
          }
        }
      else if(whichEig == "SmallestReal")
        {
        auto ival = vr(i);
        if(i == 0)
          {
          foundval = ival;
          n = 0;
          }
        else if(ival < foundval)
          {
          foundval = ival;
          n = i;
          }
        }
      else
        {
        error("Unsupported eigenvalue target, currently only support: LargestMagnitude, SmallestReal");        
        }
      }
    return n;
    }
  
template <class BigMatrixT>
std::vector<Complex>
arnoldi(const BigMatrixT& A,
        std::vector<ITensor>& phi,
        Args const& args)
    {
    int maxiter_ = args.getInt("MaxIter",10);
    int maxrestart_ = args.getInt("MaxRestart",0);
    std::string whicheig_ = args.getString("WhichEig","LargestMagnitude");
    const Real errgoal_ = args.getReal("ErrGoal",1E-6);
    const int debug_level_ = args.getInt("DebugLevel",-1);

    if(maxiter_ < 1) maxiter_ = 1;
    if(maxrestart_ < 0) maxrestart_ = 0;

    const Real Approx0 = 1E-12;
    const int Npass = args.getInt("Npass",2); // number of Gram-Schmidt passes

    const size_t nget = phi.size();
    if(nget == 0) Error("No initial vectors passed to arnoldi.");

    //if(nget > 1) Error("arnoldi currently only supports nget == 1");

    for(size_t j = 0; j < nget; ++j)
        {
        const Real nrm = norm(phi[j]);
        if(nrm == 0.0)
            Error("norm of 0 in arnoldi");
        phi[j] *= 1.0/nrm;
        }

    std::vector<Complex> eigs(nget);

    const int maxsize = A.size();

    if(phi.size() > size_t(maxsize))
        Error("arnoldi: requested more eigenvectors (phi.size()) than size of matrix (A.size())");

    if(maxsize == 1)
        {
        if(norm(phi.front()) == 0) phi.front().randomize();
        phi.front() /= norm(phi.front());
        ITensor Aphi(phi.front());
        A.product(phi.front(),Aphi);
        //eigs.front() = BraKet(Aphi,phi.front());
        gmres_details::dot(Aphi,phi.front(),eigs.front());
        return eigs;
        }

    auto actual_maxiter = std::min(maxiter_,maxsize-1);
    if(debug_level_ >= 2)
        {
        printfln("maxsize-1 = %d, maxiter = %d, actual_maxiter = %d",
                 (maxsize-1),     maxiter_ ,    actual_maxiter );
        }

    if(dim(phi.front().inds()) != size_t(maxsize))
        {
        Error("arnoldi: size of initial vector should match linear matrix size");
        }

    //Storage for Matrix that gets diagonalized 
    Matrix HR(actual_maxiter+2,actual_maxiter+2),
           HI(actual_maxiter+2,actual_maxiter+2);
    //HR = 0;
    //HI = 0;
    for(auto& el : HR) el = 0;
    for(auto& el : HI) el = 0;

    std::vector<ITensor> V(actual_maxiter+2);

    for(size_t w = 0; w < nget; ++w)
    {

    for(int r = 0; r <= maxrestart_; ++r)
        {
        Real err = 1000;
        Matrix YR,YI;
        int n = 0; //which column of Y holds the w^th eigenvector
        int niter = 0;

        //Mref holds current projection of A into V's
        MatrixRef HrefR(subMatrix(HR,0,1,0,1)),
                  HrefI(subMatrix(HI,0,1,0,1));

        V.at(0) = phi.at(w);

        for(int it = 0; it <= actual_maxiter; ++it)
            {
            const int j = it;
            A.product(V.at(j),V.at(j+1)); // V[j+1] = A*V[j]
            // "Deflate" previous eigenpairs:
            for(size_t o = 0; o < w; ++o)
                {
                //V[j+1] += (-eigs.at(o)*phi[o]*BraKet(phi[o],V[j]));
                Complex overlap_;
                gmres_details::dot(phi[o],V[j],overlap_);
                V[j+1] += (-eigs.at(o)*phi[o]*overlap_);
                }

            //Do Gram-Schmidt orthogonalization Npass times
            //Build H matrix only on the first pass
            Real nh = NAN;
            for(int pass = 1; pass <= Npass; ++pass)
                {
                for(int i = 0; i <= j; ++i)
                    {
                    //Complex h = BraKet(V.at(i),V.at(j+1));
                    Complex h;
                    gmres_details::dot(V.at(i),V.at(j+1),h);
                    if(pass == 1)
                        {
                        HR(i,j) = h.real();
                        HI(i,j) = h.imag();
                        }
                    V.at(j+1) -= h*V.at(i);
                    }
                Real nrm = norm(V.at(j+1));
                if(pass == 1) nh = nrm;

                if(nrm != 0) V.at(j+1) /= nrm;
                else         V.at(j+1).randomize();
                }

            //for(int i1 = 0; i1 <= j+1; ++i1)
            //for(int i2 = 0; i2 <= j+1; ++i2)
            //    {
            //    auto olap = BraKet(V.at(i1),V.at(i2)).real();
            //    if(fabs(olap) > 1E-12)
            //        Cout << Format(" %.2E") % BraKet(V.at(i1),V.at(i2)).real();
            //    }
            //Cout << Endl;

            //Diagonalize projected form of A to
            //obtain the w^th eigenvalue and eigenvector
            Vector D(1+j),DI(1+j);

            //TODO: eigen only takes a Matrix of Complex, not
            //the real and imaginary parts seperately.
            //Change it so that we don't have to allocate this
            //Complex matrix
            auto Hnrows = nrows(HrefR);
            auto Hncols = ncols(HrefR);
            CMatrix Href(Hnrows,Hncols);
            for(size_t irows = 0; irows < Hnrows; irows++)
              for(size_t icols = 0; icols < Hncols; icols++)
                Href(irows,icols) = Complex(HrefR(irows,icols),HrefI(irows,icols));

            eigen(Href,YR,YI,D,DI);
            n = findEig(D,DI,whicheig_); //continue to target the largest eig 
                                        //since we have 'deflated' the previous ones
            eigs.at(w) = Complex(D(n),DI(n));

            HrefR = subMatrix(HR,0,j+2,0,j+2);
            HrefI = subMatrix(HI,0,j+2,0,j+2);

            HR(1+j,j) = nh;

            //Estimate error || (A-l_j*I)*p_j || = h_{j+1,j}*[last entry of Y_j]
            //See http://web.eecs.utk.edu/~dongarra/etemplates/node216.html
            assert(nrows(YR) == size_t(1+j));
            err = nh*abs(Complex(YR(j,n),YI(j,n)));
            assert(err >= 0);

            if(debug_level_ >= 1)
                {
                if(r == 0)
                    printf("I %d e %.0E E",(1+j),err);
                else
                    printf("R %d I %d e %.0E E",r,(1+j),err);

                for(size_t j = 0; j <= w; ++j)
                    {
                    if(fabs(eigs[j].real()) > 1E-6)
                        {
                        if(fabs(eigs[j].imag()) > Approx0)
                            printf(" (%.10f,%.10f)",eigs[j].real(),eigs[j].imag());
                        else
                            printf(" %.10f",eigs[j].real());
                        }
                    else
                        {
                        if(fabs(eigs[j].imag()) > Approx0)
                            printf(" (%.5E,%.5E)",eigs[j].real(),eigs[j].imag());
                        else
                            printf(" %.5E",eigs[j].real());
                        }
                    }
                println();
                }

            ++niter;

            if(err < errgoal_) break;

            } // for loop over j

        //Cout << Endl;
        //for(int i = 0; i < niter; ++i)
        //for(int j = 0; j < niter; ++j)
        //    Cout << Format("<V[%d]|V[%d]> = %.5E") % i % j % BraKet(V.at(i),V.at(j)) << Endl;
        //Cout << Endl;

        //Compute w^th eigenvector of A
        //Cout << Format("Computing eigenvector %d") % w << Endl;
        phi.at(w) = Complex(YR(0,n),YI(0,n))*V.at(0);
        for(int j = 1; j < niter; ++j)
            {
            phi.at(w) += Complex(YR(j,n),YI(j,n))*V.at(j);
            }

        //Print(YR.Column(1+n));
        //Print(YI.Column(1+n));

        const Real nrm = norm(phi.at(w));
        if(nrm != 0)
            phi.at(w) /= nrm;
        else
            phi.at(w).randomize();

        if(err < errgoal_) break;

        //otherwise restart using the phi.at(w) computed above

        } // for loop over r

    } // for loop over w

    return eigs;
    }

template <class BigMatrixT>
Complex
arnoldi(const BigMatrixT& A,
        ITensor& vec,
        Args const& args = Args::global())
    {
    std::vector<ITensor> phi(1,vec);
    Complex res = arnoldi(A,phi,args).front();
    vec = phi.front();
    return res;
    }

template<typename VecT>
void
assembleLanczosVectors(std::vector<ITensor> const& lanczos_vectors,
                       VecT const& linear_comb,
                       double norm, ITensor& phi)
    {
    assert(lanczos_vectors.size() == linear_comb.size());
    phi = norm*linear_comb(0)*lanczos_vectors[0];
    for(int i=1; i<(int)lanczos_vectors.size(); ++i)
        phi += norm*linear_comb(i)*lanczos_vectors[i];
    }

template<typename BigMatrixT, typename ElT>
void
applyExp(BigMatrixT const& H, ITensor& phi,
         ElT tau, Args const& args)
    {
    auto tol = args.getReal("ErrGoal",1E-10);
    auto max_iter = args.getInt("MaxIter",30);
    auto debug_level = args.getInt("DebugLevel",-1);
    auto beta_tol = args.getReal("NormCutoff",1e-7);

    // Initialize Lanczos vectors
    ITensor v1 = phi;
    ITensor v0;
    ITensor w;
    Real nrm = norm(v1);
    v1 /= nrm;
    std::vector<ITensor> lanczos_vectors({v1});
    Matrix bigTmat(max_iter + 2, max_iter + 2);
    std::fill(bigTmat.begin(), bigTmat.begin()+bigTmat.size(), 0.);

    auto nmatvec = 0;

    double beta = 0;
    for (int iter=0; iter < max_iter; ++iter)
        {
        int tmat_size=iter+1;
        // Matrix-vector multiplication
        if(debug_level >= 0)
            nmatvec++;
        H.product(v1, w);

        double avnorm = norm(w);
        double alpha = real(eltC(dag(w) * v1));
        bigTmat(iter, iter) = alpha;
        w -= alpha * v1;
        if (iter > 0)
            w -= beta * v0;
        v0 = v1;
        beta = norm(w);

        // check for Lanczos sequence exhaustion
        if (std::abs(beta) < beta_tol)
            {
            // Assemble the time evolved state
            auto tmat = subMatrix(bigTmat, 0,tmat_size, 0, tmat_size);
            auto tmat_exp = expMatrix(tmat, tau);
            auto linear_comb = column(tmat_exp, 0);
            assembleLanczosVectors(lanczos_vectors, linear_comb, nrm, phi);
            break;
            }

        // update next lanczos vector
        v1 = w;
        v1 /= beta;
        lanczos_vectors.push_back(v1);
        bigTmat(iter+1, iter) = beta;
        bigTmat(iter, iter+1) = beta;

        // Convergence check
        if (iter > 0)
            {
            // Prepare extended T-matrix for exponentiation
            int tmat_ext_size = tmat_size + 2;
            auto tmat_ext = Matrix(tmat_ext_size, tmat_ext_size);
            tmat_ext = subMatrix(bigTmat, 0, tmat_ext_size, 0, tmat_ext_size);

            tmat_ext(tmat_size-1, tmat_size) = 0.;
            tmat_ext(tmat_size+1, tmat_size) = 1.;

            // Exponentiate extended T-matrix
            auto tmat_ext_exp = expMatrix(tmat_ext, tau);

            double phi1 = std::abs( nrm*tmat_ext_exp(tmat_size, 0) );
            double phi2 = std::abs( nrm*tmat_ext_exp(tmat_size + 1, 0) * avnorm );
            double error;
            if (phi1 > 10*phi2) error = phi2;
            else if (phi1 > phi2) error = (phi1*phi2)/(phi1-phi2);
            else error = phi1;
            if(debug_level >= 1)
                println("Iteration: ", iter, ", Error: ", error);
            if ((error < tol) || (iter == max_iter-1))
                {
                if (iter == max_iter-1)
                    printf("warning: applyExp not converged in %d steps\n", max_iter);

                // Assemble the time evolved state
                auto linear_comb = Vec<ElT>(tmat_ext_size);
                linear_comb = column(tmat_ext_exp, 0);
                linear_comb = subVector(linear_comb, 0, tmat_ext_size-1);
                assembleLanczosVectors(lanczos_vectors, linear_comb, nrm, phi);
                if(debug_level >= 0)
                    printf("In applyExp, number of iterations: %d\n", iter);
                break;
                }
            }  // end convergence test

        }  // Lanczos iteratrions

    if(debug_level >= 0)
        println("In applyExp, number of matrix-vector multiplies: ", nmatvec);

    return;
    } // End applyExp

template<typename BigMatrixT>
void
applyExp(BigMatrixT const& H, ITensor& phi,
         int tau, Args const& args)
    {
    applyExp(H,phi,Real(tau),args);
    return;
    }

} //namespace itensor

#endif

