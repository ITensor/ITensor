#ifndef __ITENSOR_EIGENSOLVER_H
#define __ITENSOR_EIGENSOLVER_H
#include "iqcombiner.h"

class Davidson
    {
    public:

    //Constructors ---------------

    Davidson(int niter, Real errgoal = 1E-4, int numget = 1);

    //The solve method runs the actual Davidson algorithm

    template <class SparseT, class Tensor> Real 
    solve(const SparseT& A, Tensor& phi);

    //Accessor methods ------------

    Real errgoal() const { return errgoal_; }
    void errgoal(Real val) { errgoal_ = val; }

    int numGet() const { return numget_; }
    void numGet(int val) { numget_ = val; }

    void 
    read(std::istream& s);
    void 
    write(std::ostream& s) const;

    private:

    int niter_;
    Real errgoal_;
    int numget_;

    }; //class Davidson

inline Davidson::
Davidson(int niter, Real errgoal, int numget)
    : niter_(niter),
      errgoal_(errgoal),
      numget_(numget)
    { }

template <class SparseT, class Tensor> 
inline Real Davidson::
solve(const SparseT& A, Tensor& phi)
    {
    typedef typename Tensor::IndexT IndexT;
    typedef typename Tensor::CombinerT CombinerT;

    Real lambda = 1E30;

    // p is the same as Davidson's M
    // i.e. number of states in our 
    IndexT p("p0",1);
    const Index np("np",1);
    Tensor B = phi;
    B.addindex1(p);
    for(int iter = 1; iter <= 100*numget_; ++iter)
        {
        //std::cout << boost::format("Iteration %d -------------------")%iter << std::endl;
        Tensor AB;
        A.product(B,AB);

        //Compute q which is the difference
        //between the action of A on our
        //candidate eigenstate and this
        //state times its approxiate eigenvalue
        Tensor q;
        if(p.m() == 1)
            {
            lambda = Dot(B,AB);

            //Calculate residual q
            q = B;
            q *= -lambda;
            q += AB; 
            q *= Tensor(p,1); 
            }
        else // p.m() != 1
            {
            Tensor M = AB;
            M *= conj(primeind(B,p));

            //Diagonalize M
            Tensor D,U;
            IndexT mid;
            int mink=-1,maxk=-1;
            M.symmetricDiag11(p,D,U,mid,mink,maxk);

            //lambda is the minimum eigenvalue of M
            lambda = D(mid(mink));
            //alpha pick out the corresponding eigenvector
            Tensor alpha = U * Tensor(conj(mid)(mink));
            
            phi = alpha; 
            phi *= B;

            //Calculate residual q
            q = phi;
            q *= -lambda;
            q += (AB * alpha);
            }
        std::cout << boost::format("At iter %d, lambda = %.10f")%iter%lambda << std::endl;

        //Check convergence (i.e. whether ||q|| is small)
        std::cout << boost::format("q.norm() = %.3E")%q.norm() << std::endl;
        if(q.norm() < errgoal_)
            {
            std::cout << boost::format("Davidson: %d iterations, energy = %.10f")%iter%lambda << std::endl;
            return lambda;
            }

        //Apply Davidson preconditioner
        Tensor Ad(q);
        A.diag(Ad);
        const int size = Ad.vecSize();
        Vector qv(size),dv(size);
        Ad.assignToVec(dv);
        q.assignToVec(qv);
        for(int j = 1; j <= size; ++j)
            qv(j) /= -(dv(j)-lambda+1E-33);
        Tensor xi(q);
        xi.assignFromVec(qv);


        //Perform Gram-Schmidt on xi
        //before including it in the subbasis
        Tensor d(xi);
        d *= conj(B); //m^2 d^2 p 
        d *= B;       //m^2 d^2 p
        d *= -1;
        d += xi;
        d *= 1.0/d.norm();

        //Combine d into B
        Index oldp = p;
        p = Index(nameint("p",iter),oldp.m()+1);
        Tensor newB;
        B.expandIndex(oldp,p,0,newB);

        d *= Tensor(p(p.m()));

        newB += d;

        B = newB;

        } //for(iter)

    return lambda;

    } //Davidson::solve


#endif
