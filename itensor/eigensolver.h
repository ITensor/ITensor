#ifndef __ITENSOR_EIGENSOLVER_H
#define __ITENSOR_EIGENSOLVER_H
#include "iqcombiner.h"

/* Notes on optimization:
 *
 * - Should be able to avoid re-forming the projected 
 *   A i.e. the matrix M. Instead just expand a single
 *   row and column using the new vector.
 *
 *
 */

class Davidson
    {
    public:

    //Constructors ---------------

    Davidson(int maxiter, Real errgoal = 1E-4, int numget = 1);

    //The solve method runs the actual Davidson algorithm

    template <class SparseT, class Tensor> Real 
    solve(const SparseT& A, Tensor& phi);

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

    void 
    read(std::istream& s);
    void 
    write(std::ostream& s) const;

    private:

    void
    combine(ITensor&, ITensor&, Index&) const;

    void
    combine(IQTensor&, IQTensor&, IQIndex&) const;

    int maxiter_;
    Real errgoal_;
    int numget_;

    }; //class Davidson

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
            return 1.0/(theta_-val+1E-33);
            }
    private:
        Real theta_;
    };

inline Davidson::
Davidson(int maxiter, Real errgoal, int numget)
    : maxiter_(maxiter),
      errgoal_(errgoal),
      numget_(numget)
    { }

template <class SparseT, class Tensor> 
inline Real Davidson::
solve(const SparseT& A, Tensor& phi)
    {
    typedef typename Tensor::IndexT 
    IndexT;

    const int maxsize = A.size();
    const int actual_maxiter = min(maxiter_,maxsize);
    Real lambda = 1E30, last_lambda = lambda;

    // p.m() is the same as Davidson's M
    // i.e. number of states in our current basis
    IndexT p = IQIndex("P0",Index("p0",1),QN());
    Tensor B = phi;
    B.addindex1(p);

    for(int iter = 1; iter <= actual_maxiter; ++iter)
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
            q *= conj(Tensor(p(1))); 
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
            Tensor alpha = conj(U) * Tensor(mid(mink));
            
            phi = alpha;
            phi *= B;

            //Calculate residual q
            q = phi;
            q *= -lambda;
            q += (AB * alpha);
            }
        std::cout << boost::format("At iter %d, lambda = %.10f")%iter%lambda << std::endl;

        //Check convergence (i.e. whether ||q|| is small)
        Real qnorm = q.norm();
        std::cout << boost::format("q.norm() = %.3E")%qnorm << std::endl;
        //std::cout << boost::format("fabs(lambda-last_lambda) = %.3E")%fabs(lambda-last_lambda) << std::endl;
        if( (qnorm < errgoal_ && fabs(lambda-last_lambda) < errgoal_) 
            || qnorm < 1E-12 )
            {
            std::cout << boost::format("Davidson: %d iterations, energy = %.10f")%iter%lambda << std::endl;
            return lambda;
            }

        //Apply Davidson preconditioner
        Tensor xi(q);
        {
        Tensor Ad(q);
        A.diag(Ad);
        DavidsonPrecond dp(lambda);
        Ad.mapElems(dp);
        xi /= Ad;
        }

        //Do Gram-Schmidt on xi
        //before including it in the subbasis
        Tensor d(xi);
        d *= conj(B); //m^2 d^2 p 
        d *= B;       //m^2 d^2 p
        d *= -1;
        d += xi;
        d *= 1.0/d.norm();

        //Combine d into B
        combine(d,B,p);

        last_lambda = lambda;

        } //for(iter)

    return lambda;

    } //Davidson::solve

inline void Davidson::
combine(ITensor& d, ITensor& B, Index& p) const
    {
    //Expand B's p-index
    Index oldp = p;
    p = Index(nameint("p",oldp.m()),oldp.m()+1);
    B.expandIndex(oldp,p,0);

    //Stick new p index onto d
    d *= ITensor(p(p.m()));

    //Combine them by adding
    B += d;
    }

inline void Davidson::
combine(IQTensor& d, IQTensor& B, IQIndex& P) const
    {
    if(P.nindex() != 1)
        Error("Basis IQIndex P should have a single block.");

    //Create a new version of P with a range expanded by 1
    Index oldp = P.index(1);
    Index newp = Index(nameint("p",oldp.m()),oldp.m()+1);
    IQIndex newP = IQIndex(nameint("P",oldp.m()),newp,QN());

    //Create a new IQTensor with expanded P IQIndex
    std::vector<IQIndex> iqinds;
    iqinds.reserve(B.r());
    Foreach(const IQIndex& I, B.iqinds())
        {
        if(I == P)
            iqinds.push_back(newP);
        else
            iqinds.push_back(I);
        }

    //Expand blocks of B and insert into newB
    IQTensor newB(iqinds);
    Foreach(ITensor t, B.itensors())
        {
        t.expandIndex(oldp,newp,0);
        newB.insert(t);
        }

    //Combine with d
    d *= IQTensor(newP(newP.m()));
    newB += d;

    B = newB;
    P = newP;
    }

#endif
