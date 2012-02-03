#ifndef __ITENSOR_EIGENSOLVER_H
#define __ITENSOR_EIGENSOLVER_H
#include "iqcombiner.h"

/* Notes on optimization:
 *
 * - Important not to re-form AB every time since
 *   multiplication by A is the most expensive
 *   part. Turn AB into a vector of Tensors?
 *   Or try multiplying d times A (Ad), then
 *   adding that into AB.
 *
 * - Should be able to avoid re-forming the projected 
 *   A i.e. the matrix M. Instead just expand a single
 *   row and column using the new vector.
 *
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

    void
    combine(ITensor&, ITensor&, const Index&, const Index&) const;
    void
    combine(IQTensor&, IQTensor&, const IQIndex&, const IQIndex&) const;

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
    typedef typename Tensor::IndexT 
    IndexT;

    phi *= 1.0/phi.norm();

    const int maxsize = A.size();
    const int actual_maxiter = min(maxiter_,maxsize);
    Real lambda = 1E30, last_lambda = lambda;
    Real qnorm = 1E30;

    // p.m() is the same as Davidson's M
    // i.e. number of states in our current basis
    IndexT P = IQIndex("P0",Index("p0",1),QN());
    Tensor B = phi;
    B.addindex1(P);

    Tensor AB;
    A.product(B,AB);

    for(int iter = 1; iter <= actual_maxiter; ++iter)
        {
        //Diagonalize conj(B)*A*B
        //and compute the residual q
        Tensor q;
        if(P.m() == 1)
            {
            //No need to diagonalize
            lambda = Dot(B,AB);

            //Calculate residual q
            q = B;
            q *= -lambda;
            q += AB; 
            q *= conj(Tensor(P(1))); 
            }
        else // p.m() != 1
            {
            Tensor M = AB;
            M *= conj(primeind(B,P));

            //Diagonalize M
            Tensor D,U;
            IndexT mid;
            int mink=-1,maxk=-1;
            M.symmetricDiag11(P,D,U,mid,mink,maxk);

            //lambda is the minimum eigenvalue of M
            lambda = D(mid(mink));
            //alpha picks out the corresponding eigenvector
            Tensor alpha = conj(U) * Tensor(mid(mink));
            
            //Set phi to the current best
            //eigenvector
            phi = alpha;
            phi *= B;

            //Calculate residual q
            q = phi;
            q *= -lambda;
            q += (AB * alpha);
            }

        //Check convergence
        qnorm = q.norm();
        if( (qnorm < errgoal_ && fabs(lambda-last_lambda) < errgoal_) 
            || qnorm < 1E-12 )
            {
            if(debug_level_ > 0)
                {
                std::cout << boost::format("Iter %d, lambda = %.10f")%iter%lambda << std::endl;
                std::cout << boost::format("qnorm = %.1E, lambda err = %.1E")%qnorm % fabs(lambda-last_lambda) << std::endl;
                }
            return lambda;
            }

        if(debug_level_ > 1 || (iter == 1 && debug_level_ > 0))
            {
            std::cout << boost::format("Iter %d, lambda = %.10f")%iter%lambda << std::endl;
            std::cout << boost::format("qnorm = %.1E, lambda err = %.1E")%qnorm % fabs(lambda-last_lambda) << std::endl;
            }

        //Apply Davidson preconditioner
        Tensor xi(q);
        {
        Tensor cond(q);
        A.diag(cond);
        DavidsonPrecond dp(lambda);
        cond.mapElems(dp);
        xi /= cond;
        }

        //Do Gram-Schmidt on xi
        //before including it in the subbasis
        Tensor d(xi);
        d *= conj(B); //m^2 d^2 p 
        d *= B;       //m^2 d^2 p
        d *= -1;
        d += xi;
        d *= 1.0/d.norm();

        Tensor Ad;
        A.product(d,Ad);

        IQIndex newP = IQIndex(nameint("P",P.m()),
                         Index(nameint("p",P.m()),P.m()+1),
                         QN());

        //Combine d into B
        combine(d,B,P,newP);

        //Combine Ad into AB
        combine(Ad,AB,P,newP);

        //B.scaleTo(1);

        P = newP;

        last_lambda = lambda;

        } //for(iter)

    if(debug_level_ > 0)
        {
        std::cout << boost::format("Iter %d, lambda = %.10f")%actual_maxiter%lambda << std::endl;
        std::cout << boost::format("qnorm = %.1E, lambda err = %.1E")%qnorm % fabs(lambda-last_lambda) << std::endl;
        }
    return lambda;

    } //Eigensolver::davidson

inline void Eigensolver::
combine(ITensor& d, ITensor& B, const Index& p, const Index& newp) const
    {
    //Expand B's p-index
    B.expandIndex(p,newp,0);

    //Stick new p index onto d
    d *= ITensor(newp(newp.m()));

    //Combine them by adding
    B += d;
    }

inline void Eigensolver::
combine(IQTensor& d, IQTensor& B, const IQIndex& P, const IQIndex& newP) const
    {
    if(P.nindex() != 1)
        Error("Basis IQIndex P should have a single block.");

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
        t.expandIndex(P.index(1),newP.index(1),0);
        newB.insert(t);
        }

    //Combine with d
    d *= IQTensor(newP(newP.m()));
    newB += d;

    B = newB;
    }

#endif
