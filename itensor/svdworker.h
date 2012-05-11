//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SVDWORKER_H
#define __ITENSOR_SVDWORKER_H
#include "iqcombiner.h"
#include "iqtsparse.h"

enum Direction { Fromright, Fromleft, Both, None };

template<class TensorA, class TensorB>
typename TensorA::IndexT
index_in_common(const TensorA& A, const TensorB& B, IndexType t)
    {
    typedef typename TensorA::IndexT
    IndexT;

    for(int j = 1; j <= A.r(); ++j)
        {
        const IndexT& I = A.index(j);
        if(I.type() == t && B.hasindex(I)) { return I; }
        }

    throw ITError("No common index found");
    return IndexT();
    }


//
// SVDWorker
//

class SVDWorker
    {
    public:

    //
    // Constructors
    //

    SVDWorker();

    SVDWorker(int N_);

    SVDWorker(int N_, Real cutoff, int minm, int maxm, 
              bool doRelCutoff, const LogNumber& refNorm);

    SVDWorker(std::istream& s) { read(s); }

    //
    // Site Canonical SVD 
    // (weights in both L and R, inverse weights in V)
    //

    template <class Tensor,class SparseT> 
    void 
    csvd(const Tensor& AA, Tensor& L, SparseT& V, Tensor& R)
        { 
        csvd<Tensor>(1,AA,L,V,R); 
        }

    template <class Tensor, class SparseT> 
    void 
    csvd(int b, const Tensor& AA, Tensor& L, SparseT& V, Tensor& R);

    //Object PH of type LocalOpT has to provide the deltaPhi method
    //to use the noise term feature (see localop.h)
    template <class Tensor, class SparseT, class LocalOpT>
    void 
    csvd(int b, const Tensor& AA, Tensor& L, SparseT& V, Tensor& R, 
         const LocalOpT& PH);


    //
    // Density Matrix Decomposition
    // (weights in B for dir == Fromleft, weights in A for dir == Fromright)
    //

    template <class Tensor> 
    void 
    denmatDecomp(const Tensor& AA, Tensor& A, Tensor& B, Direction dir)
        { 
        denmatDecomp<Tensor>(1,AA,A,B,dir); 
        }

    template <class Tensor> 
    void
    denmatDecomp(int b, const Tensor& AA, Tensor& A, Tensor& B, Direction dir);

    //Object PH of type LocalOpT has to provide the deltaRho method
    //to use the noise term feature (see localop.h)
    template <class Tensor, class LocalOpT>
    void
    denmatDecomp(int b, const Tensor& AA, Tensor& A, Tensor& B, Direction dir, 
                 const LocalOpT& PH);


    //
    // Singular Value Decomposition
    // (a.k.a. bond canonical SVD; weights in D)
    //

    template <class Tensor, class SparseT> 
    void 
    svd(const Tensor& AA, Tensor& U, SparseT& D, Tensor& V)
        { 
        svd<Tensor>(1,AA,U,D,V); 
        }

    template <class Tensor, class SparseT> 
    void 
    svd(int b, const Tensor& AA, Tensor& U, SparseT& D, Tensor& V);

    //Object PH of type LocalOpT has to provide the deltaPhi method
    //to use the noise term feature (see localop.h)
    template <class Tensor, class SparseT, class LocalOpT>
    void 
    svd(int b, Tensor AA, Tensor& U, SparseT& D, Tensor& V, 
        const LocalOpT& PH);


    //
    // Accessor Methods
    //

    int 
    NN() const { return N; }

    Real 
    cutoff() const { return cutoff_; }
    void 
    cutoff(Real val) { cutoff_ = val; }

    Real 
    truncerr(int b = 1) const { return truncerr_.at(b); }

    int 
    minm() const { return minm_; }
    void 
    minm(int val) { minm_ = val; }

    int 
    maxm() const { return maxm_; }
    void 
    maxm(int val) { maxm_ = val; }

    bool 
    useOrigM() const { return use_orig_m_; }
    void 
    useOrigM(bool val) { use_orig_m_ = val; }

    bool 
    showeigs() const { return showeigs_; }
    void 
    showeigs(bool val) { showeigs_ = val; }

    // If doRelCutoff_ is false,
    // refNorm_ defines an overall scale factored
    // out of the denmat before truncating.
    //
    // If doRelCutoff_ is true, refNorm_
    // is determined automatically.
    //
    bool 
    doRelCutoff() const { return doRelCutoff_; }
    void 
    doRelCutoff(bool val) { doRelCutoff_ = val; }

    // If absoluteCutoff_ == true, do a strict limit
    // on the eigenvalues of the density matrix, 
    // not trying for a truncation error
    bool 
    absoluteCutoff() const { return absoluteCutoff_; }
    void 
    absoluteCutoff(bool val) { absoluteCutoff_ = val; }

    LogNumber 
    refNorm() const { return refNorm_; }
    void 
    refNorm(const LogNumber& val) 
        { 
        if(val.sign() == 0) Error("zero refNorm");
        refNorm_ = val; 
        }

    const Vector& 
    eigsKept(int b = 1) const 
        { return eigsKept_.at(b); }
    int
    numEigsKept(int b = 1) const 
        { return eigsKept_.at(b).Length(); }

    int 
    maxEigsKept() const;

    Real 
    maxTruncerr() const;

    Real
    noise() const { return noise_; }
    void
    noise(Real val) { noise_ = val; }

    //
    // Other Methods
    //

    Real 
    diag_denmat(const ITensor& rho, Vector& D, Index& newmid, ITensor& U);
    Real 
    diag_denmat(const IQTensor& rho, Vector& D, IQIndex& newmid, IQTensor& U);

    Real 
    diag_denmat(const ITensor& rho, Vector& D, Index& newmid, 
                ITensor& C, ITensor& U);
    Real 
    diag_denmat(const IQTensor& rho, Vector& D, IQIndex& newmid, 
                IQTensor& C, IQTensor& U);

    Real 
    diag_denmat_complex(const IQTensor& rho, Vector& D, IQIndex& newmid, 
                        IQTensor& U);
    Real 
    diag_denmat_complex(const ITensor& rho, Vector& D, Index& newmid, 
                        ITensor& U);

    void 
    read(std::istream& s);
    void 
    write(std::ostream& s) const;


    void 
    svdRank2(const ITensor& A, const Index& ui, const Index& vi,
             ITensor& U, ITSparse& D, ITensor& V, int b = 1);

    void 
    svdRank2(const IQTensor& A, const IQIndex& uI, const IQIndex& vI,
             IQTensor& U, IQTSparse& D, IQTensor& V, int b = 1);

private:

    /*
    void
    diag_and_truncate(const IQTensor& rho, std::vector<Matrix>& mmatrix, 
                      std::vector<Vector>& mvector, std::vector<Real>& alleig, 
                      Real& svdtruncerr, IQIndex& newmid);
    void
    buildUnitary(const IQTensor& rho, const std::vector<Matrix>& mmatrix, 
                 const std::vector<Vector>& mvector,
                 const IQIndex& newmid, IQTensor& U);

    void
    buildCenter(const IQTensor& rho, const std::vector<Matrix>& mmatrix, 
                const std::vector<Vector>& mvector,
                const IQIndex& newmid, IQTensor& C);
                */

    ITensor 
    pseudoInverse(const ITensor& C, Real cutoff = 0);

    IQTensor 
    pseudoInverse(const IQTensor& C, Real cutoff = 0);

    /////////////////
    //
    // Data Members
    //

    int N;
    std::vector<Real> truncerr_;
    Real cutoff_;
    int minm_;
    int maxm_;
    bool use_orig_m_; 
    bool showeigs_;
    bool doRelCutoff_;
    bool absoluteCutoff_;
    LogNumber refNorm_;
    std::vector<Vector> eigsKept_;
    Real noise_;

    //
    /////////////////

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
                if(fabs(val) > cut_)
                    return 1./val;
                else
                    return 0;
                }
        private:
            Real cut_;
        };

    }; //class SVDWorker

template<class Tensor, class SparseT, class LocalOpT>
void SVDWorker::
csvd(int b, const Tensor& AA, Tensor& L, SparseT& V, Tensor& R, 
     const LocalOpT& PH)
    {
    typedef typename Tensor::IndexT 
    IndexT;
    typedef typename Tensor::CombinerT 
    CombinerT;

    Tensor UU(L),VV(R);
    SparseT D(V);
    svd(b,AA,UU,D,VV,PH);

    L = UU*D;
    R = D*VV;

    V = conj(D);
    V.pseudoInvert(0);

    } //void SVDWorker::csvd

template<class Tensor, class SparseT, class LocalOpT>
void SVDWorker::
svd(int b, Tensor AA, Tensor& U, SparseT& D, Tensor& V, 
    const LocalOpT& PH)
    {
    typedef typename Tensor::IndexT 
    IndexT;
    typedef typename Tensor::CombinerT 
    CombinerT;

    if(AA.vecSize() == 0) 
        throw ResultIsZero("denmatDecomp: AA.vecSize == 0");

    //Add in noise term
    if(noise_ > 0 && PH.isNotNull())
        {
        Real orig_norm = AA.norm();

        AA += noise_*PH.deltaPhi(AA);

        AA *= orig_norm/AA.norm();
        }

    //Combiners which transform AA
    //into a rank 2 tensor
    CombinerT Ucomb, Vcomb;
    Ucomb.doCondense(true);
    Vcomb.doCondense(true);

    //Divide up indices based on U
    //If U is null, use V instead
    const Tensor &L = (U.isNull() ? V : U);
    CombinerT &Lcomb = (U.isNull() ? Vcomb : Ucomb),
              &Rcomb = (U.isNull() ? Ucomb : Vcomb);
    for(int j = 1; j <= AA.r(); ++j) 
        { 
        const IndexT& I = AA.index(j);

        if(I == Tensor::ReImIndex()) 
            continue;

        if(L.hasindex(I))
            Lcomb.addleft(I);
        else
            Rcomb.addleft(I);
        }

    AA = Ucomb * AA * Vcomb;

    const Real saved_cutoff = cutoff_; 
    const int saved_minm = minm_,
              saved_maxm = maxm_; 
    if(use_orig_m_)
        {
        cutoff_ = -1;
        if(D.r() == 0)
            {
            Print(D);
            Error("D.r() must be > 1 to use original m option");
            }
        minm_ = D.index(1).m();
        maxm_ = D.index(1).m();
        }

    svdRank2(AA,Ucomb.right(),Vcomb.right(),U,D,V,b);

    cutoff_ = saved_cutoff; 
    minm_ = saved_minm; 
    maxm_ = saved_maxm; 

    U = conj(Ucomb) * U;
    V = V * conj(Vcomb);

    } // void SVDWorker::svd

template<class Tensor, class LocalOpT>
void SVDWorker::
denmatDecomp(int b, const Tensor& AA, Tensor& A, Tensor& B, Direction dir, 
             const LocalOpT& PH)
    {
    typedef typename Tensor::IndexT 
    IndexT;
    typedef typename Tensor::CombinerT 
    CombinerT;

    if(AA.vecSize() == 0) 
        throw ResultIsZero("denmatDecomp: AA.vecSize == 0");

    IndexT mid; 
    try {
        mid = index_in_common(A,B,Link);
        }
    catch(const ITError& e)
        {
        mid = IndexT("mid");
        }

    //If dir==None, put the O.C. on the side
    //that keeps mid's arrow the same
    bool do_edge_case = true;
    if(dir == None)
        {
        //std::cerr << boost::format("Arrow before = %s\n")%(mid.dir() == Out ? "Out" : "In");
        dir = (mid.dir() == Out ? Fromright : Fromleft);
        do_edge_case = false;
        }

    Tensor& to_orth = (dir==Fromleft ? A : B);
    Tensor& newoc   = (dir==Fromleft ? B : A);

    CombinerT comb;

    int unique_link = 0; //number of Links unique to to_orth
    for(int j = 1; j <= to_orth.r(); ++j) 
        { 
        const IndexT& I = to_orth.index(j);
        if(!(newoc.hasindex(I) || I == Tensor::ReImIndex() ))
            {
            if(I.type() == Link) ++unique_link;
            comb.addleft(I);
            }
        }

    //Check if we're at the edge
    if(unique_link == 0 && do_edge_case)
        {
        comb.init(mid.rawname());
        comb.product(AA,newoc);
        to_orth = comb; to_orth.conj();
        eigsKept_.at(b) = Vector(comb.right().m()); 
        eigsKept_.at(b) = 1.0/comb.right().m();
        return;
        }

    //Apply combiner
    comb.doCondense(true);
    comb.init(mid.rawname());
    Tensor AAc; comb.product(AA,AAc);

    const IndexT& active = comb.right();

    Tensor AAcc = conj(AAc); 
    AAcc.primeind(active); 
    Tensor rho = AAc*AAcc; 

    if(noise_ > 0 && PH.isNotNull())
        {
        rho += noise_*PH.deltaRho(rho,comb,dir);
        rho *= 1./trace(rho);
        }

    const Real saved_cutoff = cutoff_; 
    const int saved_minm = minm_,
              saved_maxm = maxm_; 
    if(use_orig_m_)
        {
        cutoff_ = -1;
        minm_ = mid.m();
        maxm_ = mid.m();
        }

    IndexT newmid;
    Tensor U;
    if(AAc.isComplex())
        truncerr_.at(b) = diag_denmat_complex(rho,eigsKept_.at(b),newmid,U);
    else
        truncerr_.at(b) = diag_denmat(rho,eigsKept_.at(b),newmid,U);

    cutoff_ = saved_cutoff; 
    minm_ = saved_minm; 
    maxm_ = saved_maxm; 

    comb.conj();
    comb.product(U,to_orth);
    newoc = conj(U) * AAc;

    } //void SVDWorker::denmatDecomp



#endif
