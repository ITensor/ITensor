#ifndef __ITENSOR_SVDWORKER_H
#define __ITENSOR_SVDWORKER_H
#include "iqcombiner.h"

template <class Tensor>
class ProjectedOp;

enum Direction { Fromright, Fromleft, Both, None };

template<class Tensor, class IndexT>
IndexT
index_in_common(const Tensor& A, const Tensor& B, IndexType t)
    {
    for(int j = 1; j <= A.r(); ++j)
        {
        const IndexT& I = A.index(j);
        if(I.type() == t && B.hasindex(I)) { return I; }
        }
    return IndexT();
    }

Index inline
index_in_common(const ITensor& A, const ITensor& B, IndexType t)
    { return index_in_common<ITensor,Index>(A,B,t); }

IQIndex inline
index_in_common(const IQTensor& A, const IQTensor& B, IndexType t)
    { return index_in_common<IQTensor,IQIndex>(A,B,t); }

//
// SVDWorker
//

class SVDWorker
    {
    public:

    //Constructors ---------------

    SVDWorker();

    SVDWorker(int N_);

    SVDWorker(int N_, Real cutoff, int minm, int maxm, 
              bool doRelCutoff, const LogNumber& refNorm);

    SVDWorker(std::istream& s) { read(s); }

    //Canonical/Bond SVD (weights in V tensor) ---------------

    template <class Tensor> 
    void 
    operator()(const Tensor& AA, Tensor& L, Tensor& V, Tensor& R)
        { 
        operator()<Tensor>(1,AA,L,V,R); 
        }

    template <class Tensor> 
    void 
    operator()(int b, const Tensor& AA, Tensor& L, Tensor& V, Tensor& R);

    template <class Tensor> 
    void 
    operator()(int b, const ProjectedOp<Tensor>& PH, 
               const Tensor& AA, Tensor& L, Tensor& V, Tensor& R);


    //Site SVD -----------------------------------------------
    //(weights in A or B tensor, for dir = Fromleft or Fromright, respectively)

    template <class Tensor> 
    void 
    operator()(const Tensor& AA, Tensor& A, Tensor& B, Direction dir)
        { 
        operator()<Tensor>(1,AA,A,B,dir); 
        }

    template <class Tensor> 
    void
    operator()(int b, const Tensor& AA, Tensor& A, Tensor& B, Direction dir);

    template <class Tensor> 
    void
    operator()(int b, const ProjectedOp<Tensor>& PH, 
               const Tensor& AA, Tensor& A, Tensor& B, Direction dir);

    //Accessors ---------------

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

private:

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

    ITensor 
    pseudoInverse(const ITensor& C, Real cutoff = 0);

    IQTensor 
    pseudoInverse(const IQTensor& C, Real cutoff = 0);

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

    }; //class SVDWorker


#endif
