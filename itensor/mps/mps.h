//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MPS_H
#define __ITENSOR_MPS_H
#include "itensor/decomp.h"
#include "itensor/mps/siteset.h"

namespace itensor {

template <class Tensor>
class MPOt;

class InitState;

//
// class MPSt
// (the lowercase t stands for "template")
// 
// Unless there is a need to use MPSt
// specifically as a templated class,
// it is recommended to use one of the 
// typedefs of MPSt:
//
//      MPS for ITensors
//    IQMPS for IQTensors
//

template <class Tensor> 
class MPSt;

using MPS = MPSt<ITensor>;
using IQMPS = MPSt<IQTensor>;

template<class Tensor>
class MPSt
    {
    protected:
    int N_;
    mutable
    std::vector<Tensor> A_;
    int l_orth_lim_,
        r_orth_lim_;
    SiteSet sites_;
    mutable
    int atb_;
    std::string writedir_;
    bool do_write_;
    public:
    using TensorT = Tensor;
    using IndexT = typename Tensor::index_type;
    using IndexValT = typename Tensor::indexval_type;
    using MPOType = MPOt<Tensor>;

    //
    // MPSt Constructors
    //

    MPSt();

    MPSt(int N);

    MPSt(SiteSet const& sites);

    MPSt(InitState const& initState);

    MPSt(MPSt const& other);

    MPSt&
    operator=(MPSt const& other);

    ~MPSt();

    //
    // MPSt Accessor Methods
    //

    int 
    N() const { return N_;}

    SiteSet const& 
    sites() const;

    explicit operator bool() const { return (not A_.empty()); }

    int 
    rightLim() const { return r_orth_lim_; }

    int 
    leftLim() const { return l_orth_lim_; }

    //Read-only access to i'th MPS tensor
    Tensor const& 
    A(int i) const;

    void
    setA(int i, Tensor const& nA) { Aref(i) = nA; }

    void
    setA(int i, Tensor && nA) { Aref(i) = std::move(nA); }

    //Returns reference to i'th MPS tensor
    //which allows reading and writing
    Tensor& 
    Aref(int i);
    Tensor& 
    Anc(int i) { return Aref(i); }

    MPSt&
    plusEq(MPSt const& R, 
           Args const& args = Args::global());

    void 
    mapprime(int oldp, int newp, IndexType type = All);

    void 
    primelinks(int oldp, int newp);

    void 
    noprimelink();

    Spectrum 
    svdBond(int b, 
            Tensor const& AA, 
            Direction dir, 
            Args const& args = Args::global());

    template<class LocalOpT>
    Spectrum 
    svdBond(int b, 
            Tensor const& AA, 
            Direction dir, 
            LocalOpT const& PH, 
            Args const& args = Args::global());

    //Move the orthogonality center to site i 
    //(leftLim() == i-1, rightLim() == i+1, orthoCenter() == i)
    void 
    position(int i, Args const& args = Args::global());

    void 
    orthogonalize(Args const& args = Args::global());

    void
    swap(MPSt& other);

    bool
    doWrite() const { return do_write_; }

    void
    doWrite(bool val, const Args& args = Args::global());

    std::string const&
    writeDir() const { return writedir_; }

    //Read from a directory containing individual tensors,
    //as created when doWrite(true) is called.
    void 
    read(std::string const& dirname);

    void 
    read(std::istream& s);

    void 
    write(std::ostream& s) const;

    protected:

    //
    //MPSt methods for writing to disk
    //

    //if doWrite(true) is called
    //setBond(b) loads bond b
    //from disk, keeping all other
    //tensors written to disk
    void
    setBond(int b) const;

    void
    setSite(int j) const;

    void
    initWrite(const Args& args = Args::global());
    void
    copyWriteDir();
    void
    cleanupWrite();

    std::string
    AFName(int j, const std::string& dirname = "") const;

    //
    //Constructor Helpers
    //

    void 
    new_tensors(std::vector<ITensor>& A_);

    void 
    random_tensors(std::vector<ITensor>& A_);

    void 
    random_tensors(std::vector<IQTensor>& A_) { }

    void 
    init_tensors(std::vector<ITensor>& A_, const InitState& initState);

    void 
    init_tensors(std::vector<IQTensor>& A_, const InitState& initState);


    private:

    friend class MPSt<ITensor>;
    friend class MPSt<IQTensor>;

    public:

    //
    // Advanced/Developer methods
    // Use with caution
    //

    void 
    rightLim(int val) { r_orth_lim_ = val; }

    void 
    leftLim(int val) { l_orth_lim_ = val; }


    //
    // Deprecated methods
    // 

    //prefer function norm(psi) instead
    Real 
    norm() const;

    //prefer function normalize(psi) instead
    Real 
    normalize();


    //prefer isOrtho(psi) instead
    bool 
    isOrtho() const { return leftLim()+1 == rightLim()-1; }

    //prefer orthoCenter(psi) instead
    int 
    orthoCenter() const;

    //prefer isComplex(psi) instead
    bool 
    isComplex() const;

    }; //class MPSt<Tensor>

template <class MPSType>
MPSType&
addAssumeOrth(MPSType      & L,
              MPSType const& R, 
              Args const& args = Args::global());

//void 
//convertToIQ(const SiteSet& sites, const std::vector<ITensor>& A, 
//            std::vector<IQTensor>& qA, QN totalq = QN(), Real cut = 1E-12);

template<class T>
MPSt<T>& 
operator*=(MPSt<T> & psi, Real a) { psi.Aref(psi.leftLim()+1) *= a; return psi; }

template<class T>
MPSt<T>& 
operator/=(MPSt<T> & psi, Real a) { psi.Aref(psi.leftLim()+1) /= a; return psi; }

template<class T>
MPSt<T>
operator*(MPSt<T> psi, Real r) { psi *= r; return psi; }

template <class T>
MPSt<T>
operator*(Real r, MPSt<T> psi) { psi *= r; return psi; }

template<class T>
MPSt<T>& 
operator*=(MPSt<T> & psi, Cplx z) { psi.Aref(psi.leftLim()+1) *= z; return psi; }

template<class T>
MPSt<T>& 
operator/=(MPSt<T> & psi, Cplx z) { psi.Aref(psi.leftLim()+1) /= z; return psi; }

template <class T>
MPSt<T>
operator*(MPSt<T> psi, Cplx z) { psi *= z; return psi; }

template <class T>
MPSt<T>
operator*(Cplx z, MPSt<T> psi) { psi *= z; return psi; }

class InitState
    {
    public:

    using Storage = std::vector<IQIndexVal>;

    using String = std::string;

    InitState(const SiteSet& sites);

    InitState(const SiteSet& sites, const String& state);

    InitState& 
    set(int i, const String& state);

    InitState& 
    setAll(const String& state);

    const IQIndexVal&
    operator()(int i) const { checkRange(i); return state_.at(i); }

    SiteSet const&
    sites() const { return sites_; }

    private:

    SiteSet sites_;
    Storage state_;

    void
    checkRange(int i) const;
    }; 


//
// Other Methods Related to MPSt
//

MPS
toMPS(IQMPS const& psi);

template<typename T>
bool
isComplex(MPSt<T> const& psi);

template<typename T>
bool
isOrtho(MPSt<T> const& psi);

template<typename T>
int
orthoCenter(MPSt<T> const& psi);

template<typename T>
Real
norm(MPSt<T> const& psi);

template<typename T>
Real
normalize(MPSt<T> & psi);

template <typename MPST>
typename MPST::IndexT 
linkInd(MPST const& psi, int b);

template <typename MPST>
typename MPST::IndexT 
rightLinkInd(MPST const& psi, int i);

template <typename MPST>
typename MPST::IndexT 
leftLinkInd(MPST const& psi, int i);

template <typename MPST>
Real
averageM(MPST const& psi);

template <typename MPST>
int
maxM(MPST const& psi);

//
// Applies a bond gate to the bond that is currently
// the OC.                                    |      |
// After calling position b, this bond is - A_[b] - A_[b+1] -
//
//      |      |
//      ==gate==
//      |      |
//  - A_[b] - A_[b+1] -
//
// Does not normalize the resulting wavefunction unless 
// Args("DoNormalize",true) is included in args.
template <class Tensor>
void 
applyGate(const Tensor& gate, 
          MPSt<Tensor>& psi,
          const Args& args = Args::global());

//Checks if A_[i] is left (left == true) 
//or right (left == false) orthogonalized
template <class Tensor>
bool 
checkOrtho(const MPSt<Tensor>& psi,
           int i, 
           bool left);

template <class T>
bool 
checkOrtho(MPSt<T> const& psi);

int 
findCenter(IQMPS const& psi);

bool inline
checkQNs(MPS const& psi) { return true; }

bool 
checkQNs(IQMPS const& psi);

QN
totalQN(IQMPS const& psi);

// Re[<psi|phi>]
template <class MPSType>
Real 
overlap(MPSType const& psi, MPSType const& phi);

// <psi|phi>
template <class MPSType>
Cplx 
overlapC(MPSType const& psi, 
         MPSType const& phi);

// <psi|phi>
template <class MPSType>
void 
overlap(MPSType const& psi,
        MPSType const& phi, 
        Real& re, Real& im);

//Computes an MPS which has the same overlap with psi_basis as psi_to_fit,
//but which differs from psi_basis only on the first site, and has same index
//structure as psi_basis. Result is stored to psi_to_fit on return.
template <class Tensor>
void 
fitWF(const MPSt<Tensor>& psi_basis, MPSt<Tensor>& psi_to_fit);

template <class Tensor>
MPSt<Tensor>
sum(MPSt<Tensor> const& L, 
    MPSt<Tensor> const& R, 
    Args const& args = Args::global());


//
// Template method for efficiently summing 
// a set of MPS's or MPO's (or any class supporting operator+=)
// Performs the sum in a tree-like fashion in an attempt to
// leave the largest summands for the last few steps
//
// Assumes terms are zero-indexed
//
template<typename MPSType>
MPSType 
sum(std::vector<MPSType> const& terms, 
    Args const& args = Args::global());

template<class Tensor>
std::ostream& 
operator<<(std::ostream& s, const MPSt<Tensor>& M);

std::ostream& 
operator<<(std::ostream& s, const InitState& state);

} //namespace itensor

#include "mps_impl.h"


#endif
