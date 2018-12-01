//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MPS_H
#define __ITENSOR_MPS_H
#include "itensor/decomp.h"
#include "itensor/mps/siteset.h"

namespace itensor {

class InitState;

class MPS
    {
    protected:
    int N_;
    mutable
    std::vector<ITensor> A_;
    int l_orth_lim_,
        r_orth_lim_;
    SiteSet sites_;
    mutable
    int atb_;
    std::string writedir_;
    bool do_write_;
    public:

    //
    // MPS Constructors
    //

    MPS();

    MPS(int N);

    MPS(SiteSet const& sites);

    MPS(InitState const& initState);

    MPS(MPS const& other);

    MPS&
    operator=(MPS const& other);

    ~MPS();

    //
    // MPS Accessor Methods
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
    ITensor const& 
    A(int i) const;

    void
    setA(int i, ITensor const& nA) { Aref(i) = nA; }

    void
    setA(int i, ITensor && nA) { Aref(i) = std::move(nA); }

    //Returns reference to i'th MPS tensor
    //which allows reading and writing
    ITensor& 
    Aref(int i);

    MPS&
    plusEq(MPS const& R, 
           Args const& args = Args::global());

    void 
    mapPrime(int oldp, int newp, TagSet const& ts = TagSet("All"));

    void 
    primelinks(int oldp, int newp);

    void 
    noprimelink();

    Spectrum 
    svdBond(int b, 
            ITensor const& AA, 
            Direction dir, 
            Args const& args = Args::global());

    template<class LocalOpT>
    Spectrum 
    svdBond(int b, 
            ITensor const& AA, 
            Direction dir, 
            LocalOpT const& PH, 
            Args const& args = Args::global());

    //Move the orthogonality center to site i 
    //(leftLim() == i-1, rightLim() == i+1, orthoCenter() == i)
    void 
    position(int i, Args args = Args::global());

    void 
    orthogonalize(Args const& args = Args::global());

    void
    swap(MPS& other);

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
    //MPS methods for writing to disk
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
    init_tensors(std::vector<ITensor>& A_, const InitState& initState);

    public:

    //
    // Advanced/Developer methods
    // Use with caution
    //

    void 
    rightLim(int val) { r_orth_lim_ = val; }

    void 
    leftLim(int val) { l_orth_lim_ = val; }

    }; //class MPS

template <typename MPSType>
MPSType&
addAssumeOrth(MPSType      & L,
              MPSType const& R, 
              Args const& args = Args::global());

//void 
//convertToIQ(const SiteSet& sites, const std::vector<ITensor>& A, 
//            std::vector<IQTensor>& qA, QN totalq = QN(), Real cut = 1E-12);

inline MPS& 
operator*=(MPS & psi, Real a) { psi.Aref(psi.leftLim()+1) *= a; return psi; }

inline MPS& 
operator/=(MPS & psi, Real a) { psi.Aref(psi.leftLim()+1) /= a; return psi; }

MPS inline
operator*(MPS psi, Real r) { psi *= r; return psi; }

MPS inline
operator*(Real r, MPS psi) { psi *= r; return psi; }

inline MPS& 
operator*=(MPS & psi, Cplx z) { psi.Aref(psi.leftLim()+1) *= z; return psi; }

inline MPS& 
operator/=(MPS & psi, Cplx z) { psi.Aref(psi.leftLim()+1) /= z; return psi; }

MPS inline
operator*(MPS psi, Cplx z) { psi *= z; return psi; }

MPS inline
operator*(Cplx z, MPS psi) { psi *= z; return psi; }

class InitState
    {
    public:
    using Storage = std::vector<IndexVal>;
    using String = std::string;

    InitState(SiteSet const& sites);

    InitState(SiteSet const& sites, String const& state);

    InitState& 
    set(int i, String const& state);

    InitState& 
    setAll(String const& state);

    IndexVal const&
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
// Other Methods Related to MPS
//

//MPS
//toMPS(IQMPS const& psi);

bool
isComplex(MPS const& psi);

bool
isOrtho(MPS const& psi);

int
orthoCenter(MPS const& psi);

Real
norm(MPS const& psi);

Real
normalize(MPS & psi);

template<typename MPSType>
Index
linkInd(MPSType const& psi, int b);

template<typename MPSType>
Index
rightLinkInd(MPSType const& psi, int i);

template<typename MPSType>
Index
leftLinkInd(MPSType const& psi, int i);

Real
averageM(MPS const& psi);

int
maxM(MPS const& psi);

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
void 
applyGate(ITensor const& gate, 
          MPS & psi,
          Args const& args = Args::global());

//Checks if A_[i] is left (left == true) 
//or right (left == false) orthogonalized
bool 
checkOrtho(MPS const& psi,
           int i, 
           bool left);

bool 
checkOrtho(MPS const& psi);

int 
findCenter(MPS const& psi);

bool 
checkQNs(MPS const& psi);

QN
totalQN(MPS const& psi);

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
void 
fitWF(MPS const& psi_basis, MPS & psi_to_fit);

template <class MPSType>
MPSType
sum(MPSType const& L, 
    MPSType const& R, 
    Args const& args = Args::global());


//
// Template method for efficiently summing 
// a set of MPS's or MPO's (or any class supporting operator+=)
// Performs the sum in a tree-like fashion in an attempt to
// leave the largest summands for the last few steps
//
// Assumes terms are zero-indexed
//
template <class MPSType>
MPSType
sum(std::vector<MPSType> const& terms, 
    Args const& args = Args::global());

std::ostream& 
operator<<(std::ostream& s, MPS const& M);

std::ostream& 
operator<<(std::ostream& s, InitState const& state);

} //namespace itensor

#include "mps_impl.h"


#endif
