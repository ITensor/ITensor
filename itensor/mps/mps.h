//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MPS_H
#define __ITENSOR_MPS_H
#include "itensor/svdalgs.h"
#include "itensor/mps/siteset.h"
//#include "itensor/mps/bondgate.h"

namespace itensor {

template <class Tensor>
class MPOt;

class InitState;

//void 
//convertToIQ(const SiteSet& sites, const std::vector<ITensor>& A, 
//            std::vector<IQTensor>& qA, QN totalq = QN(), Real cut = 1E-12);

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

template <class Tensor> class MPSt;

using MPS = MPSt<ITensor>;
using IQMPS = MPSt<IQTensor>;

template <class Tensor>
class MPSt
    {
    public:

    //
    //MPSt Constructors
    //

    MPSt();

    MPSt(const SiteSet& sites);

    MPSt(const InitState& initState);

    MPSt(const MPSt& other);

    MPSt&
    operator=(const MPSt& other);

    ~MPSt();

    //
    //MPSt Typedefs
    //

    using TensorT = Tensor;
    using IndexT = typename Tensor::index_type;
    using IndexValT = typename Tensor::indexval_type;
    using MPOType = MPOt<Tensor>;

    //
    //MPSt Accessor Methods
    //

    int 
    N() const { return N_;}

    int 
    rightLim() const { return r_orth_lim_; }
    void 
    rightLim(int val) { r_orth_lim_ = val; }

    int 
    leftLim() const { return l_orth_lim_; }
    void 
    leftLim(int val) { l_orth_lim_ = val; }

    bool 
    isOrtho() const { return (l_orth_lim_+1) == (r_orth_lim_-1); }

    int 
    orthoCenter() const;

    //Read-only access to i'th MPS tensor
    Tensor const& 
    A(int i) const;

    void
    setA(int i, Tensor const& nA) { Anc(i) = nA; }

    void
    setA(int i, Tensor && nA) { Anc(i) = std::move(nA); }

    //Returns reference to i'th MPS tensor
    //which allows reading and writing
    Tensor& 
    Anc(int i); //nc stands for non-const

    const SiteSet& 
    sites() const { return *sites_; }

    explicit operator bool() const { return valid(); }

    bool 
    valid() const { return (sites_!=0); }

    //
    //MPSt Operators
    //

    MPSt& 
    operator*=(Real a) { Anc(l_orth_lim_+1) *= a; return *this; }
    MPSt& 
    operator/=(Real a) { Anc(l_orth_lim_+1) /= a; return *this; }
    MPSt 
    operator*(Real r) const { MPSt res(*this); res *= r; return res; }

    MPSt& 
    operator*=(Complex z) { Anc(l_orth_lim_+1) *= z; return *this; }
    MPSt& 
    operator/=(Complex z) { Anc(l_orth_lim_+1) /= z; return *this; }
    MPSt 
    operator*(Complex z) const { MPSt res(*this); res *= z; return res; }

    MPSt&
    plusEq(const MPSt& R, 
           const Args& args = Global::args());

    void 
    mapprime(int oldp, int newp, IndexType type = All);

    void 
    primelinks(int oldp, int newp);

    void 
    noprimelink();

    Spectrum 
    svdBond(int b, const Tensor& AA, Direction dir, 
            const Args& args = Global::args());

    template <class LocalOpT>
    Spectrum 
    svdBond(int b, const Tensor& AA, Direction dir, 
                const LocalOpT& PH, const Args& args = Global::args());

    //Move the orthogonality center to site i 
    //(leftLim() == i-1, rightLim() == i+1, orthoCenter() == i)
    void 
    position(int i, const Args& args = Global::args());

    void 
    orthogonalize(const Args& args = Global::args());

    void 
    makeRealBasis(int j, const Args& args = Global::args());

    Real 
    norm() const;

    Real 
    normalize();

    bool 
    isComplex() const;

    //void 
    //toIQ(QN totalq, MPSt<IQTensor>& iqpsi, Real cut = 1E-12) const
    //    {
    //    iqpsi = MPSt<IQTensor>(*sites_);
    //    convertToIQ(*sites_,A_,iqpsi.A_,totalq,cut);
    //    }

    void
    swap(MPSt& other);

    bool
    doWrite() const { return do_write_; }
    void
    doWrite(bool val, const Args& args = Global::args());

    const std::string&
    writeDir() const { return writedir_; }

    //Read from a directory containing individual tensors,
    //as created when doWrite(true) is called.
    void 
    read(const std::string& dirname);

    void 
    read(std::istream& s);
    void 
    write(std::ostream& s) const;


    //
    // Deprecated methods, only for backwards compatibility
    //

    const SiteSet& 
    model() const { return *sites_; }

    protected:

    //////////////////////////

    int N_;

    mutable
    std::vector<Tensor> A_;

    int l_orth_lim_,
        r_orth_lim_;

    const SiteSet* sites_;

    mutable
    int atb_;

    std::string writedir_;

    bool do_write_;

    //////////////////////////

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
    initWrite(const Args& args = Global::args());
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

    MPSt&
    addAssumeOrth(const MPSt& R, 
                  const Args& args = Global::args());

    private:

    friend class MPSt<ITensor>;
    friend class MPSt<IQTensor>;

    }; //class MPSt<Tensor>


template <class Tensor>
MPSt<Tensor>
operator*(Real r, MPSt<Tensor> res) { res *= r; return res; }

template <class Tensor>
MPSt<Tensor>
operator*(Complex z, MPSt<Tensor> res) { res *= z; return res; }

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

    const SiteSet&
    sites() const { return *sites_; }

    private:

    const SiteSet* sites_;
    Storage state_;

    void
    checkRange(int i) const;
    }; 


//
// MPSt
// Template Methods
//

template <class Tensor>
template <class BigMatrixT>
Spectrum MPSt<Tensor>::
svdBond(int b, const Tensor& AA, Direction dir, 
        const BigMatrixT& PH, const Args& args)
    {
    setBond(b);
    if(dir == Fromleft && b-1 > l_orth_lim_)
        {
        printfln("b=%d, l_orth_lim_=%d",b,l_orth_lim_);
        Error("b-1 > l_orth_lim_");
        }
    if(dir == Fromright && b+2 < r_orth_lim_)
        {
        printfln("b=%d, r_orth_lim_=%d",b,r_orth_lim_);
        Error("b+2 < r_orth_lim_");
        }

    auto noise = args.getReal("Noise",0.);
    auto cutoff = args.getReal("Cutoff",MIN_CUT);
    auto usesvd = args.getBool("UseSVD",false);

    Spectrum res;

    if(usesvd || (noise == 0 && cutoff < 1E-12))
        {
        //Need high accuracy, use svd which calls the
        //accurate SVD method in the MatrixRef library
        Tensor D;
        res = svd(AA,A_[b],D,A_[b+1],args);

        //Normalize the ortho center if requested
        if(args.getBool("DoNormalize",false))
            {
            D *= 1./itensor::norm(D);
            }

        //Push the singular values into the appropriate site tensor
        if(dir == Fromleft) A_[b+1] *= D;
        else                A_[b]   *= D;
        }
    else
        {
        //If we don't need extreme accuracy
        //or need to use noise term
        //use density matrix approach
        res = denmatDecomp(AA,A_[b],A_[b+1],dir,PH,args);

        //Normalize the ortho center if requested
        if(args.getBool("DoNormalize",false))
            {
            Tensor& oc = (dir == Fromleft ? A_[b+1] : A_[b]);
            auto nrm = itensor::norm(oc);
            oc *= 1./nrm;
            }
        }

    if(dir == Fromleft)
        {
        l_orth_lim_ = b;
        if(r_orth_lim_ < b+2) r_orth_lim_ = b+2;
        }
    else //dir == Fromright
        {
        if(l_orth_lim_ > b-1) l_orth_lim_ = b-1;
        r_orth_lim_ = b+1;
        }

    return res;
    }

//
// Other Methods Related to MPSt
//

//
// projectOp takes the projected edge tensor W 
// of an operator and the site tensor X for the operator
// and creates the next projected edge tensor nE
//
// dir==Fromleft example:
//
//  /---A--     /---
//  |   |       |
//  E-- X -  =  nE -
//  |   |       |
//  \---A--     \---
//
template <class Tensor>
void 
projectOp(const MPSt<Tensor>& psi, int j, Direction dir, 
          const Tensor& E, const Tensor& X, Tensor& nE)
    {
    if(dir==Fromleft && j > psi.leftLim()) 
        { 
        printfln("projectOp: from left j > l_orth_lim_ (j=%d,leftLim=%d)",j,psi.leftLim());
        Error("Projecting operator at j > l_orth_lim_"); 
        }
    if(dir==Fromright && j < psi.rightLim()) 
        { 
        printfln("projectOp: from left j < r_orth_lim_ (j=%d,r_orth_lim_=%d)",j,psi.rightLim());
        Error("Projecting operator at j < r_orth_lim_"); 
        }
    nE = (E ? E*psi.A(j) : psi.A(j));
    nE *= X; 
    nE *= dag(prime(psi.A(j)));
    }


template <typename MPST>
typename MPST::IndexT 
linkInd(const MPST& psi, int b)
    { 
    return commonIndex(psi.A(b),psi.A(b+1),Link); 
    }

template <typename MPST>
typename MPST::IndexT 
rightLinkInd(const MPST& psi, int i)
    { 
    return commonIndex(psi.A(i),psi.A(i+1),Link); 
    }

template <typename MPST>
typename MPST::IndexT 
leftLinkInd(const MPST& psi, int i)
    { 
    return commonIndex(psi.A(i),psi.A(i-1),Link); 
    }

template <typename MPST>
int
averageM(MPST const& psi)
    {
    Real avgm = 0;
    for(int b = 1; b < psi.N(); ++b) avgm += linkInd(psi,b).m();
    avgm /= (psi.N()-1);
    return int(avgm);
    }

template <typename MPST>
int
maxM(MPST const& psi)
    {
    int maxM_ = 0;
    for(int b = 1; b < psi.N(); ++b) 
        {
        int mb = linkInd(psi,b).m();
        maxM_ = std::max(mb,maxM_);
        }
    return maxM_;
    }

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
          const Args& args = Global::args())
    {
    const int c = psi.orthoCenter();
    Tensor AA = psi.A(c) * psi.A(c+1) * gate;
    AA.noprime();
    psi.svdBond(c,AA,Fromleft,args);
    }

//template <class Tensor>
//void 
//applyGate(const BondGate<Tensor>& gate, 
//          MPSt<Tensor>& psi,
//          const Args& args = Global::args())
//    {
//    Tensor AA = psi.A(gate.i1()) * psi.A(gate.i1()+1) * Tensor(gate);
//    AA.noprime();
//    psi.svdBond(gate.i1(),AA,Fromleft,args);
//    }

//Checks if A_[i] is left (left == true) 
//or right (left == false) orthogonalized
template <class Tensor>
bool 
checkOrtho(const MPSt<Tensor>& psi,
           int i, 
           bool left)
    {
    using IndexT = typename Tensor::index_type;

    IndexT link = (left ? rightLinkInd(psi,i) : leftLinkInd(psi,i));
    Tensor rho = psi.A(i) * dag(prime(psi.A(i),link,4));
    Tensor Delta = makeKroneckerDelta(link,4);
    Tensor Diff = rho - Delta;

    const
    Real threshold = 1E-13;
    if(norm(Diff) < threshold) 
        {
        return true;
        }

    //Print any helpful debugging info here:
    println("checkOrtho: on line ",__LINE__," of mps.h,");
    println("checkOrtho: Tensor at position ",i," failed to be ",left?"left":"right"," ortho.");
    printfln("checkOrtho: norm(Diff) = %E",norm(Diff));
    printfln("checkOrtho: Error threshold set to %E",threshold);
    //-----------------------------

    return false;
    }

template <class Tensor>
bool 
checkOrtho(const MPSt<Tensor>& psi)
    {
    for(int i = 1; i <= psi.leftLim(); ++i)
    if(!checkOrtho(psi,i,true))
        {
        std::cout << "checkOrtho: A_[i] not left orthogonal at site i=" 
                  << i << std::endl;
        return false;
        }

    for(int i = psi.N(); i >= psi.rightLim(); --i)
    if(!checkOrtho(psi,i,false))
        {
        std::cout << "checkOrtho: A_[i] not right orthogonal at site i=" 
                  << i << std::endl;
        return false;
        }
    return true;
    }


int 
findCenter(const IQMPS& psi);

inline bool 
checkQNs(const MPS& psi) { return true; }

bool 
checkQNs(const IQMPS& psi);

QN
totalQN(const IQMPS& psi);

//
// <psi | phi>
//
template <class MPSType>
Cplx 
overlapC(MPSType const& psi, 
         MPSType const& phi)
    {
    auto N = psi.N();
    if(N != phi.N()) Error("overlap: mismatched N");

    auto l1 = linkInd(psi,1);
    auto L = phi.A(1);
    if(l1) 
        {
        L *= dag(prime(psi.A(1),l1)); 
        }
    else   
        {
        L *= dag(psi.A(1));
        }


    auto i = N;
    for(i = 2; i < N; ++i) 
        { 
        L = L * phi.A(i) * dag(prime(psi.A(i),Link)); 
        }
    L = L * phi.A(N);

    auto lNm = linkInd(psi,N-1);
    if(lNm) return (dag(prime(psi.A(N),lNm))*L).cplx();
    return (dag(psi.A(N))*L).cplx();
    }


template <class MPSType>
void 
overlap(MPSType const& psi,MPSType const& phi, Real& re, Real& im)
    {
    auto z = overlapC(psi,phi);
    re = z.real();
    im = z.imag();
    }

template <class MPSType>
Real 
overlap(MPSType const& psi, MPSType const& phi) //Re[<psi|phi>]
    {
    Real re, im;
    overlap(psi,phi,re,im);
    if(std::fabs(im) > (1E-12 * std::fabs(re)) )
        printfln("Real overlap: WARNING, dropping non-zero imaginary part (=%.5E) of expectation value.",im);
    return re;
    }

template <class MPSType>
Complex 
psiphiC(MPSType const& psi, 
        MPSType const& phi)
    {
    return overlapC(psi,phi);
    }

template <class MPSType>
void 
psiphi(MPSType const& psi,MPSType const& phi, Real& re, Real& im)
    {
    overlap(psi,phi,re,im);
    }

template <class MPSType>
Real 
psiphi(MPSType const& psi, MPSType const& phi) //Re[<psi|phi>]
    {
    return overlap(psi,phi);
    }

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
    Args const& args = Global::args())
    {
    MPSt<Tensor> res(L);
    res.plusEq(R,args);
    return res;
    }



//
// Template method for efficiently summing 
// a set of MPS's or MPO's (or any class supporting operator+=)
// Performs the sum in a tree-like fashion in an attempt to
// leave the largest summands for the last few steps
//
// Assumes terms are zero-indexed
//
template <typename MPSType>
MPSType 
sum(std::vector<MPSType> const& terms, 
    Args const& args = Global::args())
    {
    auto Nt = terms.size();
    if(Nt == 2)
        { 
        return sum(terms.at(0),terms.at(1),args);
        }
    else 
    if(Nt == 1) 
        {
        return terms.at(0);
        }
    else 
    if(Nt > 2)
        {
        //Add all MPS in pairs
        auto nsize = (Nt%2==0 ? Nt/2 : (Nt-1)/2+1);
        std::vector<MPSType> newterms(nsize); 
        for(decltype(Nt) n = 0, np = 0; n < Nt-1; n += 2, ++np)
            {
            newterms.at(np) = sum(terms.at(n),terms.at(n+1),args);
            }
        if(Nt%2 == 1) newterms.at(nsize-1) = terms.back();

        //Recursively call sum again
        return sum(newterms,args);
        }
    return MPSType();
    }

template <class Tensor>
std::ostream& 
operator<<(std::ostream& s, const MPSt<Tensor>& M);

std::ostream& 
operator<<(std::ostream& s, const InitState& state);

} //namespace itensor


#endif
