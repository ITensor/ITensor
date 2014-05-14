//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MPS_H
#define __ITENSOR_MPS_H
#include "svdalgs.h"
#include "model.h"
#include "boost/function.hpp"
#include "bondgate.h"

#define Cout std::cout
#define Endl std::endl
#define Format boost::format

namespace itensor {

template <class Tensor>
class MPOt;

class InitState;

void 
convertToIQ(const Model& model, const std::vector<ITensor>& A, 
            std::vector<IQTensor>& qA, QN totalq = QN(), Real cut = 1E-12);

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
class MPSt
    {
    public:

    //
    //MPSt Constructors
    //

    MPSt();

    MPSt(const Model& sites);

    MPSt(const InitState& initState);

    MPSt(const MPSt& other);

    MPSt&
    operator=(const MPSt& other);

    ~MPSt();

    //
    //MPSt Typedefs
    //

    typedef Tensor 
    TensorT;

    typedef typename Tensor::IndexT 
    IndexT;

    typedef typename Tensor::IndexValT 
    IndexValT;

    typedef MPOt<Tensor>
    MPOType;

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

    IQIndex 
    si(int i) const { return model_->si(i); }

    IQIndex 
    siP(int i) const { return model_->siP(i); }

    typedef typename std::vector<Tensor>::const_iterator 
    AA_it;

    //Returns pair of iterators which can be used in a 
    //Foreach loop to iterate over all MPS tensors
    const std::pair<AA_it,AA_it> 
    A() const { return std::make_pair(A_.begin()+1,A_.end()); }

    //Read-only access to i'th MPS tensor
    const Tensor& 
    A(int i) const;

    //Returns reference to i'th MPS tensor
    //which allows reading and writing
    Tensor& 
    Anc(int i); //nc stands for non-const

    const Model& 
    model() const { return *model_; }

    const Spectrum& 
    spectrum(int b) const { return spectrum_.at(b); }

    Spectrum& 
    spectrum(int b) { return spectrum_.at(b); }

    bool 
    isNull() const { return (model_==0); }
    bool 
    isNotNull() const { return (model_!=0); }

    Real 
    truncerr(int b) const { return spectrum_.at(b).truncerr(); }

    const Vector& 
    eigsKept(int b) const { return spectrum_.at(b).eigsKept(); }

    bool
    doWrite() const { return do_write_; }
    void
    doWrite(bool val, const OptSet& opts = Global::opts());

    const std::string&
    writeDir() const { return writedir_; }

    bool 
    isOrtho() const { return (l_orth_lim_+1) == (r_orth_lim_-1); }

    void 
    read(std::istream& s);
    void 
    write(std::ostream& s) const;

    //Read from a directory containing individual tensors,
    //as created when doWrite(true) is called.
    void 
    read(const std::string& dirname);

    //
    //MPSt Operators
    //

    MPSt& 
    operator*=(Real a) { Anc(l_orth_lim_+1) *= a; return *this; }
    MPSt& 
    operator/=(Real a) { Anc(l_orth_lim_+1) /= a; return *this; }
    MPSt 
    operator*(Real r) const { MPSt res(*this); res *= r; return res; }
    friend inline MPSt 
    operator*(Real r, MPSt res) { res *= r; return res; }

    MPSt& 
    operator*=(Complex z) { Anc(l_orth_lim_+1) *= z; return *this; }
    MPSt& 
    operator/=(Complex z) { Anc(l_orth_lim_+1) /= z; return *this; }
    MPSt 
    operator*(Complex z) const { MPSt res(*this); res *= z; return res; }
    friend inline MPSt 
    operator*(Complex z, MPSt res) { res *= z; return res; }

    MPSt&
    plusEq(const MPSt& R, 
           const OptSet& opts = Global::opts());

    //
    //MPSt Index Methods
    //

    void 
    mapprime(int oldp, int newp, IndexType type = All);

    void 
    primelinks(int oldp, int newp);

    void 
    noprimelink();

    IndexT 
    LinkInd(int b) const 
        { return commonIndex(A(b),A(b+1),Link); }
    IndexT 
    RightLinkInd(int i) const 
        { return commonIndex(A(i),A(i+1),Link); }
    IndexT 
    LeftLinkInd(int i)  const 
        { return commonIndex(A(i),A(i-1),Link); }

    //
    //MPSt orthogonalization methods
    //

    void 
    svdBond(int b, const Tensor& AA, Direction dir, 
            const OptSet& opts = Global::opts());

    template <class LocalOpT>
    void 
    svdBond(int b, const Tensor& AA, Direction dir, 
                const LocalOpT& PH, const OptSet& opts = Global::opts());

    //Move the orthogonality center to site i 
    //(l_orth_lim_ = i-1, r_orth_lim_ = i+1)
    void 
    position(int i, const OptSet& opts = Global::opts());

    int 
    orthoCenter() const 
        { 
        if(!isOrtho()) Error("orthogonality center not well defined.");
        return (l_orth_lim_ + 1);
        }

    void 
    orthogonalize(const OptSet& opts = Global::opts());

    void 
    makeRealBasis(int j, const OptSet& opts = Global::opts());

    //Checks if A_[i] is left (left == true) 
    //or right (left == false) orthogonalized
    bool 
    checkOrtho(int i, bool left) const;

    bool 
    checkRightOrtho(int i) const { return checkOrtho(i,false); }

    bool 
    checkLeftOrtho(int i) const { return checkOrtho(i,true); }
    
    bool 
    checkOrtho() const;


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
    // Opt DoNormalize(true) is included in opts.
    void 
    applygate(const Tensor& gate, const OptSet& opts = Global::opts());

    void 
    applygate(const BondGate<Tensor>& gate, const OptSet& opts = Global::opts());

    Real 
    norm() const;

    int
    averageM() const;

    Real 
    normalize();

    bool 
    isComplex() const
        { 
        for(int j = 1; j <= N_; ++j)
            {
            if(A_[j].isComplex()) return true;
            }
        return false;
        }

    void 
    toIQ(QN totalq, MPSt<IQTensor>& iqpsi, Real cut = 1E-12) const
        {
        iqpsi = MPSt<IQTensor>(*model_);
        iqpsi.spectrum_ = spectrum_;
        convertToIQ(*model_,A_,iqpsi.A_,totalq,cut);
        }

    void
    swap(MPSt& other);

    protected:

    //////////////////////////
    //
    //Data Members

    int N_;

    mutable
    std::vector<Tensor> A_;

    int l_orth_lim_,
        r_orth_lim_;

    const Model* model_;

    std::vector<Spectrum> spectrum_;

    mutable
    int atb_;

    std::string writedir_;

    bool do_write_;

    //
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
    initWrite(const OptSet& opts = Global::opts());
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
                  const OptSet& opts = Global::opts());

    private:

    friend class MPSt<ITensor>;
    friend class MPSt<IQTensor>;

    }; //class MPSt<Tensor>
typedef MPSt<ITensor> MPS;
typedef MPSt<IQTensor> IQMPS;

class InitState
    {
    public:

    typedef std::vector<IQIndexVal>
    Storage;

    typedef std::string
    String;

    InitState(const Model& sites)
        : 
        sites_(&sites), 
        state_(1+sites.N())
        { 
        for(int n = 1; n <= sites_->N(); ++n)
            {
            state_[n] = sites_->si(n)(1);
            }
        }

    InitState(const Model& sites, const String& state)
        : 
        sites_(&sites), 
        state_(1+sites.N())
        { 
        setAll(state);
        }

    InitState& 
    set(int i, const String& state)
        { 
        checkRange(i);
        state_.at(i) = sites_->st(i,state);
        return *this;
        }

    InitState& 
    setAll(const String& state)
        { 
        for(int n = 1; n <= sites_->N(); ++n)
            {
            state_[n] = sites_->st(n,state);
            }
        return *this;
        }

    const IQIndexVal&
    operator()(int i) const { checkRange(i); return state_.at(i); }

    const Model&
    model() const { return *sites_; }

    private:

    const Model* sites_;
    Storage state_;

    void
    checkRange(int i) const
        {
        if(i > sites_->N() || i < 1) 
            {
            Cout << "i = " << i << Endl;
            Cout << "Valid range is 1 to " << sites_->N() << Endl;
            Error("i out of range");
            }
        }
    }; 


//
// MPSt
// Template Methods
//

template <class Tensor>
template <class BigMatrixT>
void MPSt<Tensor>::
svdBond(int b, const Tensor& AA, Direction dir, 
        const BigMatrixT& PH, const OptSet& opts)
    {
    setBond(b);

    if(dir == Fromleft && b-1 > l_orth_lim_)
        {
        Cout << Format("b=%d, l_orth_lim_=%d")
                %b%l_orth_lim_ << Endl;
        Error("b-1 > l_orth_lim_");
        }
    if(dir == Fromright && b+2 < r_orth_lim_)
        {
        Cout << Format("b=%d, r_orth_lim_=%d")
                %b%r_orth_lim_ << Endl;
        Error("b+2 < r_orth_lim_");
        }

    const Real noise = opts.getReal("Noise",0.);
    const Real cutoff = opts.getReal("Cutoff",MIN_CUT);

    if(opts.getBool("UseSVD",false) || (noise == 0 && cutoff < 1E-12))
        {
        //Need high accuracy, use svd which calls the
        //accurate SVD method in the MatrixRef library
        Tensor D;
        /*
        if(Global::debug1())
            {
            Cout << "Calling svdBond SVD" << Endl;
            Cout << "with spectrum:\n" << spectrum_.at(b) << Endl;
            }
            */
        spectrum_.at(b) = svd(AA,A_[b],D,A_[b+1],opts);

        //Normalize the ortho center if requested
        if(opts.getBool("DoNormalize",false))
            {
            D *= 1./D.norm();
            }

        //Push the singular values into the appropriate site tensor
        if(dir == Fromleft)
            A_[b+1] *= D;
        else
            A_[b] *= D;
        }
    else
        {
        //If we don't need extreme accuracy
        //or need to use noise term
        //use density matrix approach
        //Cout << "Calling svdBond denmatDecomp" << Endl;
        spectrum_.at(b) = 
        denmatDecomp(AA,A_[b],A_[b+1],dir,PH,opts);

        //Normalize the ortho center if requested
        if(opts.getBool("DoNormalize",false))
            {
            Tensor& oc = (dir == Fromleft ? A_[b+1] : A_[b]);
            Real norm = oc.norm();
            oc *= 1./norm;
            }
        }

    if(dir == Fromleft)
        {
        l_orth_lim_ = b;
        if(r_orth_lim_ < b+2) 
            {
            r_orth_lim_ = b+2;
            }
        }
    else //dir == Fromright
        {
        if(l_orth_lim_ > b-1) 
            {
            l_orth_lim_ = b-1;
            }
        r_orth_lim_ = b+1;
        }
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
        Cout << Format("projectOp: from left j > l_orth_lim_ (j=%d,leftLim=%d)")%j%psi.leftLim() << Endl;
        Error("Projecting operator at j > l_orth_lim_"); 
        }
    if(dir==Fromright && j < psi.rightLim()) 
        { 
        Cout << Format("projectOp: from left j < r_orth_lim_ (j=%d,r_orth_lim_=%d)")%j%psi.rightLim() << Endl;
        Error("Projecting operator at j < r_orth_lim_"); 
        }
    nE = (E.isNull() ? psi.A(j) : E * psi.A(j));
    nE *= X; 
    nE *= conj(prime(psi.A(j)));
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
Complex 
psiphiC(const MPSType& psi, const MPSType& phi)
    {
    typedef typename MPSType::TensorT
    Tensor;
    typedef typename Tensor::IndexT
    IndexT;

    const int N = psi.N();
    if(N != phi.N()) Error("psiphi: mismatched N");

    IndexT l1 = psi.LinkInd(1);
    Tensor L = phi.A(1);
    if(l1.isNull())
        L *= conj(psi.A(1));
    else
        L *= conj(prime(psi.A(1),l1)); 

    for(int i = 2; i < N; ++i) 
        { 
        L = L * phi.A(i) * conj(prime(psi.A(i),Link)); 
        }
    L = L * phi.A(N);

    Complex z;

    IndexT lNm = psi.LinkInd(N-1);
    if(lNm.isNull())
        z = BraKet(psi.A(N),L);
    else
        z = BraKet(prime(psi.A(N),lNm),L);

    return z;
    }

template <class MPSType>
void 
psiphi(const MPSType& psi, const MPSType& phi, Real& re, Real& im)
    {
    Complex z = psiphiC(psi,phi);
    re = z.real();
    im = z.imag();
    }

template <class MPSType>
Real 
psiphi(const MPSType& psi, const MPSType& phi) //Re[<psi|phi>]
    {
    Real re, im;
    psiphi(psi,phi,re,im);
    if(fabs(im) > (1E-12 * fabs(re)) )
        Cout << Format("Real psiphi: WARNING, dropping non-zero imaginary part (=%.5E) of expectation value.") % im << Endl;
    return re;
    }

//Computes an MPS which has the same overlap with psi_basis as psi_to_fit,
//but which differs from psi_basis only on the first site, and has same index
//structure as psi_basis. Result is stored to psi_to_fit on return.
template <class Tensor>
void 
fitWF(const MPSt<Tensor>& psi_basis, MPSt<Tensor>& psi_to_fit);

template <class Tensor>
MPSt<Tensor>
sum(const MPSt<Tensor>& L, 
    const MPSt<Tensor>& R, 
    const OptSet& opts = Global::opts())
    {
    MPSt<Tensor> res(L);
    res.plusEq(R,opts);
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
sum(const std::vector<MPSType>& terms, 
    const OptSet& opts = Global::opts())
    {
    const int Nt = terms.size();
    if(Nt == 2)
        { 
        return sum(terms.at(0),terms.at(1),opts);
        }
    else 
    if(Nt == 1) 
        {
        return terms.at(0);
        }
    else 
    if(Nt > 2)
        {
        //Add all MPS's in pairs
        const int nsize = (Nt%2==0 ? Nt/2 : (Nt-1)/2+1);
        std::vector<MPSType> newterms(nsize); 
        for(int n = 0, np = 0; n < Nt-1; n += 2, ++np)
            {
            newterms.at(np) = sum(terms.at(n),terms.at(n+1),opts);
            }
        if(Nt%2 == 1) newterms.at(nsize-1) = terms.back();

        //Recursively call sum again
        return sum(newterms,opts);
        }
    return MPSType();
    }

template <class Tensor>
std::ostream& 
operator<<(std::ostream& s, const MPSt<Tensor>& M)
    {
    s << "\n";
    for(int i = 1; i <= M.N(); ++i) 
        s << M.A(i) << "\n";
    return s;
    }

std::ostream& 
operator<<(std::ostream& s, const InitState& state);

}; //namespace itensor

#undef Cout
#undef Endl
#undef Format

#endif
