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
#ifndef __ITENSOR_MPS_H
#define __ITENSOR_MPS_H
#include "itensor/decomp.h"
#include "itensor/mps/siteset.h"

namespace itensor {

// Some forward definitions
class InitState;

class MPS
    {
    protected:
    int N_;
    mutable
    std::vector<ITensor> A_;
    int l_orth_lim_,
        r_orth_lim_;
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

    MPS(SiteSet const& sites, int m = 1);

    MPS(IndexSet const& sites, int m = 1);

    MPS(InitState const& initState);

    MPS(MPS const& other);

    MPS&
    operator=(MPS const& other);

    ~MPS();

    //
    // MPS Accessor Methods
    //

    int 
    length() const { return N_;}

    // Deprecated in favor of length()
    int 
    N() const
        {
        Global::warnDeprecated(".N() is deprecated in favor of length(MPS)");
        return this->length();
        }

    explicit operator bool() const { return (not A_.empty()); }

    int 
    rightLim() const { return r_orth_lim_; }

    int 
    leftLim() const { return l_orth_lim_; }

    // Read-only access to i'th MPS tensor
    // Is 1-indexed
    ITensor const&
    operator()(int i) const;

    //Returns reference to i'th MPS tensor
    //which allows reading and writing
    //Is 1-indexed
    ITensor& 
    ref(int i);

    void
    set(int i, ITensor const& nA) { ref(i) = nA; }

    void
    set(int i, ITensor && nA) { ref(i) = std::move(nA); }

    Real
    normalize();

    MPS&
    plusEq(MPS const& R, 
           Args const& args = Args::global());

    // Dagger all MPS tensors
    MPS&
    dag()
      {
      for(auto i : range1(N_))
          A_[i].dag();
      return *this;
      }

    MPS&
    replaceSiteInds(IndexSet const& sites);

    MPS&
    replaceLinkInds(IndexSet const& links);

    //
    //MPS Index Methods
    //

    MPS&
    setTags(TagSet const& ts, IndexSet const& is)
        {
        if(do_write_)
            Error("setTags not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].setTags(ts,is);
        return *this;
        }

    template<typename... VarArgs>
    MPS&
    setTags(VarArgs&&... vargs)
        {
        if(do_write_)
            Error("setTags not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].setTags(std::forward<VarArgs>(vargs)...);
        return *this;
        }

    MPS&
    noTags(IndexSet const& is)
        {
        if(do_write_)
            Error("noTags not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].noTags(is);
        return *this;
        }

    template<typename... VarArgs>
    MPS& 
    noTags(VarArgs&&... vargs)
        {
        if(do_write_)
            Error("noTags not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].noTags(std::forward<VarArgs>(vargs)...);
        return *this;
        }

    MPS&
    addTags(TagSet const& ts, IndexSet const& is)
        {
        if(do_write_)
            Error("addTags not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].addTags(ts,is);
        return *this;
        }

    template<typename... VarArgs>
    MPS&
    addTags(VarArgs&&... vargs)
        {
        if(do_write_)
            Error("addTags not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].addTags(std::forward<VarArgs>(vargs)...);
        return *this;
        }

    MPS&
    removeTags(TagSet const& ts, IndexSet const& is)
        {
        if(do_write_)
            Error("removeTags not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].removeTags(ts,is);
        return *this;
        }

    template<typename... VarArgs>
    MPS&
    removeTags(VarArgs&&... vargs)
        {
        if(do_write_)
            Error("removeTags not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].removeTags(std::forward<VarArgs>(vargs)...);
        return *this;
        }

    MPS&
    replaceTags(TagSet const& ts1, TagSet const& ts2, IndexSet const& is)
        {
        if(do_write_)
            Error("replaceTags not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].replaceTags(ts1,ts2,is);
        return *this;
        }

    template<typename... VarArgs>
    MPS&
    replaceTags(VarArgs&&... vargs)
        {
        if(do_write_)
            Error("replaceTags not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].replaceTags(std::forward<VarArgs>(vargs)...);
        return *this;
        }

    MPS&
    swapTags(TagSet const& ts1, TagSet const& ts2, IndexSet const& is)
        {
        if(do_write_)
            Error("swapTags not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].swapTags(ts1,ts2,is);
        return *this;
        }

    template<typename... VarArgs>
    MPS&
    swapTags(VarArgs&&... vargs)
        {
        if(do_write_)
            Error("swapTags not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].swapTags(std::forward<VarArgs>(vargs)...);
        return *this;
        }

    MPS&
    prime(int plev, IndexSet const& is)
        {
        if(do_write_)
            Error("prime not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].prime(plev,is);
        return *this;
        }

    MPS&
    prime(IndexSet const& is)
        {
        if(do_write_)
            Error("prime not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].prime(is);
        return *this;
        }

    template<typename... VarArgs>
    MPS&
    prime(VarArgs&&... vargs)
        {
        if(do_write_)
            Error("prime not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].prime(std::forward<VarArgs>(vargs)...);
        return *this;
        }

    MPS&
    setPrime(int plev, IndexSet const& is)
        {
        if(do_write_)
            Error("setPrime not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].setPrime(plev,is);
        return *this;
        }

    template<typename... VarArgs>
    MPS&
    setPrime(VarArgs&&... vargs)
        {
        if(do_write_)
            Error("setPrime not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].setPrime(std::forward<VarArgs>(vargs)...);
        return *this;
        }

    MPS&
    mapPrime(int plevold, int plevnew, IndexSet const& is)
        {
        if(do_write_)
            Error("mapPrime not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].mapPrime(plevold,plevnew,is);
        return *this;
        }

    template<typename... VarArgs>
    MPS&
    mapPrime(VarArgs&&... vargs)
        {
        if(do_write_)
            Error("mapPrime not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].mapPrime(std::forward<VarArgs>(vargs)...);
        return *this;
        }

    MPS&
    swapPrime(int plevold, int plevnew, IndexSet const& is)
        {
        if(do_write_)
            Error("swapPrime not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].swapPrime(plevold,plevnew,is);
        return *this;
        }

    template<typename... VarArgs>
    MPS&
    swapPrime(VarArgs&&... vargs)
        {
        if(do_write_)
            Error("swapPrime not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].swapPrime(std::forward<VarArgs>(vargs)...);
        return *this;
        }

    MPS&
    noPrime(IndexSet const& is)
        {
        if(do_write_)
            Error("noPrime not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].noPrime(is);
        return *this;
        }

    template<typename... VarArgs>
    MPS&
    noPrime(VarArgs&&... vargs)
        {
        if(do_write_)
            Error("noPrime not supported if doWrite(true)");
        for(int i = 1; i <= N_; ++i)
            A_[i].noPrime(std::forward<VarArgs>(vargs)...);
        return *this;
        }

    // Randomize the tensors of the MPS
    // Right now, only supports randomizing dim = 1 MPS
    MPS&
    randomize(Args const& args = Args::global());

    Spectrum 
    svdBond(int b, 
            ITensor const& AA, 
            Direction dir, 
            Args args = Args::global());

    template<class LocalOpT>
    Spectrum 
    svdBond(int b, 
            ITensor const& AA, 
            Direction dir, 
            LocalOpT const& PH, 
            Args args = Args::global());

    //Move the orthogonality center to site i 
    //(leftLim() == i-1, rightLim() == i+1, orthoCenter() == i)
    MPS& 
    position(int i, Args args = Args::global());

    MPS& 
    orthogonalize(Args args = Args::global());

    void
    swap(MPS & other);

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

    //Returns an unsafe reference to i'th MPS tensor
    //which allows reading and writing
    //Modifying this reference may break orthogonality!
    ITensor&
    uref(int i);

    //
    // Deprecated
    //

    void
    setA(int i, ITensor const& nA) { this->set(i,nA); }

    void
    setA(int i, ITensor && nA) { this->set(i,nA); }

    ITensor const& 
    A(int i) const;

    ITensor& 
    Aref(int i);

    }; //class MPS

template <typename MPSType>
MPSType&
addAssumeOrth(MPSType      & L,
              MPSType const& R, 
              Args const& args = Args::global());

MPS& 
operator*=(MPS & x, Real a);

MPS& 
operator/=(MPS & x, Real a);

MPS
operator*(MPS x, Real r);

MPS
operator*(Real r, MPS x);

MPS& 
operator*=(MPS & x, Cplx z);

MPS& 
operator/=(MPS & x, Cplx z);

MPS
operator*(MPS x, Cplx z);

MPS
operator*(Cplx z, MPS x);

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

MPS
dag(MPS A);

//
// MPS tag functions
//

MPS
setTags(MPS A, TagSet const& ts, IndexSet const& is);

template<typename... VarArgs>
MPS
setTags(MPS A,
        VarArgs&&... vargs)
    {
    A.setTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

MPS
noTags(MPS A, IndexSet const& is);

template<typename... VarArgs>
MPS
noTags(MPS A,
       VarArgs&&... vargs)
    {
    A.noTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

MPS
addTags(MPS A, TagSet const& ts, IndexSet const& is);

template<typename... VarArgs>
MPS
addTags(MPS A,
        VarArgs&&... vargs)
    {
    A.addTags(std::forward<VarArgs>(vargs)...);
    return A;
    }
    
MPS
removeTags(MPS A, TagSet const& ts, IndexSet const& is);

template<typename... VarArgs>
MPS
removeTags(MPS A,
           VarArgs&&... vargs)
    {
    A.removeTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

MPS
replaceTags(MPS A, TagSet const& ts1, TagSet const& ts2, IndexSet const& is);

template<typename... VarArgs>
MPS
replaceTags(MPS A,
            VarArgs&&... vargs)
    {
    A.replaceTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

MPS
swapTags(MPS A, TagSet const& ts1, TagSet const& ts2, IndexSet const& is);

template<typename... VarArgs>
MPS
swapTags(MPS A,
         VarArgs&&... vargs)
    {
    A.swapTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

MPS
prime(MPS A, int plev, IndexSet const& is);

MPS
prime(MPS A, IndexSet const& is);

template<typename... VarArgs>
MPS
prime(MPS A,
      VarArgs&&... vargs)
    {
    A.prime(std::forward<VarArgs>(vargs)...);
    return A; 
    }

MPS
setPrime(MPS A, int plev, IndexSet const& is);

template<typename... VarArgs>
MPS
setPrime(MPS A,
         VarArgs&&... vargs)
    {
    A.setPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

MPS
mapPrime(MPS A, int plevold, int plevnew, IndexSet const& is);

template<typename... VarArgs>
MPS
mapPrime(MPS A,
         VarArgs&&... vargs)
    {
    A.mapPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

MPS
swapPrime(MPS A, int plevold, int plevnew, IndexSet const& is);

template<typename... VarArgs>
MPS
swapPrime(MPS A,
          VarArgs&&... vargs)
    {
    A.swapPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

MPS
noPrime(MPS A, IndexSet const& is);

template<typename... VarArgs>
MPS
noPrime(MPS A,
        VarArgs&&... vargs)
    {
    A.noPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

//
// Other Methods Related to MPS
//

int 
length(MPS const& W);

bool
hasQNs(InitState const& initstate);

template <class MPSType>
bool
hasQNs(MPSType const& x);

//Create a random MPS
MPS
randomMPS(SiteSet const& sites,
          int m = 1,
          Args const& args = Args::global());

MPS
randomMPS(SiteSet const& sites,
          Args const& args);

MPS
randomMPS(InitState const& initstate,
          Args const& args = Args::global());

//Remove the QNs of each tensor of the MPS
template <class MPSType>
MPSType
removeQNs(MPSType const& x);

bool
isComplex(MPS const& x);

bool
isOrtho(MPS const& x);

int
orthoCenter(MPS const& x);

int
rightLim(MPS const& x);

int
leftLim(MPS const& x);

Real
norm(MPS const& x);

bool
hasSiteInds(MPS const& x, IndexSet const& sites);

template <class MPSType>
IndexSet
siteInds(MPSType const& W, int b);

IndexSet
siteInds(MPS const& x);

MPS
replaceSiteInds(MPS x, IndexSet const& sites);

MPS
replaceLinkInds(MPS x, IndexSet const& links);

Index
siteIndex(MPS const& x, int j);

template<typename MPSType>
Index
leftLinkIndex(MPSType const& x, int b);

template<typename MPSType>
Index
rightLinkIndex(MPSType const& x, int b);

template<typename MPSType>
Index
linkIndex(MPSType const& x, int b);

template<typename MPSType>
IndexSet
linkInds(MPSType const& x, int b);

template<typename MPSType>
IndexSet
linkInds(MPSType const& x);

template<typename MPSType>
Real
averageLinkDim(MPSType const& x);

template<typename MPSType>
int
maxLinkDim(MPSType const& x);

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
          MPS & x,
          Args const& args = Args::global());

//Checks if A_[i] is left (left == true) 
//or right (left == false) orthogonalized
template <typename MPSType>
bool 
checkOrtho(MPSType const& A,
           int i, 
           bool left,
           Real threshold = 1E-13);

template <typename MPSType>
bool
checkOrtho(MPSType const& A,
           Real threshold = 1E-13);

int 
findCenter(MPS const& x);

bool 
checkQNs(MPS const& x);

QN
totalQN(MPS const& x);

// Re[<x|y>]
Real 
inner(MPS const& x, MPS const& y);

// <x|y>
Cplx 
innerC(MPS const& x, 
       MPS const& y);

// <x|y>
void 
inner(MPS const& x,
      MPS const& y, 
      Real& re, Real& im);

//Computes an MPS which has the same overlap with x_basis as x_to_fit,
//but which differs from x_basis only on the first site, and has same index
//structure as x_basis. Result is stored to x_to_fit on return.
void 
fitWF(MPS const& x_basis, MPS & x_to_fit);

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

    
//
//  Try and verify the SiteSet and MPS are mutually consistent.
//     We say "try" because even if d matches they could still correspond to different
//     Hilbert spaces or different basis sets (less likely).
//  
bool inline
checkConsistent(const MPS& psi, 
                const SiteSet& sites
                )
{
    bool consistent=sites.length()==psi.length();
    for (int j=1;j<=sites.length();j++)
       consistent = consistent && sites(j).dim()==siteIndex(psi,j).dim();
    return consistent;
}
//-------------------------------------------------------------------------------------------
//
//  Template implementation of expect function for Real and Complex types and fixed site list.
//
template <class T> std::vector<std::vector<T>> 
expectT(const MPS& _psi, 
        const SiteSet& sites,
        const std::vector<string>& vops,
        std::vector<int> site_list 
        )
{
    assert(checkConsistent(_psi,sites));
//
// Work with copy because we need to move the orth-center
//        
    MPS psi=_psi; 
    if (!isOrtho(psi)) psi.orthogonalize();
    psi.normalize(); //Is this expensive if we are already orthogonalized?

    std::vector<std::vector<T>> ex;
    for (auto i:site_list)
    {
        psi.position(i); //Set the ortho centre.
        std::vector<T> exi; //row of data for site i
        for (auto str_op:vops)
        {
            auto e=psi(i) * sites.op(str_op,i) * dag(prime(psi(i), sites(i)));
            exi.push_back(eltT<T>(e));
        }
        ex.push_back(exi); //push row into the table.
    }

    return ex;
}

void fixRange(detail::RangeHelper<int>&r,int N);
//
//  Convert range to explicit list in one place
//
template <class T> std::vector<std::vector<T>> 
expectT(const MPS& _psi, 
        const SiteSet& sites,
        const std::vector<string>& vops,
        detail::RangeHelper<int> site_range=range1(0) //fake default because we don't have access to sites.length in the function signature.
        )
{
    fixRange(site_range,_psi.length());
    std::vector<int> site_list;
    for (auto i:site_range) site_list.push_back(i);
    return expectT<T>(_psi,sites,vops,site_list);
}

//  1D containers for returning arrays of numbers from expect().
typedef std::vector<Real   > VecR;  
typedef std::vector<Complex> VecC;  
//  2D containers for returning tables of numbers from expect() and correlation functions
typedef std::vector<VecR> VecVecR;  
typedef std::vector<VecC> VecVecC;  

//
//  User versions of expect (Real) and expectC(Complex).  These just function forward to the template than does the work.
//  We need hand code all 4 combinations of {Real,Complex}(X){range,list}
//
inline VecVecR 
expect (const MPS& _psi, 
        const SiteSet& sites,
        const std::vector<string>& vops,
        detail::RangeHelper<int> site_range=range1(0) //fake default because we don't have access to sites.length in the function signature.
        )
{
    return expectT<Real>(_psi,sites,vops,site_range);
}


inline VecVecC 
expectC(const MPS& _psi, 
        const SiteSet& sites,
        const std::vector<string>& vops,
        detail::RangeHelper<int> site_range=range1(0) //fake default because we don't have access to sites.length in the function signature.
        )
{
    return expectT<Complex>(_psi,sites,vops,site_range);
}

inline VecVecR 
expect (const MPS& _psi, 
        const SiteSet& sites,
        const std::vector<string>& vops,
        std::vector<int> site_list 
        )
{
    return expectT<Real>(_psi,sites,vops,site_list);
}

inline VecVecC 
expectC(const MPS& _psi, 
        const SiteSet& sites,
        const std::vector<string>& vops,
        std::vector<int> site_list 
        )
{
    return expectT<Complex>(_psi,sites,vops,site_list);
}

//
//  Single operator versions
//
template <typename T> 
std::vector<T> 
getFirstColumn(const std::vector<std::vector<T>>& m)
{
   std::vector<T> col0;
   for (auto r:m)
        col0.push_back(r[0]);
   return col0;
}

inline VecR 
expect (const MPS& _psi, 
        const SiteSet& sites,
        const std::string& vop,
        detail::RangeHelper<int> site_range=range1(0) //fake default because we don't have access to sites.length in the function signature.
        )
{
    std::vector<string> vops;
    vops.push_back(vop);
    return getFirstColumn(expectT<Real>(_psi,sites,vops,site_range));
}

inline VecC 
expectC(const MPS& _psi, 
        const SiteSet& sites,
        const std::string& vop,
        detail::RangeHelper<int> site_range=range1(0) //fake default because we don't have access to sites.length in the function signature.
        )
{
    std::vector<string> vops;
    vops.push_back(vop);
    return getFirstColumn(expectT<Complex>(_psi,sites,vops,site_range));
}

inline VecR 
expect (const MPS& _psi, 
        const SiteSet& sites,
        const std::string& vop,
        std::vector<int> site_list 
        )
{
    std::vector<string> vops;
    vops.push_back(vop);
    return getFirstColumn(expectT<Real>(_psi,sites,vops,site_list));
}

inline VecC 
expectC(const MPS& _psi, 
        const SiteSet& sites,
        const std::string& vop,
        std::vector<int> site_list 
        )
{
    std::vector<string> vops;
    vops.push_back(vop);
    return getFirstColumn(expectT<Complex>(_psi,sites,vops,site_list));
}


//-------------------------------------------------------------------------------------------
//
//  Template implementation of correlationMatrix function for Real and Complex types.
//
using std::conj;
inline Real conj(Real r) {return r;} //dummy conj so we can compile generic function

//  Template function that does all the work.  See mps.cc for implementation.
template <class T> std::vector<std::vector<T>>
correlationMatrixT(const MPS& _psi,
                   const SiteSet& sites,
                   const string& op1,
                   const string& op2,
                   detail::RangeHelper<int> site_range, //No defaults
                   Args const& args //No defaults
                   );
//
//  Real and Complex wrappers for various combinations of default args.
//
inline VecVecR
correlationMatrix(const MPS& psi,
                  const SiteSet& sites,
                  const string& op1,
                  const string& op2,
                  detail::RangeHelper<int> site_range,
                  Args const& args = Args::global()
                 )
{
    return correlationMatrixT<Real>(psi,sites,op1,op2,site_range,args);
}


inline VecVecR
correlationMatrix(const MPS& psi,
                  const SiteSet& sites,
                  const string& op1,
                  const string& op2,
                  Args const& args = Args::global()
                 )
{
    return correlationMatrixT<Real>(psi,sites,op1,op2,range1(0),args);
}

inline VecVecC
correlationMatrixC(const MPS& psi,
                   const SiteSet& sites,
                   const string& op1,
                   const string& op2,
                   detail::RangeHelper<int> site_range,
                  Args const& args = Args::global()
                  )
{
    return correlationMatrixT<Complex>(psi,sites,op1,op2,site_range,args);
}
inline VecVecC
correlationMatrixC(const MPS& psi,
                   const SiteSet& sites,
                   const string& op1,
                   const string& op2,
                   Args const& args = Args::global()
                  )
{
    return correlationMatrixT<Complex>(psi,sites,op1,op2,range1(0),args);
}


std::ostream& 
operator<<(std::ostream& s, MPS const& M);

std::ostream& 
operator<<(std::ostream& s, InitState const& state);

//
// Deprecated
//

int
maxM(MPS const& x);

Real
averageM(MPS const& x);

// Deprecated in favor of x.normalize()
Real
normalize(MPS & x);

template <class MPSType>
Real 
overlap(MPSType const& psi, MPSType const& phi);

template <class MPSType>
Cplx 
overlapC(MPSType const& psi, 
         MPSType const& phi);

template <class MPSType>
void 
overlap(MPSType const& psi,
        MPSType const& phi, 
        Real& re, Real& im);

Spectrum
orthMPS(ITensor& A1, ITensor& A2, Direction dir, Args const& args);

#ifdef ITENSOR_USE_HDF5
void
h5_write(h5::group parent, std::string const& name, MPS const& M);
void
h5_read(h5::group parent, std::string const& name, MPS & M);
#endif

} //namespace itensor

#include "mps_impl.h"


#endif
