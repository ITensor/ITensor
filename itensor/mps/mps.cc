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
#include "itensor/mps/mps.h"
#include "itensor/mps/mpo.h"
#include "itensor/mps/localop.h"
#include "itensor/util/print_macro.h"
#include "itensor/util/str.h"
#include "itensor/tensor/algs.h"

namespace itensor {

using std::istream;
using std::ostream;
using std::cout;
using std::endl;
using std::vector;
using std::find;
using std::pair;
using std::make_pair;
using std::string;

void
new_tensors(std::vector<ITensor>& A,
            SiteSet const& sites,
            int m = 1);

void
new_tensors(std::vector<ITensor>& A,
            IndexSet const& sites,
            int m = 1);

//
// class MPS
//

//
// Constructors
//

MPS::
MPS() 
    : 
    N_(0), 
    atb_(1),
    writedir_("./"),
    do_write_(false)
    { }

MPS::
MPS(int N)
    : 
    N_(N), 
    A_(N+2), //idmrg may use A_[0] and A[N+1]
    l_orth_lim_(0),
    r_orth_lim_(N+1),
    atb_(1),
    writedir_("./"),
    do_write_(false)
    { 
    }

MPS::
MPS(SiteSet const& sites,
    int m)
    : 
    N_(sites.length()), 
    A_(sites.length()+2), //idmrg may use A_[0] and A[N+1]
    l_orth_lim_(0),
    r_orth_lim_(sites.length()+1),
    atb_(1),
    writedir_("./"),
    do_write_(false)
    { 
    new_tensors(A_,sites,m);
    }

MPS::
MPS(IndexSet const& sites,
    int m)
    :
    N_(sites.length()),
    A_(sites.length()+2), //idmrg may use A_[0] and A[N+1]
    l_orth_lim_(0),
    r_orth_lim_(sites.length()+1),
    atb_(1),
    writedir_("./"),
    do_write_(false)
    {
    new_tensors(A_,sites,m);
    }

MPS::
MPS(InitState const& initState)
    : 
    N_(initState.sites().length()),
    A_(initState.sites().length()+2), //idmrg may use A_[0] and A[N+1]
    l_orth_lim_(0),
    r_orth_lim_(2),
    atb_(1),
    writedir_("./"),
    do_write_(false)
    { 
    init_tensors(A_,initState);
    }

MPS::
MPS(MPS const& other)
    : 
    N_(other.N_),
    A_(other.A_),
    l_orth_lim_(other.l_orth_lim_),
    r_orth_lim_(other.r_orth_lim_),
    atb_(other.atb_),
    writedir_(other.writedir_),
    do_write_(other.do_write_)
    { 
    copyWriteDir();
    }

MPS& MPS::
operator=(MPS const& other)
    { 
    N_ = other.N_;
    A_ = other.A_;
    l_orth_lim_ = other.l_orth_lim_;
    r_orth_lim_ = other.r_orth_lim_;
    atb_ = other.atb_;
    writedir_ = other.writedir_;
    do_write_ = other.do_write_;

    copyWriteDir();
    return *this;
    }

MPS::
~MPS()
    {
    cleanupWrite();
    }

// TODO: Add note in error to call randomMPS(InitState,dim) instead, once it is implemented
MPS& MPS::
randomize(Args const& args)
    {
    if(maxLinkDim(*this)>1) Error(".randomize() not currently supported on MPS with bond dimension greater than 1."); 
    for(auto i : range1(N_)) ref(i).randomize(args);
    return *this;
    }

Real MPS::
normalize()
    {
    auto nrm = norm(*this);
    if(std::fabs(nrm) < 1E-20) Error("Zero norm");
    *this /= nrm;
    return nrm;
    }

MPS&
operator*=(MPS & x, Real a) { x.ref(x.leftLim()+1) *= a; return x; }

MPS&
operator/=(MPS & x, Real a) { x.ref(x.leftLim()+1) /= a; return x; }

MPS
operator*(MPS x, Real r) { x *= r; return x; }

MPS
operator*(Real r, MPS x) { x *= r; return x; }

MPS&
operator*=(MPS & x, Cplx z) { x.ref(x.leftLim()+1) *= z; return x; }

MPS&
operator/=(MPS & x, Cplx z) { x.ref(x.leftLim()+1) /= z; return x; }

MPS
operator*(MPS x, Cplx z) { x *= z; return x; }

MPS
operator*(Cplx z, MPS x) { x *= z; return x; }

int
length(MPS const& W)
    {
    return W.length();
    }

Matrix
randomOrthog(int n, int m)
    {
    auto r = std::max(n,m);
    auto M = randn(r,r);
    Matrix Q,R;
    QR(M,Q,R,{"PositiveDiagonal",true});
    if(m < n)
        {
        reduceCols(Q,m);
        }
    else if(n < m)
        {
        reduceCols(Q,n);
        return Matrix(transpose(Q));
        }
    return Q;
    }

MPS
randomCircuitMPS(SiteSet const& s, int m, Args const& args)
    {
    auto N = length(s);
    auto M = MPS(N);
    auto l = vector<Index>(N+1);

    //Make N'th MPS tensor
    int chi = dim(s(N));
    l[N-1] = Index(chi,format("Link,l=%d",N-1));
    auto O = randomOrthog(chi,dim(s(N)));
    M.ref(N) = matrixITensor(O,l[N-1],s(N));

    for(int j = N-1; j > 1; j -= 1)
        {
        //Make j'th MPS tensor
        auto prev_chi = chi;
        chi *= dim(s(j));
        chi = std::min(m,chi);
        l[j-1] = Index(chi,format("Link,l=%d",j-1));
        O = randomOrthog(chi,prev_chi*dim(s(j)));
        auto [C,c] = combiner(s(j),l[j]);
        M.ref(j) = matrixITensor(O,l[j-1],c);
        M.ref(j) *= C;
        }

    //Make 1st MPS tensor
    O = randomOrthog(1,dim(s(1))*chi);
    auto [C,c] = combiner(s(1),l[1]);
    l[0] = Index(1,"Link,l=0");
    M.ref(1) = matrixITensor(O,l[0],c);
    M.ref(1) *= C;
    M.ref(1) *= setElt(l[0](1));

    M.leftLim(0);
    M.rightLim(2);
    return M;
    }

MPS
randomMPS(SiteSet const& sites, int m, Args const& args)
    {
    if(hasQNs(sites))
    {
    Error("randomMPS(SiteSet) with QN conservation is ambiguous, use randomMPS(InitState) instead.");
    }

    return randomCircuitMPS(sites,m,args);
    }

MPS
randomMPS(SiteSet const& sites, Args const& args)
    {
    return randomMPS(sites,1,args);
    }

//TODO: implement for m>1 in terms of random gates
MPS
randomMPS(InitState const& initstate, int m, Args const& args)
    {
    if(m>1) Error("randomMPS(InitState,m>1) not currently supported.");
    auto psi = MPS(initstate);
    psi.randomize(args);
    return psi;
    }

MPS
randomMPS(InitState const& initstate, Args const& args)
    {
    return randomMPS(initstate,1,args);
    }

ITensor const& MPS::
operator()(int i) const
    { 
    if(i < 0) i = N_+i+1;
    setSite(i);
    return A_.at(i); 
    }

ITensor& MPS::
ref(int i)
    { 
    if(i < 0) i = N_+i+1;
    setSite(i);
    if(i <= l_orth_lim_) l_orth_lim_ = i-1;
    if(i >= r_orth_lim_) r_orth_lim_ = i+1;
    return A_.at(i); 
    }

ITensor& MPS::
uref(int i)
    {
    if(do_write_)
        Error("replaceTags not supported if doWrite(true)");
    return A_.at(i);
    }

// Deprecated
ITensor const& MPS::
A(int i) const
    { 
    return this->operator()(i);
    }

// Deprecated
ITensor& MPS::
Aref(int i)
    { 
    return this->ref(i);
    }

void MPS::
doWrite(bool val, Args const& args) 
    { 
    if(val == do_write_) return;

    if(val == true)
        {
        initWrite(args); 
        }
    else
        {
        read(writedir_);
        cleanupWrite();
        }
    }


void MPS::
read(std::istream & s)
    {
    itensor::read(s,N_);
    A_.resize(N_+2);
    for(auto j : range(A_))
        {
        itensor::read(s,A_[j]);
        }
    itensor::read(s,l_orth_lim_);
    itensor::read(s,r_orth_lim_);
    }

void MPS::
write(std::ostream& s) const
    {
    if(do_write_)
        Error("MPS::write not yet supported if doWrite(true)");

    itensor::write(s,length());
    for(auto j : range(A_.size()))
        {
        itensor::write(s,A_[j]);
        }
    itensor::write(s,leftLim());
    itensor::write(s,rightLim());
    }

void MPS::
read(std::string const& dirname)
    {
    l_orth_lim_ = 0;
    r_orth_lim_ = N_+1;

    //std::string dname_ = dirname;
    //if(dname_[dname_.length()-1] != '/')
    //    dname_ += "/";

    for(auto j : range(A_.size()))
        {
    	readFromFile(AFName(j,dirname),A_.at(j));
        }
    }


string MPS::
AFName(int j, string const& dirname) const
    { 
    if(dirname == "")
        {
        return format("%s/A_%03d",writedir_,j);
        }
    else
        {
        return format("%s/A_%03d",dirname,j);
        }
    }

void MPS::
setBond(int b) const
    {
    if(b == atb_) return;
    if(!do_write_)
        {
        atb_ = b;
        return;
        }
    if(b < 1 || b >= N_) return;

    //
    //Shift atb_ (location of bond that is loaded into RAM)
    //to requested value b, writing any non-Null tensors to
    //disk along the way
    //
    while(b > atb_)
        {
        if(A_.at(atb_))
            {
            writeToFile(AFName(atb_),A_.at(atb_));
            A_.at(atb_) = ITensor();
            }
        if(A_.at(atb_+1))
            {
            writeToFile(AFName(atb_+1),A_.at(atb_+1));
            if(atb_+1 != b) A_.at(atb_+1) = ITensor();
            }
        ++atb_;
        }
    while(b < atb_)
        {
        if(A_.at(atb_))
            {
            writeToFile(AFName(atb_),A_.at(atb_));
            if(atb_ != b+1) A_.at(atb_) = ITensor();
            }
        if(A_.at(atb_+1))
            {
            writeToFile(AFName(atb_+1),A_.at(atb_+1));
            A_.at(atb_+1) = ITensor();
            }
        --atb_;
        }
    assert(atb_ == b);
    //
    //Load tensors at bond b into RAM if
    //they aren't loaded already
    //
    if(!A_.at(b))
        {
        readFromFile(AFName(b),A_.at(b));
        }

    if(!A_.at(b+1))
        {
        readFromFile(AFName(b+1),A_.at(b+1));
        }

    //if(b == 1)
        //{
        //writeToFile(writedir_+"/sites",*sites_);
        //std::ofstream inf((format("%s/info")%writedir_).str().c_str());
        //    inf.write((char*) &l_orth_lim_,sizeof(l_orth_lim_));
        //    inf.write((char*) &r_orth_lim_,sizeof(r_orth_lim_));
        //    svd_.write(inf);
        //inf.close();
        //}
    }

void MPS::
setSite(int j) const
    {
    if(!do_write_)
        {
        atb_ = (j > atb_ ? j-1 : j);
        return;
        }
    if(j < 1 || j > N_) return;

    if(j < atb_)
        {
        //Cout << Format("j=%d < atb_=%d, calling setBond(%d)")
        //        % j % atb_ % j << Endl;
        setBond(j);
        }
    else
    if(j > atb_+1)
        {
        //Cout << Format("j=%d > atb_+1=%d, calling setBond(%d)")
        //        % j % (atb_+1) % (j-1) << Endl;
        setBond(j-1);
        }

    //otherwise the set bond already
    //contains this site
    }


void
new_tensors(std::vector<ITensor>& A,
            SiteSet const& sites,
            int m)
    {
    auto N = length(sites);
    auto a = std::vector<Index>(N+1);
    if(hasQNs(sites))
        {
        if(m==1) for(auto i : range1(N)) a[i] = Index(QN(),m,format("Link,l=%d",i));
        else Error("Cannot create QN conserving MPS with bond dimension greater than 1 from a SiteSet");
        }
    else
        {
        for(auto i : range1(N)) a[i] = Index(m,format("Link,l=%d",i));
        }
    A[1] = ITensor(sites(1),a[1]);
    for(int i = 2; i < N; i++)
        { 
        A[i] = ITensor(dag(a[i-1]),sites(i),a[i]); 
        }
    A[N] = ITensor(dag(a[N-1]),sites(N));
    }

void
new_tensors(std::vector<ITensor>& A,
            IndexSet const& sites,
            int m)
    {
    auto N = length(sites);
    auto a = std::vector<Index>(N+1);
    if(hasQNs(sites))
        {
        if(m==1) for(auto i : range1(N)) a[i] = Index(QN(),m,format("Link,l=%d",i));
        else Error("Cannot create QN conserving MPS with bond dimension greater than 1 from an IndexSet");
        }
    else
        {
        for(auto i : range1(N)) a[i] = Index(m,format("Link,l=%d",i));
        }
    A[1] = ITensor(sites(1),a[1]);
    for(int i = 2; i < N; i++)
        {
        A[i] = ITensor(dag(a[i-1]),sites(i),a[i]);
        }
    A[N] = ITensor(dag(a[N-1]),sites(N));
    }

void MPS::
init_tensors(std::vector<ITensor>& A_, InitState const& initState)
    { 
    auto a = std::vector<Index>(N_+1);
    if(hasQNs(initState))
        {
        auto qa = std::vector<QN>(N_+1); //qn[i] = qn on i^th bond
        for(auto i : range1(N_)) qa[0] -= qn(initState(i))*In;
        //Taking OC to be at the leftmost site,
        //compute the QuantumNumbers of all the Links.
        for(auto i : range1(N_))
            {
            //Taking the divergence to be zero,solve for qa[i]
            qa[i] = Out*(-qa[i-1]*In - qn(initState(i)));
            }
        for(auto i : range1(N_)) a[i] = Index(qa[i],1,format("Link,l=%d",i));
        }
    else
        {
        for(auto i : range1(N_)) a[i] = Index(1,format("Link,l=%d",i));
        }

    A_[1] = setElt(initState(1),a[1](1));
    for(auto i : range(2,N_))
        {
        A_[i] = setElt(itensor::dag(a[i-1])(1),initState(i),a[i](1));
        }
    A_[N_] = setElt(itensor::dag(a[N_-1])(1),initState(N_));
    }


MPS& MPS::
plusEq(MPS const& R,
       Args const& args)
    {
    //cout << "calling new orthog in sum" << endl;
    if(!itensor::isOrtho(*this))
        {
        try { 
            orthogonalize();
            }
        catch(ResultIsZero const& rz) 
            { 
            *this = R;
            return *this;
            }
        }

    if(!itensor::isOrtho(R))
        {
        MPS oR(R);
        try { 
            oR.orthogonalize(); 
            }
        catch(ResultIsZero const& rz) 
            { 
            return *this;
            }
        return addAssumeOrth(*this,oR,args);
        }

    return addAssumeOrth(*this,R,args);
    }

Spectrum MPS::
svdBond(int b, ITensor const& AA, Direction dir, Args args)
    {
    return svdBond(b,AA,dir,LocalOp(),args);
    }

struct SqrtInv
    {
    Real
    operator()(Real val) const 
        { 
        if(val == 0) return 0;
        return 1./std::sqrt(std::fabs(val)); 
        }
    };

struct Sqrt
    {
    Real
    operator()(Real val) const { return std::sqrt(std::fabs(val)); }
    };

Spectrum
orthMPS(ITensor& A1, ITensor& A2, Direction dir, Args const& args)
    {
    ITensor& L = (dir == Fromleft ? A1 : A2);
    ITensor& R = (dir == Fromleft ? A2 : A1);

    auto bnd = commonIndex(L,R);
    if(!bnd) return Spectrum();

    if(args.getBool("Verbose",false))
        {
        Print(inds(L));
        }

    ITensor A,B(bnd);
    ITensor D;
    auto spec = svd(L,A,D,B,args);

    L = A;
    R *= (D*B);

    return spec;
    }


MPS& MPS::
position(int i, Args args)
    {
    if(not *this) Error("position: MPS is default constructed");

    if(args.getBool("DoSVDBond",false))
        {
        while(l_orth_lim_ < i-1)
            {
            if(l_orth_lim_ < 0) l_orth_lim_ = 0;
            setBond(l_orth_lim_+1);
            auto WF = operator()(l_orth_lim_+1) * operator()(l_orth_lim_+2);
            auto original_link_tags = tags(linkIndex(*this,l_orth_lim_+1));
            svdBond(l_orth_lim_+1,WF,Fromleft,{args,"LeftTags=",original_link_tags});
            }
        while(r_orth_lim_ > i+1)
            {
            if(r_orth_lim_ > N_+1) r_orth_lim_ = N_+1;
            setBond(r_orth_lim_-2);
            auto WF = operator()(r_orth_lim_-2) * operator()(r_orth_lim_-1);
            auto original_link_tags = tags(linkIndex(*this,r_orth_lim_-2));
            svdBond(r_orth_lim_-2,WF,Fromright,{args,"RightTags=",original_link_tags});
            }
        }
    else //use orthMPS
        {
        while(l_orth_lim_ < i-1)
            {
            if(l_orth_lim_ < 0) l_orth_lim_ = 0;
            setBond(l_orth_lim_+1);

            // Current bond
            auto b = l_orth_lim_+1;

            // Store the original tags for link b so that it can
            // be put back onto the newly introduced link index
            auto original_link_tags = tags(linkIndex(*this,b));
            orthMPS(ref(b),ref(b+1),Fromleft,{args,"LeftTags=",original_link_tags});

            ++l_orth_lim_;
            if(r_orth_lim_ < l_orth_lim_+2) r_orth_lim_ = l_orth_lim_+2;
            }
        while(r_orth_lim_ > i+1)
            {
            if(r_orth_lim_ > N_+1) r_orth_lim_ = N_+1;
            setBond(r_orth_lim_-2);

            // Current bond
            auto b = r_orth_lim_-2;

            // Store the original tags for link b so that it can
            // be put back onto the newly introduced link index
            auto original_link_tags = tags(linkIndex(*this,b));
            args.add("LeftTags=",original_link_tags);
            if(original_link_tags == TagSet("Link,V,0"))
                {
                // Using LeftTags=Link,V will conflict with
                // default right tags of svd (within orthMPS)
                args.add("LeftTags=","Link");
                }
            orthMPS(ref(b),ref(b+1),Fromright,args);

            --r_orth_lim_;
            if(l_orth_lim_ > r_orth_lim_-2) l_orth_lim_ = r_orth_lim_-2;
            }
        }
    return *this;
    }


MPS& MPS::
orthogonalize(Args args)
    {
    if( args.defined("Maxm") )
      {
      if( args.defined("MaxDim") )
        {
        Global::warnDeprecated("Args Maxm and MaxDim are both defined. Maxm is deprecated in favor of MaxDim, MaxDim will be used.");
        }
      else
        {
        Global::warnDeprecated("Arg Maxm is deprecated in favor of MaxDim.");
        args.add("MaxDim",args.getInt("Maxm"));
        }
      }

    if(doWrite()) Error("Cannot call orthogonalize when doWrite()==true");

    auto& psi = *this;
    auto N = N_;

    auto cutoff = args.getReal("Cutoff",1E-13);
    auto dargs = Args{"Cutoff",cutoff};
    auto maxdim_set = args.defined("MaxDim");
    if(maxdim_set) dargs.add("MaxDim",args.getInt("MaxDim"));

    // Truncate blocks of degenerate singular values
    dargs.add("RespectDegenerate",args.getBool("RespectDegenerate",true));

    int rand_plev = 14741;

    //Build environment tensors from the left
    auto E = vector<ITensor>(N+1);

    auto psic = itensor::dag(psi);

    //TODO: use sim(linkInds)
    //That would require changing the requirements of diagPosSemiDef to
    //allow more general ITensors
    psic.replaceLinkInds(itensor::prime(itensor::linkInds(psic),rand_plev));

    E[1] = psi(1)*psic(1); 

    for(int j = 2; j < N; ++j)
        E[j] = E[j-1] * psi(j) * psic(j);
    auto rho = E[N-1] * psi(N) * itensor::prime(psic(N),rand_plev,siteInds(psic,N));

    auto original_tags = tags(linkIndex(psi,N-1));
    auto [U,D] = diagPosSemiDef(rho,{dargs,"Tags=",original_tags});

    //O is partial inner of previous and new MPS
    auto O = U * psi(N) * psi(N-1);
    psi.ref(N) = itensor::dag(U);

    for(int j = N-1; j > 1; --j)
        {
        if(not maxdim_set)
            {
            //Infer maxdim from bond dim of original MPS
            //i.e. upper bound on rank of rho
            auto ci = commonIndex(O,E.at(j-1));
            auto maxdim = (ci) ? dim(ci) : 1l;
            dargs.add("MaxDim",maxdim);
            }
        rho = E.at(j-1) * O * itensor::dag(itensor::prime(O,rand_plev));
        original_tags = tags(linkIndex(psi,j-1));
        std::tie(U,D) = diagPosSemiDef(rho,{dargs,"Tags=",original_tags});
        O *= U;
        O *= psi(j-1);
        psi.ref(j) = itensor::dag(U);
        }
    psi.ref(1) = O;

    l_orth_lim_ = 0;
    r_orth_lim_ = 2;
    return *this;
    }

namespace detail
  {
  // Make an order-2 delta tensor with dense storage
  ITensor
  denseDelta(Index const& i, Index const& j)
      {
      auto del = ITensor(i,j);
      for( auto ii : range1(minDim(del)) )
        del.set(ii,ii,1);
      return del;
      }
  }

template <typename MPSType>
bool
checkOrtho(MPSType const& A,
           int i,
           Direction dir,
           Real threshold)
    {
    auto left = (dir==Fromleft);

    auto Adag = dag(A);
    Adag.replaceLinkInds(prime(linkInds(Adag)));

    auto lout_dag = (left ? leftLinkIndex(Adag,i) : rightLinkIndex(Adag,i));
    auto rho = A(i) * prime(Adag(i),-1,lout_dag);

    auto lin = (left ? rightLinkIndex(A,i) : leftLinkIndex(A,i));
    auto lin_dag = (left ? rightLinkIndex(Adag,i) : leftLinkIndex(Adag,i));
    auto Delta = detail::denseDelta(lin, lin_dag);

    auto Diff = rho - Delta;

    if(norm(Diff) < threshold)
        return true;

    //Print any helpful debugging info here:
    println("checkOrtho: on line ",__LINE__," of mps.h,");
    println("checkOrtho: Tensor at position ",i," failed to be ",left?"left":"right"," ortho.");
    printfln("checkOrtho: norm(Diff) = %E",norm(Diff));
    printfln("checkOrtho: Error threshold set to %E",threshold);
    //-----------------------------

    return false;
    }
template bool checkOrtho<MPS>(MPS const& A, int i, Direction dir, Real threshold);
template bool checkOrtho<MPO>(MPO const& A, int i, Direction dir, Real threshold);

template <typename MPSType>
bool
checkOrtho(MPSType const& A,
           Real threshold)
    {
    for(int i = 1; i <= A.leftLim(); ++i)
    if(!checkOrtho(A,i,Fromleft,threshold))
        {
        std::cout << "checkOrtho: A_[i] not left orthogonal at site i="
                  << i << std::endl;
        return false;
        }

    for(int i = length(A); i >= A.rightLim(); --i)
    if(!checkOrtho(A,i,Fromright,threshold))
        {
        std::cout << "checkOrtho: A_[i] not right orthogonal at site i="
                  << i << std::endl;
        return false;
        }
    return true;
    }
template bool checkOrtho<MPS>(MPS const& A, Real threshold);
template bool checkOrtho<MPO>(MPO const& A, Real threshold);

void
applyGate(ITensor const& gate, 
          MPS & psi,
          Args const& args)
    {
    auto fromleft = args.getBool("Fromleft",true);
    const int c = orthoCenter(psi);
    ITensor AA = psi(c) * psi(c+1) * gate;
    AA.noPrime();
    if(fromleft) psi.svdBond(c,AA,Fromleft,args);
    else         psi.svdBond(c,AA,Fromright,args);
    }

void MPS::
initWrite(Args const& args)
    {
    if(!do_write_)
        {
        std::string write_dir_parent = args.getString("WriteDir","./");
        writedir_ = mkTempDir("psi",write_dir_parent);

        //Write all null tensors to disk immediately because
        //later logic assumes null means written to disk
        for(size_t j = 0; j < A_.size(); ++j)
            {
            if(!A_.at(j)) writeToFile(AFName(j),A_.at(j));
            }

        if(args.getBool("WriteAll",false))
            {
            for(int j = 0; j < int(A_.size()); ++j)
                {
                if(!A_.at(j)) continue;
                writeToFile(AFName(j),A_.at(j));
                if(j < atb_ || j > atb_+1)
                    {
                    A_[j] = ITensor{};
                    }
                }
            }

        do_write_ = true;
        }
    }

void MPS::
copyWriteDir()
    {
    if(do_write_)
        {
        string old_writedir = writedir_;
        string global_write_dir = Args::global().getString("WriteDir","./");
        writedir_ = mkTempDir("psi",global_write_dir);

        string cmdstr = "cp -r " + old_writedir + "/* " + writedir_;
        println("Copying MPS with doWrite()==true. Issuing command: ",cmdstr);
        system(cmdstr.c_str());
        }
    }

void MPS::
cleanupWrite()
    {
    if(do_write_)
        {
        const string cmdstr = "rm -fr " + writedir_;
        system(cmdstr.c_str());
        do_write_ = false;
        }   
    }

void MPS::
swap(MPS& other)
    {
    if(N_ != other.N_)
        Error("Require same system size to swap MPS");
    A_.swap(other.A_);
    std::swap(l_orth_lim_,other.l_orth_lim_);
    std::swap(r_orth_lim_,other.r_orth_lim_);
    std::swap(atb_,other.atb_);
    std::swap(writedir_,other.writedir_);
    std::swap(do_write_,other.do_write_);
    }

InitState::
InitState(SiteSet const& sites)
    : 
    sites_(sites), 
    state_(1+length(sites))
    { 
    for(int n = 1; n <= length(sites_); ++n)
        {
        state_[n] = sites_(n)(1);
        }
    }

InitState::
InitState(SiteSet const& sites, String const& state)
    : 
    sites_(sites), 
    state_(1+length(sites))
    { 
    setAll(state);
    }

InitState& InitState::
set(int i, const String& state)
    { 
    checkRange(i);
    state_.at(i) = sites_(i,state);
    return *this;
    }

InitState& InitState::
setAll(String const& state)
    { 
    for(int n = 1; n <= length(sites_); ++n)
        {
        state_[n] = sites_(n,state);
        }
    return *this;
    }

bool
isComplex(MPS const& psi)
    {
    for(auto j : range1(length(psi)))
        {
        if(itensor::isComplex(psi(j))) return true;
        }
    return false;
    }

bool
isOrtho(MPS const& psi)
    {
    return psi.leftLim()+1 == psi.rightLim()-1;
    }

int
orthoCenter(MPS const& psi)
    {
    if(!isOrtho(psi)) Error("orthogonality center not well defined.");
    return (psi.leftLim() + 1);
    }

int
rightLim(MPS const& x)
    {
    return x.rightLim();
    }

int
leftLim(MPS const& x)
    {
    return x.leftLim();
    }

Real
norm(MPS const& psi)
    {
    if(not isOrtho(psi)) Error("\
MPS must have well-defined ortho center to compute norm; \
call .position(j) or .orthogonalize() to set ortho center");
    return itensor::norm(psi(orthoCenter(psi)));
    }

template<typename MPSType>
Real
averageLinkDim(MPSType const& x)
    {
    Real avgdim = 0.;
    auto N = length(x);
    for( auto b : range1(N-1) )
        {
        avgdim += dim(linkIndex(x,b));
        }
    avgdim /= (N-1);
    return avgdim;
    }
template Real averageLinkDim<MPS>(MPS const& x);
template Real averageLinkDim<MPO>(MPO const& x);

Real
averageM(MPS const& psi)
    {
    Global::warnDeprecated("averageM(MPS) is deprecated in favor of averageLinkDim(MPS)");
    return averageLinkDim(psi);
    }

template<typename MPSType>
int
maxLinkDim(MPSType const& x)
    {
    int maxdim = 0;
    for( auto b : range1(length(x)-1) )
        {
        int mb = dim(linkIndex(x,b));
        maxdim = std::max(mb,maxdim);
        }
    return maxdim;
    }
template int maxLinkDim<MPS>(MPS const& x);
template int maxLinkDim<MPO>(MPO const& x);

int
maxM(MPS const& psi)
    {
    Global::warnDeprecated("maxM(MPS) is deprecated in favor of maxLinkDim(MPS)");
    return maxLinkDim(psi);
    }

Real
normalize(MPS & psi)
    {
    Global::warnDeprecated("normalize(MPS) is deprecated in favor of .normalize()");
    auto nrm = psi.normalize();
    return nrm;
    }

MPS
dag(MPS A)
    {
    A.dag();
    return A;
    }

MPS
setTags(MPS A, TagSet const& ts, IndexSet const& is)
    {
    A.setTags(ts,is);
    return A;
    }

MPS
noTags(MPS A, IndexSet const& is)
    {
    A.noTags(is);
    return A;
    }

MPS
addTags(MPS A, TagSet const& ts, IndexSet const& is)
    {
    A.addTags(ts,is);
    return A;
    }
  
MPS
removeTags(MPS A, TagSet const& ts, IndexSet const& is)
    {
    A.removeTags(ts,is);
    return A;
    }

MPS
replaceTags(MPS A, TagSet const& ts1, TagSet const& ts2, IndexSet const& is)
    {
    A.replaceTags(ts1,ts2,is);
    return A;
    }

MPS
swapTags(MPS A, TagSet const& ts1, TagSet const& ts2, IndexSet const& is)
    {
    A.swapTags(ts1,ts2,is);
    return A;
    }

MPS
prime(MPS A, int plev, IndexSet const& is)
    {
    A.prime(plev,is);
    return A;
    }

MPS
prime(MPS A, IndexSet const& is)
    {
    A.prime(is);
    return A;
    }

MPS
setPrime(MPS A, int plev, IndexSet const& is)
    {
    A.setPrime(plev,is);
    return A;
    }

MPS
mapPrime(MPS A, int plevold, int plevnew, IndexSet const& is)
    {
    A.mapPrime(plevold,plevnew,is);
    return A;
    }

MPS
swapPrime(MPS A, int plevold, int plevnew, IndexSet const& is)
    {
    A.swapPrime(plevold,plevnew,is);
    return A;
    }

MPS
noPrime(MPS A, IndexSet const& is)
    {
    A.noPrime(is);
    return A;
    }

bool
hasSiteInds(MPS const& x, IndexSet const& sites)
    {
    auto N = length(x);
    if( N!=length(sites) ) Error("In hasSiteInds(MPS,IndexSet), lengths of MPS and IndexSet of site indices don't match");
    for( auto n : range1(N) )
      {
      if( !hasIndex(x(n),sites(n)) ) return false;
      }
    return true;
    }

template <typename MPSType>
IndexSet
siteInds(MPSType const& W, int b)
    {
    return uniqueInds(W(b),{W(b-1),W(b+1)});
    }
template IndexSet siteInds<MPS>(MPS const& W, int b);
template IndexSet siteInds<MPO>(MPO const& W, int b);

Index
siteIndex(MPS const& psi, int j)
    {
    return uniqueIndex(psi(j),{psi(j-1),psi(j+1)});
    }

template <typename MPSType>
Index
rightLinkIndex(MPSType const& psi, int i)
    {
    if( i == length(psi) ) return Index();
    return commonIndex(psi(i),psi(i+1));
    }
template Index rightLinkIndex<MPS>(MPS const& W, int i);
template Index rightLinkIndex<MPO>(MPO const& W, int i);

// Note: for ITensors with QNs, this is different
// from rightLinkIndex(psi,i-1) since indices 
// evaluate equal even in arrow directions are different
template <typename MPSType>
Index
leftLinkIndex(MPSType const& psi, int i)
    {
    if( i == 1 ) return Index();
    return commonIndex(psi(i),psi(i-1));
    }
template Index leftLinkIndex<MPS>(MPS const& W, int i);
template Index leftLinkIndex<MPO>(MPO const& W, int i);

// This is a shorthand for rightLinkIndex
template <typename MPSType>
Index
linkIndex(MPSType const& psi, int i)
    {
    return rightLinkIndex(psi,i);
    }
template Index linkIndex<MPS>(MPS const& W, int i);
template Index linkIndex<MPO>(MPO const& W, int i);

template <typename MPSType>
IndexSet
linkInds(MPSType const& x, int i)
    {
    if(i == 1) return IndexSet(rightLinkIndex(x,i));
    else if( i == length(x) ) return IndexSet(leftLinkIndex(x,i));
    return IndexSet(leftLinkIndex(x,i),rightLinkIndex(x,i));
    }
template IndexSet linkInds<MPS>(MPS const& x, int i);
template IndexSet linkInds<MPO>(MPO const& x, int i);

template <typename MPSType>
IndexSet
linkInds(MPSType const& x)
    {
    auto N = length(x);
    auto inds = IndexSetBuilder(N-1);
    for( auto n : range1(N-1) )
      {
      auto s = linkIndex(x,n);
      if(!s) Error("MPS or MPO has missing/null link index");
      inds.nextIndex(std::move(s));
      }
    return inds.build();
    }
template IndexSet linkInds<MPS>(MPS const& x);
template IndexSet linkInds<MPO>(MPO const& x);

template <class MPSType>
bool
hasQNs(MPSType const& x)
    {
    for(auto i : range1(length(x)))
        if(not hasQNs(x(i))) return false;
    return true;
    }
template bool hasQNs<MPS>(MPS const& x);
template bool hasQNs<MPO>(MPO const& x);

IndexSet
siteInds(MPS const& x)
    {
    auto N = length(x);
    auto inds = IndexSetBuilder(N);
    for( auto n : range1(N) )
      {
      auto s = siteIndex(x,n);
      inds.nextIndex(std::move(s));
      }
    return inds.build();
    }

MPS& MPS::
replaceSiteInds(IndexSet const& sites)
    {
    auto& x = *this;
    auto N = itensor::length(x);
    if( itensor::length(sites)!=N ) Error("In replaceSiteInds(MPS,IndexSet), number of site indices not equal to number of MPS tensors");
    auto sx = itensor::siteInds(x);
    if( equals(sx,sites) ) return x;
    for( auto n : range1(N) )
      {
      auto sn = sites(n);
      A_[n].replaceInds({sx(n)},{sn});
      }
    return x;
    }

MPS
replaceSiteInds(MPS x, IndexSet const& sites)
    {
    x.replaceSiteInds(sites);
    return x;
    }

MPS& MPS::
replaceLinkInds(IndexSet const& links)
    {
    auto& x = *this;
    auto N = itensor::length(x);
    if( N==1 ) return x;
    if( itensor::length(links)!=(N-1) ) Error("In replaceLinkInds(MPS,IndexSet), number of link indices input is not equal to the number of links of the MPS");
    auto lx = itensor::linkInds(x);
    if( equals(lx,links) ) return x;
    for( auto n : range1(N-1) )
      {
      if( n == 1 ) A_[n].replaceInds({lx(n)},{links(n)});
      else A_[n].replaceInds({lx(n-1),lx(n)},{links(n-1),links(n)});
      }
    A_[N].replaceInds({lx(N-1)},{links(N-1)});
    return x;
    }

MPS
replaceLinkInds(MPS x, IndexSet const& links)
    {
    x.replaceLinkInds(links);
    return x;
    }

void InitState::
checkRange(int i) const
    {
    if(i > length(sites_) || i < 1) 
        {
        println("i = ",i);
        println("Valid range is 1 to ",length(sites_));
        Error("i out of range");
        }
    }

bool
hasQNs(InitState const& initstate)
    {
    for(auto i : range1(length(initstate.sites())))
        if(not hasQNs(initstate(i))) return false;
    return true;
    }

int 
periodicWrap(int j, int N)
    {
    if(j < 1)
        while(j < 1) j += N;
    else
    if(j > N)
        while(j > N) j -= N;
    return j;
    }


int 
findCenter(MPS const& psi)
    {
    for(int j = 1; j <= length(psi); ++j) 
        {
        auto& A = psi(j);
        if(A.order() == 0) Error("Zero order tensor in MPS");
        bool allSameDir = true;
        auto it = A.inds().begin();
        Arrow dir = (*it).dir();
        for(++it; it != A.inds().end(); ++it)
            {
            if((*it).dir() != dir)
                {
                allSameDir = false;
                break;
                }
            }

        //Found the ortho. center
        if(allSameDir) return j;
        }
    return -1;
    }


std::ostream& 
operator<<(std::ostream& s, MPS const& M)
    {
    s << "\n";
    for(int i = 1; i <= length(M); ++i) 
        {
        s << M(i) << "\n";
        }
    return s;
    }

std::ostream& 
operator<<(std::ostream& s, InitState const& state)
    {
    s << "\n";
    for(int i = 1; i <= length(state.sites()); ++i) 
        {
        s << state(i) << "\n";
        }
    return s;
    }

template <class MPSType>
MPSType
removeQNs(MPSType const& psi)
    {
    int N = length(psi);
    MPSType res;
    res = MPSType(N);
    for(int j = 0; j <= N+1; ++j)
        {
        res.ref(j) = removeQNs(psi(j));
        }
    res.leftLim(psi.leftLim());
    res.rightLim(psi.rightLim());
    return res;
    }
template MPS removeQNs<MPS>(MPS const& psi);
template MPO removeQNs<MPO>(MPO const& psi);

template <class MPSType>
MPSType
sum(MPSType const& L, 
    MPSType const& R, 
    Args const& args)
    {
    auto res = L;
    res.plusEq(R,args);
    return res;
    }
template MPS sum<MPS>(MPS const& L, MPS const& R, Args const& args);
template MPO sum<MPO>(MPO const& L, MPO const& R, Args const& args);

template <class MPSType>
MPSType
sum(std::vector<MPSType> const& terms, 
    Args const& args)
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
template MPS sum<MPS>(std::vector<MPS> const& terms, Args const& args);
template MPO sum<MPO>(std::vector<MPO> const& terms, Args const& args);

Cplx
innerC(MPS const& psi, 
       MPS const& phi)
    {
    auto N = length(psi);
    if(N != length(phi)) Error("inner: mismatched N");

    auto psidag = dag(psi);
    psidag.replaceSiteInds(siteInds(phi));
    psidag.replaceLinkInds(sim(linkInds(psidag)));

    auto L = phi(1) * psidag(1);
    if(N == 1) return eltC(L);
    for(auto i : range1(2,N) ) 
        L = L * phi(i) * psidag(i);
    return eltC(L);
    }

void
inner(MPS const& psi, MPS const& phi, Real& re, Real& im)
    {
    auto z = innerC(psi,phi);
    re = real(z);
    im = imag(z);
    }

Real
inner(MPS const& psi, MPS const& phi) //Re[<psi|phi>]
    {
    if(isComplex(psi) || isComplex(phi)) Error("Cannot use inner(...) with complex MPS/MPO, use innerC(...) instead");
    Real re, im;
    inner(psi,phi,re,im);
    return re;
    }

//
// Deprecated
//


template <class MPSType>
Cplx
overlapC(MPSType const& psi, 
         MPSType const& phi)
    {
    Global::warnDeprecated("overlap is deprecated in favor of inner/trace");
    auto N = length(psi);
    if(N != length(phi)) Error("overlap: mismatched N");

    auto rand_plev = 4351345;

    auto psidag = psi;
    //psidag.dag().primeLinks(rand_plev);
    psidag.dag();
    psidag.replaceLinkInds(prime(linkInds(psidag),rand_plev));

    auto L = phi(1) * psidag(1);
    if(N == 1) return eltC(L);
    for(auto i : range1(2,N) ) 
        L = L * phi(i) * psidag(i);
    return eltC(L);
    }
template Cplx overlapC<MPS>(MPS const& psi, MPS const& phi);
template Cplx overlapC<MPO>(MPO const& psi, MPO const& phi);


template <typename MPSType>
void
overlap(MPSType const& psi, MPSType const& phi, Real& re, Real& im)
    {
    Global::warnDeprecated("overlap is deprecated in favor of inner/trace");
    auto z = overlapC(psi,phi);
    re = real(z);
    im = imag(z);
    }
template void overlap<MPS>(MPS const& psi, MPS const& phi, Real& re, Real& im);
template void overlap<MPO>(MPO const& psi, MPO const& phi, Real& re, Real& im);

template <typename MPSType>
Real
overlap(MPSType const& psi, MPSType const& phi) //Re[<psi|phi>]
    {
    Global::warnDeprecated("overlap is deprecated in favor of inner/trace");
    Real re, im;
    overlap(psi,phi,re,im);
    if(std::fabs(im) > (1E-12 * std::fabs(re)) )
        printfln("Real overlap: WARNING, dropping non-zero imaginary part (=%.5E) of expectation value.",im);
    return re;
    }
template Real overlap<MPS>(MPS const& psi, MPS const& phi);
template Real overlap<MPO>(MPO const& psi, MPO const& phi);

} //namespace itensor
