//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include <map>
#include "itensor/mps/mps.h"
#include "itensor/mps/localop.h"
#include "itensor/util/print_macro.h"

namespace itensor {

using std::map;
using std::istream;
using std::ostream;
using std::cout;
using std::endl;
using std::vector;
using std::find;
using std::pair;
using std::make_pair;
using std::string;

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
MPS(SiteSet const& sites)
    : 
    N_(sites.N()), 
    A_(sites.N()+2), //idmrg may use A_[0] and A[N+1]
    l_orth_lim_(0),
    r_orth_lim_(sites.N()+1),
    sites_(sites), 
    atb_(1),
    writedir_("./"),
    do_write_(false)
    { 
    random_tensors(A_);
    }

MPS::
MPS(InitState const& initState)
    : 
    N_(initState.sites().N()),
    A_(initState.sites().N()+2), //idmrg may use A_[0] and A[N+1]
    l_orth_lim_(0),
    r_orth_lim_(2),
    sites_(initState.sites()), 
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
    sites_(other.sites_),
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
    sites_ = other.sites_;
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

ITensor const& MPS::
A(int i) const
    { 
    if(i < 0) i = N_+i+1;
    setSite(i);
    return A_.at(i); 
    }

ITensor& MPS::
Aref(int i)
    { 
    if(i < 0) i = N_+i+1;
    setSite(i);
    if(i <= l_orth_lim_) l_orth_lim_ = i-1;
    if(i >= r_orth_lim_) r_orth_lim_ = i+1;
    return A_.at(i); 
    }

SiteSet const& MPS::
sites() const 
    { 
    if(not sites_) Error("MPS SiteSet is default-initialized");
    return sites_; 
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
    //Check that tensors read from disk were constructed
    //using the same sites
    auto s1 = findtype(A_.at(1),Site);
    s1.noprime();
    if(sites_ && s1 != Index(sites_(1)))
        {
        Print(A_.at(1).inds());
        Print(s1);
        Print(Index(sites_(1)));
        Error("Tensors read from disk not compatible with SiteSet passed to constructor.");
        }
    itensor::read(s,l_orth_lim_);
    itensor::read(s,r_orth_lim_);
    }

void MPS::
write(std::ostream& s) const
    {
    if(do_write_)
        Error("MPS::write not yet supported if doWrite(true)");

    itensor::write(s,N());
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


void MPS::
new_tensors(std::vector<ITensor>& A_)
    {
    std::vector<Index> a(N_+1);
    for(int i = 1; i <= N_; ++i)
        { 
        a[i] = Index(nameint("a",i)); 
        }
    A_[1] = ITensor(sites()(1),a[1]);
    for(int i = 2; i < N_; i++)
        { 
        A_[i] = ITensor(dag(a[i-1]),sites()(i),a[i]); 
        }
    A_[N_] = ITensor(dag(a[N_-1]),sites()(N_));
    }

void MPS::
random_tensors(std::vector<ITensor>& A_)
    { 
    new_tensors(A_); 
    for(int i = 1; i <= N_; ++i)
        {
        randomize(A_[i]); 
        }
    }

void MPS::
init_tensors(std::vector<ITensor>& A_, InitState const& initState)
    { 
    std::vector<Index> a(N_+1);
    for(auto i : range1(N_)) a[i] = Index(nameint("a",i));

    A_[1] = setElt(IndexVal(initState(1)),a[1](1));
    for(auto i : range(2,N_))
        {
        A_[i] = setElt(dag(a[i-1])(1),IndexVal(initState(i)),a[i](1));
        }
    A_[N_] = setElt(dag(a[N_-1])(1),IndexVal(initState(N_)));
    }


//template <class Tensor>
//void MPS::
//init_tensors(std::vector<IQTensor>& A_, const InitState& initState)
//    {
//    auto qa = std::vector<QN>(N_+1); //qn[i] = qn on i^th bond
//    for(auto i : range1(N_)) qa[0] -= initState(i).qn()*In;
//
//    //Taking OC to be at the leftmost site,
//    //compute the QuantumNumbers of all the Links.
//    for(auto i : range1(N_))
//        {
//        //Taking the divergence to be zero,solve for qa[i]
//        qa[i] = Out*(-qa[i-1]*In - initState(i).qn());
//        }
//
//    auto a = std::vector<IQIndex>(N_+1);
//    for(auto i : range1(N_))
//        { 
//        a[i] = IQIndex(nameint("L",i),Index(nameint("l",i)),qa[i]); 
//        }
//
//    A_[1] = setElt(initState(1),a[1](1));
//    for(auto i : range(2,N_))
//        {
//        A_[i] = setElt(dag(a[i-1])(1),initState(i),a[i](1)); 
//        }
//    A_[N_] = setElt(dag(a[N_-1])(1),initState(N_));
//    }
//template
//void MPS<IQTensor>::
//init_tensors(std::vector<IQTensor>& A_, const InitState& initState);



//template <>
//MPS<IQTensor>& MPS<IQTensor>::operator+=(const MPS<IQTensor>& other)
//    {
//    if(do_write_)
//        Error("operator+= not supported if doWrite(true)");
//
//    primelinks(0,4);
//
//    //Create new link indices
//    vector<IQIndex> nlinks(N);
//    for(int b = 1; b < N_; ++b)
//        {
//        IQIndex l1 = linkInd(*this,b);
//        IQIndex l2 = linkInd(other,b);
//        vector<IndexQN> iq(l1.indices());
//        iq.insert(iq.begin(),l2.indices().begin(),l2.indices().end());
//        nlinks.at(b) = IQIndex(l2,iq);
//        }
//    //Create new A tensors
//    vector<IQTensor> nA(N+1);
//    nA[1] = IQTensor(si(1),nlinks[1]);
//    for(int j = 2; j < N_; ++j)
//        nA[j] = IQTensor(dag(nlinks[j-1]),si(j),nlinks[j]);
//    nA[N] = IQTensor(dag(nlinks[N-1]),si(N));
//
//    for(int j = 1; j <= N_; ++j)
//        {
//        Foreach(const ITensor& t, A(j).blocks())
//            { nA[j].insert(t); }
//        Foreach(const ITensor& t, other.A(j).blocks())
//            { nA[j].insert(t); }
//        }
//
//    A.swap(nA);
//
//    orthogonalize();
//
//    return *this;
//    }
//
//template <class Tensor>
//MPS& MPS::operator+=(const MPS& other)
//    {
//    if(do_write_)
//        Error("operator+= not supported if doWrite(true)");
//
//    primelinks(0,4);
//
//    vector<Tensor> first(N), second(N);
//    for(int i = 1; i < N_; ++i)
//        {
//        Index l1 = rightLinkInd(*this,i);
//        Index l2 = rightLinkInd(other,i);
//        Index r(l1);
//        plussers(l1,l2,r,first[i],second[i]);
//        }
//
//    Anc(1) = A(1) * first[1] + other.A(1) * second[1];
//    for(int i = 2; i < N_; ++i)
//        {
//        Anc(i) = dag(first[i-1]) * A(i) * first[i] 
//                  + dag(second[i-1]) * other.A(i) * second[i];
//        }
//    Anc(N) = dag(first[N-1]) * A(N) + dag(second[N-1]) * other.A(N);
//
//    noprimelink();
//
//    orthogonalize();
//
//    return *this;
//    }
//template
//MPS<ITensor>& MPS<ITensor>::operator+=(const MPS<ITensor>& other);

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



//
//MPS Index Methods
//

void MPS::
mapprime(int oldp, int newp, IndexType type)
    { 
    if(do_write_)
        Error("mapprime not supported if doWrite(true)");
    for(int i = 1; i <= N_; ++i) 
        A_[i].mapprime(oldp,newp,type); 
    }

void MPS::
primelinks(int oldp, int newp)
    { 
    if(do_write_)
        Error("primelinks not supported if doWrite(true)");
    for(int i = 1; i <= N_; ++i) 
        A_[i].mapprime(oldp,newp,Link); 
    }

void MPS::
noprimelink()
    { 
    if(do_write_)
        Error("noprimelink not supported if doWrite(true)");
    for(int i = 1; i <= N_; ++i) 
        A_[i].noprime(Link); 
    }

Spectrum MPS::
svdBond(int b, ITensor const& AA, Direction dir, Args const& args)
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
orthMPS(ITensor& A1, ITensor& A2, Direction dir, const Args& args)
    {
    ITensor& L = (dir == Fromleft ? A1 : A2);
    ITensor& R = (dir == Fromleft ? A2 : A1);

    auto bnd = commonIndex(L,R,Link);
    if(!bnd) return Spectrum();

    if(args.getBool("Verbose",false))
        {
        Print(L.inds());
        }

    ITensor A,B(bnd);
    ITensor D;
    auto spec = svd(L,A,D,B,args);

    L = A;
    R *= (D*B);

    return spec;
    }


void MPS::
position(int i, Args const& args)
    {
    if(not *this) Error("position: MPS is default constructed");

    if(args.getBool("DoSVDBond",false))
        {
        while(l_orth_lim_ < i-1)
            {
            if(l_orth_lim_ < 0) l_orth_lim_ = 0;
            setBond(l_orth_lim_+1);
            auto WF = A(l_orth_lim_+1) * A(l_orth_lim_+2);
            svdBond(l_orth_lim_+1,WF,Fromleft,args);
            }
        while(r_orth_lim_ > i+1)
            {
            if(r_orth_lim_ > N_+1) r_orth_lim_ = N_+1;
            setBond(r_orth_lim_-2);
            auto WF = A(r_orth_lim_-2) * A(r_orth_lim_-1);
            svdBond(r_orth_lim_-2,WF,Fromright,args);
            }
        }
    else //use orthMPS
        {
        while(l_orth_lim_ < i-1)
            {
            if(l_orth_lim_ < 0) l_orth_lim_ = 0;
            setBond(l_orth_lim_+1);
            orthMPS(Aref(l_orth_lim_+1),Aref(l_orth_lim_+2),Fromleft,args);
            ++l_orth_lim_;
            if(r_orth_lim_ < l_orth_lim_+2) r_orth_lim_ = l_orth_lim_+2;
            }
        while(r_orth_lim_ > i+1)
            {
            if(r_orth_lim_ > N_+1) r_orth_lim_ = N_+1;
            setBond(r_orth_lim_-2);
            orthMPS(Aref(r_orth_lim_-2),Aref(r_orth_lim_-1),Fromright,args);
            --r_orth_lim_;
            if(l_orth_lim_ > r_orth_lim_-2) l_orth_lim_ = r_orth_lim_-2;
            }
        }
    }

int MPS::
orthoCenter() const 
    { 
    if(!itensor::isOrtho(*this)) Error("orthogonality center not well defined.");
    return (leftLim() + 1);
    }

void MPS::
orthogonalize(Args const& args)
    {
    if(doWrite()) Error("Cannot call orthogonalize when doWrite()==true");

    auto cutoff = args.getReal("Cutoff",1E-13);
    auto dargs = Args{"Cutoff",cutoff};
    auto maxm_set = args.defined("Maxm");
    if(maxm_set) dargs.add("Maxm",args.getInt("Maxm"));

    int plev = 14741;

    //Build environment tensors from the left
    auto E = vector<ITensor>(N_+1);
    auto ci = commonIndex(A_.at(1),A_.at(2),Link);
    E.at(1) = A_.at(1)*dag(prime(A_.at(1),ci,plev));
    for(int j = 2; j < N_; ++j)
        {
        E.at(j) = E.at(j-1) * A_.at(j) * dag(prime(A_.at(j),Link,plev));
        }

    auto rho = E.at(N_-1) * A_.at(N_) * dag(prime(A_.at(N_),plev));
    ITensor U,D;
    diagHermitian(rho,U,D,{dargs,"IndexType=",Link});

    //O is partial overlap of previous and new MPS
    auto O = U * A_.at(N_) * A_.at(N_-1);
    A_.at(N_) = dag(U);

    for(int j = N_-1; j > 1; --j)
        {
        if(not maxm_set)
            {
            //Infer maxm from bond dim of original MPS
            //i.e. upper bound on rank of rho
            auto ci = commonIndex(O,E.at(j-1));
            auto maxm = (ci) ? ci.m() : 1l;
            dargs.add("Maxm",maxm);
            }
        rho = E.at(j-1) * O * dag(prime(O,plev));
        auto spec = diagHermitian(rho,U,D,{dargs,"IndexType=",Link});
        O *= U;
        O *= A_.at(j-1);
        A_.at(j) = dag(U);
        }
    A_.at(1) = O;

    l_orth_lim_ = 0;
    r_orth_lim_ = 2;
    }

//Methods for use internally by checkOrtho
ITensor
makeKroneckerDelta(Index const& i, int plev)
    {
    return delta(i,prime(i,plev));
    }
//IQTensor
//makeKroneckerDelta(const IQIndex& I, int plev)
//    {
//    IQTensor D(I,prime(I,plev));
//
//    for(int j = 1; j <= I.nindex(); ++j)
//        {
//        D += makeKroneckerDelta(I.index(j),plev);
//        }
//    return D;
//    }

//template <class Tensor>
//bool MPS::
//checkOrtho(int i, bool left) const
//    {
//    setSite(i);
//    IndexT link = (left ? rightLinkInd(*this,i) : leftLinkInd(*this,i));
//
//    Tensor rho = A(i) * dag(prime(A(i),link,4));
//
//    Tensor Delta = makeKroneckerDelta(link,4);
//
//    Tensor Diff = rho - Delta;
//
//    const
//    Real threshold = 1E-13;
//    //cout << format("i = %d, Diff.norm() = %.4E")
//    //        % i
//    //        % Diff.norm()
//    //        << endl;
//    if(Diff.norm() < threshold) 
//        {
//        return true;
//        }
//
//    //Print any helpful debugging info here:
//    cout << "checkOrtho: on line " << __LINE__ 
//         << " of mps.h," << endl;
//    cout << "checkOrtho: Tensor at position " << i 
//         << " failed to be " << (left ? "left" : "right") 
//         << " ortho." << endl;
//    cout << "checkOrtho: Diff.norm() = " << format("%E") 
//         % Diff.norm() << endl;
//    cout << "checkOrtho: Error threshold set to " 
//              << format("%E") % threshold << endl;
//    //-----------------------------
//
//    return false;
//    }
//template
//bool MPS<ITensor>::checkOrtho(int i, bool left) const;
//template
//bool MPS<IQTensor>::checkOrtho(int i, bool left) const;

//template <class Tensor>
//bool MPS::
//checkOrtho() const
//    {
//    for(int i = 1; i <= l_orth_lim_; ++i)
//    if(!checkLeftOrtho(i))
//        {
//        cout << "checkOrtho: A_[i] not left orthogonal at site i=" 
//                  << i << endl;
//        return false;
//        }
//
//    for(int i = N(); i >= r_orth_lim_; --i)
//    if(!checkRightOrtho(i))
//        {
//        cout << "checkOrtho: A_[i] not right orthogonal at site i=" 
//                  << i << endl;
//        return false;
//        }
//    return true;
//    }
//template
//bool MPS<ITensor>::checkOrtho() const;
//template
//bool MPS<IQTensor>::checkOrtho() const;

bool
checkOrtho(MPS const& psi,
           int i, 
           bool left)
    {
    Index link = (left ? rightLinkInd(psi,i) : leftLinkInd(psi,i));
    ITensor rho = psi.A(i) * dag(prime(psi.A(i),link,4));
    ITensor Delta = delta(link, prime(link, 4));
    ITensor Diff = rho - Delta;

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

bool
checkOrtho(MPS const& psi)
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


//template <class Tensor>
//void MPS::
//applygate(const Tensor& gate, const Args& args)
//    {
//    setBond(l_orth_lim_+1);
//    Tensor AA = A_.at(l_orth_lim_+1) * A_.at(l_orth_lim_+2) * gate;
//    AA.noprime();
//    svdBond(l_orth_lim_+1,AA,Fromleft,args);
//    }
//template
//void MPS<ITensor>::applygate(const ITensor& gate,const Args& args);
//template
//void MPS<IQTensor>::applygate(const IQTensor& gate,const Args& args);

//template <class Tensor>
//void MPS::
//applygate(const BondGate<Tensor>& gate, 
//          const Args& args)
//    {
//    const int gate_b = std::min(gate.i(),gate.j());
//    setBond(gate_b);
//    Tensor AA = A_.at(gate_b) * A_.at(gate_b+1) * Tensor(gate);
//    AA.noprime();
//    svdBond(gate_b,AA,Fromleft,args);
//    }
//template
//void MPS<ITensor>::applygate(const BondGate<ITensor>& gate,const Args& args);
//template
//void MPS<IQTensor>::applygate(const BondGate<IQTensor>& gate,const Args& args);

void
applyGate(ITensor const& gate, 
          MPS & psi,
          Args const& args)
    {
    auto fromleft = args.getBool("Fromleft",true);
    const int c = psi.orthoCenter();
    ITensor AA = psi.A(c) * psi.A(c+1) * gate;
    AA.noprime();
    if(fromleft) psi.svdBond(c,AA,Fromleft,args);
    else         psi.svdBond(c,AA,Fromright,args);
    }

Real MPS::
norm() const 
    { 
    return itensor::norm(*this);
    }

Real MPS::
normalize()
    {
    return itensor::normalize(*this);
    }

bool MPS::
isComplex() const
    { 
    for(auto j : range1(N_))
        {
        if(itensor::isComplex(A_[j])) return true;
        }
    return false;
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

        writeToFile(writedir_+"/sites",sites_);

        do_write_ = true;
        }
    }

void MPS::
copyWriteDir()
    {
    if(do_write_)
        {
        string old_writedir = writedir_;
        string global_write_dir = Global::args().getString("WriteDir","./");
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
    std::swap(sites_,other.sites_);
    std::swap(atb_,other.atb_);
    std::swap(writedir_,other.writedir_);
    std::swap(do_write_,other.do_write_);
    }

InitState::
InitState(SiteSet const& sites)
    : 
    sites_(sites), 
    state_(1+sites.N())
    { 
    for(int n = 1; n <= sites_.N(); ++n)
        {
        state_[n] = sites_(n)(1);
        }
    }

InitState::
InitState(SiteSet const& sites, String const& state)
    : 
    sites_(sites), 
    state_(1+sites.N())
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
    for(int n = 1; n <= sites_.N(); ++n)
        {
        state_[n] = sites_(n,state);
        }
    return *this;
    }

void InitState::
checkRange(int i) const
    {
    if(i > sites_.N() || i < 1) 
        {
        println("i = ",i);
        println("Valid range is 1 to ",sites_.N());
        Error("i out of range");
        }
    }

////Auxilary method for convertToIQ
//long 
//collapseCols(Vector const& Diag, Matrix& M)
//    {
//    long nr = Diag.size(), 
//         nc = long(sumels(Diag));
//    assert(nr != 0);
//    if(nc == 0) return nc;
//    M = Matrix(nr,nc);
//    long c = 0;
//    for(long r = 1; r <= nr; ++r)
//        if(Diag(r) == 1) { M(r,++c) = 1; }
//    return nc;
//    }

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

//void 
//convertToIQ(const SiteSet& sites, const vector<ITensor>& A, 
//            vector<IQTensor>& qA, QN totalq, Real cut)
//    {
//    const int N = sites.N();
//    qA.resize(A.size());
//    const bool is_mpo = hasindex(A[1],sites.siP(1));
//    const int Dim = sites.si(1).m();
//    if(sites.si(2).m() != Dim)
//        Error("convertToIQ assumes uniform site dimension");
//    const int PDim = (is_mpo ? Dim : 1);
//
//    // If MPO, set all tensors to identity ops initially
//    if(is_mpo)
//        {
//        for(int j = 1; j <= N; ++j)
//            qA.at(j) = sites.op("Id",j);
//        }
//
//    const int fullrank = (is_mpo ? 4 : 3);
//    int start = 1, end = N;
//
//    for(int j = 1; j <= N; ++j)
//        if(A[j].r() == fullrank)
//            if(A.at(periodicWrap(j-1,N)).r() < fullrank) 
//                {
//                start = periodicWrap(j-1,N);
//                //cout << "Got start at " << start << "\n";
//                break;
//                }
//
//    for(int j = 1; j <= N; ++j)
//        if(A[j].r() == fullrank)
//            if(A.at(periodicWrap(j+1,N)).r() < fullrank) 
//                {
//                end = periodicWrap(j+1,N);
//                //cout << "Got end at " << end << "\n";
//                break;
//                }
//
//    //cout << "Converting to IQ with (start, end) = " << start SP end << endl;
//
//    vector<IQIndex> linkind(N+1);
//
//    map<QN,Vector> qD; //Diags of compressor matrices by QN
//
//    using qt_vt = map<QN,vector<ITensor> >::value_type;
//    map<QN,vector<ITensor> > qt; //ITensor blocks by QN
//
//    using qC_vt = map<QN,ITensor>::value_type;
//    map<QN,ITensor> qC; //Compressor ITensors by QN
//
//    ITensor block;
//    vector<ITensor> nblock;
//    vector<IndexQN> iq;
//
//    QN q;
//
//    qC[totalq] = ITensor(); //Represents Virtual index
//    //First value of prev_q below set to totalq
//
//    const int show_s = 0;
//
//    Index bond, prev_bond;
//    int Send = (end < start ? N+end : end); 
//    for(int S = start; S <= Send; ++S)
//        {
//        int s = periodicWrap(S,N);
//        int sprev = periodicWrap(S-1,N);
//        int snext = periodicWrap(S+1,N);
//
//        qD.clear(); 
//        qt.clear();
//
//        if(S > start) prev_bond = commonIndex(A[sprev],A[s],Link);
//        if(S < Send) bond = commonIndex(A[s],A[snext],Link);
//
//        if(s == show_s) { PrintData(A[s]); }
//
//        for(const qC_vt& x : qC) 
//        for(int n = 1; n <= Dim;  ++n)
//        for(int u = 1; u <= PDim; ++u)
//            {
//            //Each compressor represents a particular
//            //QN channel given by prev_q
//            const QN& prev_q = x.first; 
//            //Alias previous compressor ITensor to comp
//            const ITensor& comp = x.second; 
//
//            q = (is_mpo ? prev_q+sites.si(s).qn(n)-sites.si(s).qn(u) 
//                        : prev_q-sites.si(s).qn(n));
//
//            //For the last site, only keep blocks 
//            //compatible with specified totalq i.e. q=0 here
//            if(S == Send && q != QN()) continue;
//
//            //Set Site indices of A[s] and its previous Link Index
//            block = A[s];
//            if(S != start) block *= dag(comp);
//            block *= Index(sites.si(s))(n);
//            if(is_mpo) block *= Index(sites.siP(s))(u);
//
//            //Initialize D Vector (D records which values of
//            //the right Link Index to keep for the current QN q)
//            auto count = qD.count(q);
//            Vector& D = qD[q];
//            if(count == 0) 
//                { 
//                D.resize(bond.m()); 
//                for(auto& el : D) el = 0; 
//                }
//
//            if(s == show_s)
//                {
//                println("For n = ",n);
//                printfln("Got a block with norm %.10f",norm(block));
//                println("bond.m() = ",bond.m());
//                PrintData(block);
//                if(s != 1) PrintData(comp);
//                }
//
//            bool keep_block = false;
//            if(S == Send) 
//                { keep_block = true; }
//            else
//                {
//                if(bond.m() == 1 && norm(block) != 0) 
//                    { 
//                    for(auto& el : D) el = 1; 
//                    keep_block = true; 
//                    }
//                else
//                    {
//                    ITensor summed_block;
//                    if(S==start) 
//                        { summed_block = block; }
//                    else
//                        {
//                        //Here we sum over the previous link index
//                        //which is already ok, analyze the one to the right
//                        assert(comp.r()==2);
//                        auto ci = comp.inds().begin();
//                        const Index& new_ind = (*ci==prev_bond ? *(ci+1) : *ci);
//                        summed_block = diag(1,new_ind) * block;
//                        }
//                    //summed_block.print("summed_block");
//
//                    Real rel_cut = -1;
//                    const ITensor& sb = summed_block;
//                    for(int j = 1; j <= bond.m(); ++j)
//                        { rel_cut = std::max(std::fabs(sb.real(bond(j))),rel_cut); }
//                    assert(rel_cut >= 0);
//                    //Real rel_cut = summed_block.norm()/summed_block.vecSize();
//                    rel_cut *= cut;
//                    //cout << "rel_cut == " << rel_cut << "\n";
//
//                    if(rel_cut > 0)
//                    for(int j = 1; j <= bond.m(); ++j)
//                        {
//                        if(std::fabs(sb.real(bond(j))) > rel_cut) 
//                            { 
//                            D(j) = 1; 
//                            keep_block = true; 
//                            }
//                        }
//                    }
//                } //else (S != Send)
//
//            if(keep_block)
//                {
//                qD[q] = D;
//
//                IndexSet newinds(block.inds());
//                if(is_mpo) 
//                    {
//                    newinds.addindex(dag(sites.si(s)(n).indexqn()));
//                    newinds.addindex(sites.siP(s)(u).indexqn());
//                    }
//                else 
//                    { 
//                    newinds.addindex(sites.si(s)(n).indexqn()); 
//                    }
//
//                if(s==show_s)
//                    {
//                    PrintData(block);
//                    cout << "D = " << D << "\n";
//                    }
//
//                qt[q].push_back(ITensor(std::move(newinds),std::move(block.store())));
//
//                }
//            }
//
//        qC.clear();
//
//        for(const qt_vt& x : qt)
//            {
//            const vector<ITensor>& blks = x.second;
//            if(blks.size() != 0)
//                {
//                q = x.first; 
//                if(S == Send) 
//                    { for(const ITensor& t : blks) nblock.push_back(t); }
//                else
//                    {
//                    Matrix M; 
//                    auto mm = collapseCols(qD[q],M);
//                    if(s==show_s)
//                        {
//                        println("Adding block, mm = ",mm);
//                        Print(q);
//                        cout << "qD[q] = " << qD[q] << "\n";
//                        cout << "M = \n" << M << "\n";
//                        int count = 0;
//                        for(const ITensor& t : blks) 
//                            {
//                            printfln("t%02d",++count," ",t);
//                            }
//                        }
//                    string qname = format("ql%d(%+d:%d)",s,q.sz(),q.Nf());
//                    Index qbond(qname,mm);
//                    auto compressor = matrixTensor(std::move(M),bond,qbond);
//                    for(const ITensor& t : blks) nblock.push_back(t * compressor);
//                    iq.push_back(IndexQN(qbond,q));
//                    qC[q] = compressor;
//                    }
//                }
//            }
//
//        if(S != Send) 
//            { 
//            if(iq.empty()) 
//                {
//                cout << "At site " << s << "\n";
//                Error("convertToIQ: no compatible QNs to put into Link.");
//                }
//            linkind[s] = IQIndex(nameint("qL",s),std::move(iq)); 
//            }
//        if(S == start)
//            {
//            qA.at(s) = (is_mpo ? IQTensor(dag(sites.si(s)),sites.siP(s),linkind.at(s)) 
//                            : IQTensor(sites.si(s),linkind[s]));
//            }
//        else 
//        if(S == Send)
//            {
//            qA.at(s) = (is_mpo ? IQTensor(dag(linkind[sprev]),dag(sites.si(s)),sites.siP(s)) 
//                            : IQTensor(dag(linkind[sprev]),sites.si(s)));
//            }
//        else
//            {
//            qA.at(s) = (is_mpo ? IQTensor(dag(linkind[sprev]),dag(sites.si(s)),sites.siP(s),linkind[s]) 
//                            : IQTensor(dag(linkind[sprev]),sites.si(s),linkind[s]));
//            }
//
//        for(const ITensor& nb : nblock) 
//            { qA.at(s) += nb; } 
//        nblock.clear();
//
//        if(s==show_s)
//            {
//            printfln("qA[%d]",s,qA[s]);
//            Error("Stopping");
//            }
//
//        } //for loop over s
//
//    } //void convertToIQ

/*
template <class Tensor> 
template <class IQMPSType> 
void MPS::convertToIQ(IQMPSType& iqpsi, QN totalq, Real cut) const
{
    assert(sites_ != 0);
    const SiteSet& sst = *sites_;

    iqpsi = IQMPSType(sst,maxm,cutoff);

    if(!A_[1].hasindex(si(1))) Error("convertToIQ: incorrect primelevel for conversion");
    bool is_mpo = A_[1].hasindex(prime(si(1)));
    const int Dim = si(1).m();
    const int PDim = (is_mpo ? Dim : 1);

    vector<IQIndex> linkind(N);

    using qD_vt = map<QN,Vector>::value_type;
    map<QN,Vector> qD; //Diags of compressor matrices by QN
    using qt_vt = map<QN,vector<ITensor> >::value_type;
    map<QN,vector<ITensor> > qt; //ITensor blocks by QN
    using qC_vt = map<QN,ITensor>::value_type;
    map<QN,ITensor> qC; //Compressor ITensors by QN
    ITensor block;
    vector<ITensor> nblock;
    vector<IndexQN> iq;

    QN q;

    qC[totalq] = ITensor(); //Represents Virtual index
    //First value of prev_q below set to totalq

    const int show_s = 0;

    Index bond, prev_bond;
    for(int s = 1; s <= N; ++s)
    {

        qD.clear(); qt.clear();
        if(s > 1) prev_bond = linkInd(*this,s-1); 
        if(s < N) bond = linkInd(*this,s);

        if(s == show_s) 
        {
            PrintData(A_[s]);
        }

        Foreach(const qC_vt& x, qC) {
        const QN& prev_q = x.first; const ITensor& comp = x.second; 
        for(int n = 1; n <= Dim;  ++n)
        for(int u = 1; u <= PDim; ++u)
        {
            q = (is_mpo ? prev_q+si(s).qn(n)-si(s).qn(u) : prev_q-si(s).qn(n));

            //For the last site, only keep blocks 
            //compatible with specified totalq i.e. q=0 here
            if(s == N && q != QN()) continue;

            //Set Site indices of A_[s] and its previous Link Index
            block = A_[s];
            if(s != 1) block *= dag(comp);
            block *= si(s)(n);
            if(is_mpo) block *= siP(s)(u);

            //Initialize D Vector (D records which values of
            //the right Link Index to keep for the current QN q)
            int count = qD.count(q);
            Vector& D = qD[q];
            if(count == 0) { D.ReDimension(bond.m()); D = 0; }

            if(s == show_s)
            {
                cout << format("For n = %d\n")%n;
                cout << format("Got a block with norm %.10f\n")%block.norm();
                cout << format("bond.m() = %d\n")%bond.m();
                PrintData(block);
                if(s != 1) PrintData(comp);
            }

            bool keep_block = false;
            if(s == N) keep_block = true;
            else
            {
                if(bond.m() == 1 && block.norm() != 0) { D = 1; keep_block = true; }
                else
                {
                    ITensor summed_block;
                    if(s==1) summed_block = block;
                    else
                    {
                        //Here we sum over the previous link index
                        //which is already ok, analyze the one to the right
                        assert(comp.r()==2);
                        Index new_ind = (comp.index(1)==prev_bond ? comp.index(2) : comp.index(1));
                        summed_block = ITensor(new_ind,1) * block;
                    }
                    //cout << format("s = %d, bond=")%s << bond << "\n";
                    //summed_block.print("summed_block");

                    Real rel_cut = -1;
                    for(int j = 1; j <= bond.m(); ++j)
                    { rel_cut = std::max(std::fabs(summed_block.val1(j)),rel_cut); }
                    assert(rel_cut >= 0);
                    //Real rel_cut = summed_block.norm()/summed_block.vecSize();
                    rel_cut *= cut;
                    //cout << "rel_cut == " << rel_cut << "\n";

                    if(rel_cut > 0)
                    for(int j = 1; j <= bond.m(); ++j)
                    if(std::fabs(summed_block.val1(j)) > rel_cut) 
                    { D(j) = 1; keep_block = true; }
                }
            } //else (s != N)

            //if(!keep_block && q == totalq) { D(1) = 1; keep_block = true; }

            if(keep_block)
            {
                qD[q] = D;

                if(is_mpo) 
                {
                block.addindex(dag(si(s)(n).index()));
                block.addindex(siP(s)(u).index());
                }
                else { block.addindex(si(s)(n).index()); }

                qt[q].push_back(block);

                if(s==show_s)
                {
                block.print("Kept block",ShowData);
                cout << "D = " << D << "\n";
                }
            }
        }}

        qC.clear();

        Foreach(const qt_vt& x, qt)
            {
            const vector<ITensor>& blks = x.second;
            if(blks.size() != 0)
                {
                q = x.first; 
                if(s == N) 
                    { Foreach(const ITensor& t, blks) nblock.push_back(t); }
                else
                    {
                    Matrix M; 
                    auto mm = collapseCols(qD[q],M);
                    if(s==show_s)
                        {
                        cout << format("Adding block, mm = %d\n")%mm;
                        q.print("q");
                        cout << "qD[q] = " << qD[q] << "\n";
                        cout << "M = \n" << M << "\n";
                        int count = 0;
                        Foreach(const ITensor& t, blks) 
                        t.print((format("t%02d")%(++count)).str(),ShowData);
                        }
                    //string qname = (format("ql%d(%+d:%d:%s)")%s%q.sz()%q.Nf()%(q.Nfp() == 0 ? "+" : "-")).str();
                    string qname = (format("ql%d(%+d:%d)")%s%q.sz()%q.Nf()).str();
                    Index qbond(qname,mm);
                    ITensor compressor(bond,qbond,M);
                    Foreach(const ITensor& t, blks) nblock.push_back(t * compressor);
                    iq.push_back(IndexQN(qbond,q));
                    qC[q] = compressor;
                    }
                }
            }

        if(s != N) 
        { 
            if(iq.empty()) 
            {
                cout << "At site " << s << "\n";
                Error("convertToIQ: no compatible QNs to put into Link.");
            }
            linkind[s] = IQIndex(nameint("qL",s),iq); iq.clear(); 
        }
        if(s == 1)
        {
            iqpsi.Anc(s) = (is_mpo ? IQTensor(dag(si(s)),siP(s),linkind[s]) : IQTensor(si(s),linkind[s]));
        }
        else if(s == N)
        {
            iqpsi.Anc(s) = (is_mpo ? IQTensor(dag(linkind[s-1]),dag(si(s)),siP(s)) 
                                    : IQTensor(dag(linkind[s-1]),si(s)));
        }
        else
        {
            iqpsi.Anc(s) = (is_mpo ? IQTensor(dag(linkind[s-1]),dag(si(s)),siP(s),linkind[s]) 
                                    : IQTensor(dag(linkind[s-1]),si(s),linkind[s]));
        }

        Foreach(const ITensor& nb, nblock) { iqpsi.Anc(s) += nb; } nblock.clear();

        if(0) //try to get this working ideally
        if(!is_mpo && s > 1) 
        {
            IQTensor AA = iqpsi.bondTensor(s-1);
            iqpsi.doSVD(s-1,AA,Fromleft);
        }

        if(s==show_s)
        {
        iqpsi.A(s).print((format("qA[%d]")%s).str(),ShowData);
        Error("Stopping");
        }

    } //for loop over s

    assert(checkQNs(iqpsi));

} //void convertToIQ(IQMPSType& iqpsi) const
*/

int 
findCenter(MPS const& psi)
    {
    for(int j = 1; j <= psi.N(); ++j) 
        {
        auto& A = psi.A(j);
        if(A.r() == 0) Error("Zero rank tensor in MPS");
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
    for(int i = 1; i <= M.N(); ++i) 
        {
        s << M.A(i) << "\n";
        }
    return s;
    }

std::ostream& 
operator<<(std::ostream& s, InitState const& state)
    {
    s << "\n";
    for(int i = 1; i <= state.sites().N(); ++i) 
        {
        s << state(i) << "\n";
        }
    return s;
    }

//MPS
//toMPS(IQMPS const& psi)
//    {
//    int N = psi.N();
//    MPS res;
//    if(psi.sites()) res = MPS(psi.sites());
//    else            res = MPS(N);
//    for(int j = 0; j <= N+1; ++j)
//        {
//        res.Aref(j) = ITensor(psi.A(j));
//        }
//    res.leftLim(psi.leftLim());
//    res.rightLim(psi.rightLim());
//    return res;
//    }

MPS
sum(MPS const& L, 
    MPS const& R, 
    Args const& args)
    {
    auto res = L;
    res.plusEq(R,args);
    return res;
    }

MPS
sum(std::vector<MPS> const& terms, 
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
        std::vector<MPS> newterms(nsize); 
        for(decltype(Nt) n = 0, np = 0; n < Nt-1; n += 2, ++np)
            {
            newterms.at(np) = sum(terms.at(n),terms.at(n+1),args);
            }
        if(Nt%2 == 1) newterms.at(nsize-1) = terms.back();

        //Recursively call sum again
        return sum(newterms,args);
        }
    return MPS();
    }

Cplx
overlapC(MPS const& psi, 
         MPS const& phi)
    {
    auto N = psi.N();
    if(N != phi.N()) Error("overlap: mismatched N");

    auto l1 = linkInd(psi,1);
    auto L = phi.A(1);
    if(l1) L *= dag(prime(psi.A(1),l1)); 
    else   L *= dag(psi.A(1));

    if(N == 1) return L.cplx();

    for(decltype(N) i = 2; i < N; ++i) 
        { 
        L = L * phi.A(i) * dag(prime(psi.A(i),Link)); 
        }
    L = L * phi.A(N);

    auto lNm = linkInd(psi,N-1);
    if(lNm) return (dag(prime(psi.A(N),lNm))*L).cplx();
    return (dag(psi.A(N))*L).cplx();
    }

} //namespace itensor
