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
// class MPSt
//

//
// Constructors
//

template <class Tensor>
MPSt<Tensor>::
MPSt() 
    : 
    N_(0), 
    sites_(nullptr),
    atb_(1),
    writedir_("./"),
    do_write_(false)
    { }
template MPSt<ITensor>::
MPSt();
template MPSt<IQTensor>::
MPSt();

template <class Tensor>
MPSt<Tensor>::
MPSt(const SiteSet& sites)
    : 
    N_(sites.N()), 
    A_(sites.N()+2), //idmrg may use A_[0] and A[N+1]
    l_orth_lim_(0),
    r_orth_lim_(sites.N()+1),
    sites_(&sites), 
    atb_(1),
    writedir_("./"),
    do_write_(false)
    { 
    random_tensors(A_);
    }
template MPSt<ITensor>::
MPSt(const SiteSet& sites);
template MPSt<IQTensor>::
MPSt(const SiteSet& sites);

template <class Tensor>
MPSt<Tensor>::
MPSt(const InitState& initState)
    : 
    N_(initState.sites().N()),
    A_(initState.sites().N()+2), //idmrg may use A_[0] and A[N+1]
    l_orth_lim_(0),
    r_orth_lim_(2),
    sites_(&(initState.sites())), 
    atb_(1),
    writedir_("./"),
    do_write_(false)
    { 
    init_tensors(A_,initState);
    }
template MPSt<ITensor>::
MPSt(const InitState& initState);
template MPSt<IQTensor>::
MPSt(const InitState& initState);

//template <class Tensor>
//MPSt<Tensor>::
//MPSt(const SiteSet& sites, std::istream& s)
//    : 
//    N_(sites.N()), 
//    A_(sites.N()+2), //idmrg may use A_[0] and A[N+1]
//    sites_(&sites),
//    atb_(1),
//    writedir_("."),
//    do_write_(false)
//    { 
//    read(s); 
//    }
//template MPSt<ITensor>::
//MPSt(const SiteSet& sites, std::istream& s);
//template MPSt<IQTensor>::
//MPSt(const SiteSet& sites, std::istream& s);

template <class Tensor>
MPSt<Tensor>::
MPSt(const MPSt& other)
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
template MPSt<ITensor>::
MPSt(const MPSt<ITensor>&);
template MPSt<IQTensor>::
MPSt(const MPSt<IQTensor>&);

template <class Tensor>
MPSt<Tensor>& MPSt<Tensor>::
operator=(const MPSt& other)
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
template MPSt<ITensor>& MPSt<ITensor>::
operator=(const MPSt<ITensor>&);
template MPSt<IQTensor>& MPSt<IQTensor>::
operator=(const MPSt<IQTensor>&);

template <class Tensor>
MPSt<Tensor>::
~MPSt()
    {
    cleanupWrite();
    }
template MPSt<ITensor>::~MPSt();
template MPSt<IQTensor>::~MPSt();

template <class Tensor>
const Tensor& MPSt<Tensor>::
A(int i) const
    { 
    if(i < 0) i = N_+i+1;
    setSite(i);
    return A_.at(i); 
    }
template
const ITensor& MPSt<ITensor>::A(int i) const;
template
const IQTensor& MPSt<IQTensor>::A(int i) const;

template <class Tensor>
Tensor& MPSt<Tensor>::
Anc(int i)
    { 
    if(i < 0) i = N_+i+1;
    setSite(i);
    if(i <= l_orth_lim_) l_orth_lim_ = i-1;
    if(i >= r_orth_lim_) r_orth_lim_ = i+1;
    return A_.at(i); 
    }
template
ITensor& MPSt<ITensor>::Anc(int i);
template
IQTensor& MPSt<IQTensor>::Anc(int i);

template <class Tensor>
void MPSt<Tensor>::
doWrite(bool val, const Args& args) 
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
template void MPSt<ITensor>::
doWrite(bool val, const Args& args);
template void MPSt<IQTensor>::
doWrite(bool val, const Args& args);


template <class Tensor>
void MPSt<Tensor>::
read(std::istream& s)
    {
    if(not sites_) Error("Can't read to default constructed MPS");
    for(auto j : range(A_))
        {
        itensor::read(s,A_[j]);
        }
    //Check that tensors read from disk were constructed
    //using the same sites
    auto s1 = findtype(A_.at(1),Site);
    s1.noprime();
    if(s1 != IndexT(sites_->si(1)))
        {
        Print(A_.at(1).inds());
        Print(s1);
        Print(IndexT(sites_->si(1)));
        Error("Tensors read from disk not compatible with SiteSet passed to constructor.");
        }
    itensor::read(s,l_orth_lim_);
    itensor::read(s,r_orth_lim_);
    }
template
void MPSt<ITensor>::read(std::istream& s);
template
void MPSt<IQTensor>::read(std::istream& s);


template <class Tensor>
void MPSt<Tensor>::
write(std::ostream& s) const
    {
    if(do_write_)
        Error("MPSt::write not yet supported if doWrite(true)");

    for(auto j : range(A_.size()))
        {
        itensor::write(s,A_[j]);
        }
    itensor::write(s,leftLim());
    itensor::write(s,rightLim());
    }
template
void MPSt<ITensor>::write(std::ostream& s) const;
template
void MPSt<IQTensor>::write(std::ostream& s) const;

template <class Tensor>
void MPSt<Tensor>::
read(const std::string& dirname)
    {
    if(sites_ == 0)
        Error("Can't read to default constructed MPS, must specify SiteSet");

    l_orth_lim_ = 0;
    r_orth_lim_ = N_+1;

    //std::string dname_ = dirname;
    //if(dname_[dname_.length()-1] != '/')
    //    dname_ += "/";

    for(size_t j = 0; j < A_.size(); ++j) 
        {
    	readFromFile(AFName(j,dirname),A_.at(j));
        }
    }
template
void MPSt<ITensor>::read(const std::string& dirname);
template
void MPSt<IQTensor>::read(const std::string& dirname);


template <class Tensor>
string MPSt<Tensor>::
AFName(int j, const string& dirname) const
    { 
    if(dirname == "")
        return format("%s/A_%03d",writedir_,j);
    else
        return format("%s/A_%03d",dirname,j);
    }
template
string MPSt<ITensor>::AFName(int j, const string&) const;
template
string MPSt<IQTensor>::AFName(int j, const string&) const;

template <class Tensor>
void MPSt<Tensor>::
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
            A_.at(atb_) = Tensor();
            }
        if(A_.at(atb_+1))
            {
            writeToFile(AFName(atb_+1),A_.at(atb_+1));
            if(atb_+1 != b) A_.at(atb_+1) = Tensor();
            }
        ++atb_;
        }
    while(b < atb_)
        {
        if(A_.at(atb_))
            {
            writeToFile(AFName(atb_),A_.at(atb_));
            if(atb_ != b+1) A_.at(atb_) = Tensor();
            }
        if(A_.at(atb_+1))
            {
            writeToFile(AFName(atb_+1),A_.at(atb_+1));
            A_.at(atb_+1) = Tensor();
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
template
void MPSt<ITensor>::setBond(int b) const;
template
void MPSt<IQTensor>::setBond(int b) const;

template <class Tensor>
void MPSt<Tensor>::
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
template
void MPSt<ITensor>::setSite(int j) const;
template
void MPSt<IQTensor>::setSite(int j) const;


template <class Tensor>
void MPSt<Tensor>::
new_tensors(std::vector<ITensor>& A_)
    {
    std::vector<Index> a(N_+1);
    for(int i = 1; i <= N_; ++i)
        { a[i] = Index(nameint("a",i)); }
    A_[1] = ITensor(sites()(1),a[1]);
    for(int i = 2; i < N_; i++)
        { A_[i] = ITensor(dag(a[i-1]),sites()(i),a[i]); }
    A_[N_] = ITensor(dag(a[N_-1]),sites()(N_));
    }
template
void MPSt<ITensor>::new_tensors(std::vector<ITensor>& A_);
template
void MPSt<IQTensor>::new_tensors(std::vector<ITensor>& A_);

template <class Tensor>
void MPSt<Tensor>::
random_tensors(std::vector<ITensor>& A_)
    { 
    new_tensors(A_); 
    for(int i = 1; i <= N_; ++i)
        {
        randomize(A_[i]); 
        }
    }
template
void MPSt<ITensor>::random_tensors(std::vector<ITensor>& A_);
template
void MPSt<IQTensor>::random_tensors(std::vector<ITensor>& A_);

template <class Tensor>
void MPSt<Tensor>::
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
template
void MPSt<ITensor>::
init_tensors(std::vector<ITensor>& A_, const InitState& initState);


template <class Tensor>
void MPSt<Tensor>::
init_tensors(std::vector<IQTensor>& A_, const InitState& initState)
    {
    auto qa = std::vector<QN>(N_+1); //qn[i] = qn on i^th bond
    for(auto i : range1(N_)) qa[0] -= initState(i).qn()*In;

    //Taking OC to be at the leftmost site,
    //compute the QuantumNumbers of all the Links.
    for(auto i : range1(N_))
        {
        //Taking the divergence to be zero,solve for qa[i]
        qa[i] = Out*(-qa[i-1]*In - initState(i).qn());
        }

    auto a = std::vector<IQIndex>(N_+1);
    for(auto i : range1(N_))
        { 
        a[i] = IQIndex(nameint("L",i),Index(nameint("l",i)),qa[i]); 
        }

    A_[1] = setElt(initState(1),a[1](1));
    for(auto i : range(2,N_))
        {
        A_[i] = setElt(dag(a[i-1])(1),initState(i),a[i](1)); 
        }
    A_[N_] = setElt(dag(a[N_-1])(1),initState(N_));
    }
template
void MPSt<IQTensor>::
init_tensors(std::vector<IQTensor>& A_, const InitState& initState);


void 
plussers(Index const& l1, 
         Index const& l2, 
         Index      & sumind, 
         ITensor    & first, 
         ITensor    & second)
    {
    auto m = l1.m()+l2.m();
    if(m <= 0) m = 1;
    sumind = Index(sumind.rawname(),m);

    first = delta(l1,sumind);
    auto S = Matrix(l2.m(),sumind.m());
    for(auto i : range(l2.m()))
        {
        S(i,l1.m()+i) = 1;
        }
    second = matrixTensor(std::move(S),l2,sumind);
    }

void 
plussers(IQIndex const& l1, IQIndex const& l2, 
         IQIndex& sumind, 
         IQTensor& first, IQTensor& second)
    {
    map<Index,Index> l1map, l2map;
    vector<IndexQN> iq;
    for(IndexQN const& x : l1)
        {
        Index jj(x.index.rawname(),x.m(),x.type());
        l1map[x.index] = jj;
        iq.push_back(IndexQN(jj,x.qn));
        }
    for(IndexQN const& x : l2)
        {
        Index jj(x.index.rawname(),x.m(),x.type());
        l2map[x.index] = jj;
        iq.push_back(IndexQN(jj,x.qn));
        }
    sumind = IQIndex(sumind.rawname(),std::move(iq),sumind.dir(),sumind.primeLevel());
    first = IQTensor(dag(l1),sumind);
    for(IndexQN const& il1 : l1)
        {
        Index& s1 = l1map[il1.index];
        auto t = delta(il1.index,s1);
        first += t;
        }
    second = IQTensor(dag(l2),sumind);
    for(IndexQN const& il2 : l2)
        {
        Index& s2 = l2map[il2.index];
        auto t = delta(il2.index,s2);
        second += t;
        }
    }

//template <>
//MPSt<IQTensor>& MPSt<IQTensor>::operator+=(const MPSt<IQTensor>& other)
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
//MPSt<Tensor>& MPSt<Tensor>::operator+=(const MPSt<Tensor>& other)
//    {
//    if(do_write_)
//        Error("operator+= not supported if doWrite(true)");
//
//    primelinks(0,4);
//
//    vector<Tensor> first(N), second(N);
//    for(int i = 1; i < N_; ++i)
//        {
//        IndexT l1 = rightLinkInd(*this,i);
//        IndexT l2 = rightLinkInd(other,i);
//        IndexT r(l1);
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
//MPSt<ITensor>& MPSt<ITensor>::operator+=(const MPSt<ITensor>& other);

template <class Tensor>
MPSt<Tensor>& MPSt<Tensor>::
plusEq(const MPSt<Tensor>& R,
       const Args& args)
    {
    //cout << "calling new orthog in sum" << endl;
    if(!itensor::isOrtho(*this))
        {
        try { 
            orthogonalize(); 
            }
        catch(const ResultIsZero& rz) 
            { 
            *this = R;
            return *this;
            }
        }

    if(!itensor::isOrtho(R))
        {
        MPSt<Tensor> oR(R);
        try { 
            oR.orthogonalize(); 
            }
        catch(const ResultIsZero& rz) 
            { 
            return *this;
            }
        return addAssumeOrth(oR,args);
        }

    return addAssumeOrth(R,args);
    }
template
MPSt<ITensor>& MPSt<ITensor>::
plusEq(const MPSt<ITensor>& R, const Args& args);
template
MPSt<IQTensor>& MPSt<IQTensor>::
plusEq(const MPSt<IQTensor>& R, const Args& args);

//
// Adds two MPSs but doesn't attempt to
// orthogonalize them first
//
template <class Tensor>
MPSt<Tensor>& MPSt<Tensor>::
addAssumeOrth(MPSt<Tensor> const& R,
              Args const& args)
    {
    auto& L = *this;
    primelinks(0,4);

    vector<Tensor> first(N_), 
                   second(N_);
    for(auto i : range1(N_-1))
        {
        auto l1 = rightLinkInd(*this,i);
        auto l2 = rightLinkInd(R,i);
        auto r = l1;
        plussers(l1,l2,r,first[i],second[i]);
        }

    Anc(1) = L.A(1) * first[1] + R.A(1) * second[1];
    for(auto i : range1(2,N_-1))
        {
        Anc(i) = dag(first[i-1]) * L.A(i) * first[i] 
                     + dag(second[i-1]) * R.A(i) * second[i];
        }
    Anc(N_) = dag(first[N_-1]) * L.A(N_) + dag(second[N_-1]) * R.A(N_);

    noprimelink();

    orthogonalize(args);

    return *this;
    }
template
MPSt<ITensor>& MPSt<ITensor>::
addAssumeOrth(const MPSt<ITensor>& R, const Args& args);
template
MPSt<IQTensor>& MPSt<IQTensor>::
addAssumeOrth(const MPSt<IQTensor>& R, const Args& args);


//
//MPSt Index Methods
//

template <class Tensor>
void MPSt<Tensor>::
mapprime(int oldp, int newp, IndexType type)
    { 
    if(do_write_)
        Error("mapprime not supported if doWrite(true)");
    for(int i = 1; i <= N_; ++i) 
        A_[i].mapprime(oldp,newp,type); 
    }
template
void MPSt<ITensor>::mapprime(int oldp, int newp, IndexType type);
template
void MPSt<IQTensor>::mapprime(int oldp, int newp, IndexType type);

template <class Tensor>
void MPSt<Tensor>::
primelinks(int oldp, int newp)
    { 
    if(do_write_)
        Error("primelinks not supported if doWrite(true)");
    for(int i = 1; i <= N_; ++i) 
        A_[i].mapprime(oldp,newp,Link); 
    }
template
void MPSt<ITensor>::primelinks(int oldp, int newp);
template
void MPSt<IQTensor>::primelinks(int oldp, int newp);

template <class Tensor>
void MPSt<Tensor>::
noprimelink()
    { 
    if(do_write_)
        Error("noprimelink not supported if doWrite(true)");
    for(int i = 1; i <= N_; ++i) 
        A_[i].noprime(Link); 
    }
template
void MPSt<ITensor>::noprimelink();
template
void MPSt<IQTensor>::noprimelink();

template<class Tensor> 
Spectrum MPSt<Tensor>::
svdBond(int b, const Tensor& AA, Direction dir, const Args& args)
    {
    return svdBond(b,AA,dir,LocalOp<Tensor>(),args);
    }
template Spectrum MPSt<ITensor>::
svdBond(int b, const ITensor& AA, Direction dir, const Args& args);
template Spectrum MPSt<IQTensor>::
svdBond(int b, const IQTensor& AA, Direction dir, const Args& args);


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

template<class Tensor>
Spectrum
orthMPS(Tensor& A1, Tensor& A2, Direction dir, const Args& args)
    {
    using IndexT = typename Tensor::index_type;

    Tensor& L = (dir == Fromleft ? A1 : A2);
    Tensor& R = (dir == Fromleft ? A2 : A1);

    IndexT bnd = commonIndex(L,R,Link);
    if(!bnd) return Spectrum();

    if(args.getBool("Verbose",false))
        {
        Print(L.inds());
        }

    Tensor A,B(bnd);
    Tensor D;
    Spectrum spec = svd(L,A,D,B,args);

    L = A;
    R *= (D*B);

    //Older density matrix implementation
    //Doesn't flip arrows appropriately

    //Tensor rho = prime(L,bnd)*dag(L);

    //Tensor U;
    //Tensor D;
    //diagHermitian(rho,U,D,spec,args);


    //Tensor Di = D;
    //Di.mapElems(SqrtInv());
    //D.mapElems(Sqrt());

    //const
    //Tensor siRho = dag(U)*Di*prime(U),
    //       sRho = dag(U)*D*prime(U);

    //L *= siRho;
    //L.noprime();

    //R = prime(R,bnd)*sRho;

    return spec;
    }
template Spectrum
orthMPS(ITensor& A1, ITensor& A2, Direction dir, const Args& args);
template Spectrum
orthMPS(IQTensor& A1, IQTensor& A2, Direction dir, const Args& args);


template<class Tensor> 
void MPSt<Tensor>::
position(int i, const Args& args)
    {
    if(!(*this)) Error("position: MPS is default constructed");

    if(args.getBool("DoSVDBond",false))
        {
        while(l_orth_lim_ < i-1)
            {
            if(l_orth_lim_ < 0) l_orth_lim_ = 0;
            setBond(l_orth_lim_+1);
            Tensor WF = A(l_orth_lim_+1) * A(l_orth_lim_+2);
            svdBond(l_orth_lim_+1,WF,Fromleft,args);
            }
        while(r_orth_lim_ > i+1)
            {
            if(r_orth_lim_ > N_+1) r_orth_lim_ = N_+1;
            setBond(r_orth_lim_-2);
            Tensor WF = A(r_orth_lim_-2) * A(r_orth_lim_-1);
            svdBond(r_orth_lim_-2,WF,Fromright,args);
            }
        }
    else //use orthMPS
        {
        while(l_orth_lim_ < i-1)
            {
            if(l_orth_lim_ < 0) l_orth_lim_ = 0;
            setBond(l_orth_lim_+1);
            orthMPS(Anc(l_orth_lim_+1),Anc(l_orth_lim_+2),Fromleft,args);
            ++l_orth_lim_;
            if(r_orth_lim_ < l_orth_lim_+2) r_orth_lim_ = l_orth_lim_+2;
            }
        while(r_orth_lim_ > i+1)
            {
            if(r_orth_lim_ > N_+1) r_orth_lim_ = N_+1;
            setBond(r_orth_lim_-2);
            orthMPS(Anc(r_orth_lim_-2),Anc(r_orth_lim_-1),Fromright,args);
            --r_orth_lim_;
            if(l_orth_lim_ > r_orth_lim_-2) l_orth_lim_ = r_orth_lim_-2;
            }
        }
    }
template void MPSt<ITensor>::
position(int b, const Args& args);
template void MPSt<IQTensor>::
position(int b, const Args& args);

template <class Tensor>
int MPSt<Tensor>::
orthoCenter() const 
    { 
    if(!itensor::isOrtho(*this)) Error("orthogonality center not well defined.");
    return (leftLim() + 1);
    }
template
int MPSt<ITensor>::orthoCenter() const;
template
int MPSt<IQTensor>::orthoCenter() const;

template <class Tensor>
void MPSt<Tensor>::
orthogonalize(const Args& args)
    {
    //Do a half-sweep to the right, orthogonalizing each bond
    //but lower the cutoff since the basis to the right
    //might not be ortho: don't want to over truncate
    l_orth_lim_ = 0;
    r_orth_lim_ = N()+1;
    //Use smaller cutoff to orthogonalize w/ minimal truncation
    auto orig_cut = args.getReal("Cutoff",MIN_CUT);
    position(N_,{args,"Cutoff",0.1*orig_cut});
    //Now basis is ortho, ok to truncate
    position(1,args);
    }
template
void MPSt<ITensor>::orthogonalize(const Args& args);
template
void MPSt<IQTensor>::orthogonalize(const Args& args);

template <class Tensor>
void MPSt<Tensor>::
makeRealBasis(int j, const Args& args)
    {
    if(!(*this)) Error("position: MPS is default constructed");
    l_orth_lim_ = 0;
    while(l_orth_lim_ < j-1)
        {
        setBond(l_orth_lim_+1);
        Tensor WF = A(l_orth_lim_+1) * A(l_orth_lim_+2);
        orthoDecomp(WF,A_[l_orth_lim_+1],A_[l_orth_lim_+2],Fromleft,args);
        ++l_orth_lim_;
        }
    r_orth_lim_ = N_+1;
    while(r_orth_lim_ > j+1)
        {
        setBond(r_orth_lim_-2);
        Tensor WF = A(r_orth_lim_-2) * A(r_orth_lim_-1);
        orthoDecomp(WF,A_[r_orth_lim_-2],A_[r_orth_lim_-1],Fromright,args);
        --r_orth_lim_;
        }
    }
template
void MPSt<ITensor>::makeRealBasis(int j, const Args& args);
template
void MPSt<IQTensor>::makeRealBasis(int j, const Args& args);

//Methods for use internally by checkOrtho
ITensor
makeKroneckerDelta(const Index& i, int plev)
    {
    return delta(i,prime(i,plev));
    }
IQTensor
makeKroneckerDelta(const IQIndex& I, int plev)
    {
    IQTensor D(I,prime(I,plev));

    for(int j = 1; j <= I.nindex(); ++j)
        {
        D += makeKroneckerDelta(I.index(j),plev);
        }
    return D;
    }

//template <class Tensor>
//bool MPSt<Tensor>::
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
//bool MPSt<ITensor>::checkOrtho(int i, bool left) const;
//template
//bool MPSt<IQTensor>::checkOrtho(int i, bool left) const;

//template <class Tensor>
//bool MPSt<Tensor>::
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
//bool MPSt<ITensor>::checkOrtho() const;
//template
//bool MPSt<IQTensor>::checkOrtho() const;


//template <class Tensor>
//void MPSt<Tensor>::
//applygate(const Tensor& gate, const Args& args)
//    {
//    setBond(l_orth_lim_+1);
//    Tensor AA = A_.at(l_orth_lim_+1) * A_.at(l_orth_lim_+2) * gate;
//    AA.noprime();
//    svdBond(l_orth_lim_+1,AA,Fromleft,args);
//    }
//template
//void MPSt<ITensor>::applygate(const ITensor& gate,const Args& args);
//template
//void MPSt<IQTensor>::applygate(const IQTensor& gate,const Args& args);

//template <class Tensor>
//void MPSt<Tensor>::
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
//void MPSt<ITensor>::applygate(const BondGate<ITensor>& gate,const Args& args);
//template
//void MPSt<IQTensor>::applygate(const BondGate<IQTensor>& gate,const Args& args);

template <class Tensor>
Real MPSt<Tensor>::
norm() const 
    { 
    if(itensor::isOrtho(*this))
        {
        return itensor::norm(A(itensor::orthoCenter(*this)));
        }
    return std::sqrt(psiphi(*this,*this)); 
    }
template Real MPSt<ITensor>::
norm() const;
template Real MPSt<IQTensor>::
norm() const;

template <class Tensor>
Real MPSt<Tensor>::
normalize()
    {
    auto norm_ = norm();
    if(std::fabs(norm_) < 1E-20) Error("Zero norm");
    *this /= norm_;
    return norm_;
    }
template
Real MPSt<ITensor>::
normalize();
template
Real MPSt<IQTensor>::
normalize();

template <class Tensor>
bool MPSt<Tensor>::
isComplex() const
    { 
    for(auto j : range1(N_))
        {
        if(itensor::isComplex(A_[j])) return true;
        }
    return false;
    }
template
bool MPSt<ITensor>::isComplex() const;
template
bool MPSt<IQTensor>::isComplex() const;

template <class Tensor>
void MPSt<Tensor>::
initWrite(const Args& args)
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
                    A_[j] = Tensor();
                }
            }

        writeToFile(writedir_+"/sites",*sites_);

        do_write_ = true;
        }
    }
template
void MPSt<ITensor>::initWrite(const Args&);
template
void MPSt<IQTensor>::initWrite(const Args&);

template <class Tensor>
void MPSt<Tensor>::
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
template
void MPSt<ITensor>::copyWriteDir();
template
void MPSt<IQTensor>::copyWriteDir();


template <class Tensor>
void MPSt<Tensor>::
cleanupWrite()
    {
    if(do_write_)
        {
        const string cmdstr = "rm -fr " + writedir_;
        system(cmdstr.c_str());
        do_write_ = false;
        }   
    }
template
void MPSt<ITensor>::cleanupWrite();
template
void MPSt<IQTensor>::cleanupWrite();

template<class Tensor>
void MPSt<Tensor>::
swap(MPSt<Tensor>& other)
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
template
void MPSt<ITensor>::swap(MPSt<ITensor>& other);
template
void MPSt<IQTensor>::swap(MPSt<IQTensor>& other);

InitState::
InitState(const SiteSet& sites)
    : 
    sites_(&sites), 
    state_(1+sites.N())
    { 
    for(int n = 1; n <= sites_->N(); ++n)
        {
        state_[n] = sites_->si(n)(1);
        }
    }

InitState::
InitState(const SiteSet& sites, const String& state)
    : 
    sites_(&sites), 
    state_(1+sites.N())
    { 
    setAll(state);
    }

InitState& InitState::
set(int i, const String& state)
    { 
    checkRange(i);
    state_.at(i) = sites_->st(i,state);
    return *this;
    }

InitState& InitState::
setAll(const String& state)
    { 
    for(int n = 1; n <= sites_->N(); ++n)
        {
        state_[n] = sites_->st(n,state);
        }
    return *this;
    }

void InitState::
checkRange(int i) const
    {
    if(i > sites_->N() || i < 1) 
        {
        println("i = ",i);
        println("Valid range is 1 to ",sites_->N());
        Error("i out of range");
        }
    }

//Auxilary method for convertToIQ
long
collapseCols(const Vector& Diag, Matrix& M)
{
    long nr = Diag.size(),
            nc = long(sumels(Diag));
    assert(nr != 0);
    if(nc == 0) return nc;
    M = Matrix(nr,nc);
    long c = 0;
    for(long r = 0; r < nr; ++r)
        if(Diag(r) == 1) { M(r,c) = 1; ++c;}
    return nc;
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

void
convertToIQ(const SiteSet& sites, const vector<ITensor>& A,
             vector<IQTensor>& qA, QN totalq, Real cut)
{
    const int N = sites.N();
    qA.resize(A.size());
    const bool is_mpo = hasindex(A[1],sites.siP(1));
    const int Dim = sites.si(1).m();
    if(sites.si(2).m() != Dim)
        Error("convertToIQ assumes uniform site dimension");
    const int PDim = (is_mpo ? Dim : 1);

    // If MPO, set all tensors to identity ops initially
    if(is_mpo)
    {
        for(int j = 1; j <= N; ++j)
            qA.at(j) = sites.op("Id",j);
    }

    const int fullrank = (is_mpo ? 4 : 3);
    int start = 1, end = N;

    for(int j = 1; j <= N; ++j)
        if(A[j].r() == fullrank)
        if(A.at(periodicWrap(j-1,N)).r() < fullrank)
        {
            start = periodicWrap(j-1,N);
            //cout << "Got start at " << start << "\n";
            break;
        }

    for(int j = 1; j <= N; ++j)
        if(A[j].r() == fullrank)
        if(A.at(periodicWrap(j+1,N)).r() < fullrank)
        {
            end = periodicWrap(j+1,N);
            //cout << "Got end at " << end << "\n";
            break;
        }

//    cout << "Converting to IQ with (start, end) = " << start<< " "<< end << endl;

    vector<IQIndex> linkind(N+1);

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
    int Send = (end < start ? N+end : end);

    for(int S = start; S <= Send; ++S)
    {
        PrintData(S);
        int s = periodicWrap(S,N);
        int sprev = periodicWrap(S-1,N);
        int snext = periodicWrap(S+1,N);

        qD.clear();
        qt.clear();

        if(S > start) prev_bond = commonIndex(A[sprev],A[s],Link);
        if(S < Send) bond = commonIndex(A[s],A[snext],Link);

        if(s == show_s) { PrintData(A[s]); }

        for(const qC_vt& x : qC)
            for(int n = 1; n <= Dim;  ++n)
                for(int u = 1; u <= PDim; ++u)
                {
                    //Each compressor represents a particular
                    //QN channel given by prev_q
                    const QN& prev_q = x.first;
                    //Alias previous compressor ITensor to comp
                    const ITensor& comp = x.second;

                    q = (is_mpo ? prev_q+sites.si(s).qn(n)-sites.si(s).qn(u)
                                : prev_q-sites.si(s).qn(n));

                    //For the last site, only keep blocks
                    //compatible with specified totalq i.e. q=0 here
                    if(S == Send && q != QN()) continue;

                    //Set Site indices of A[s] and its previous Link Index
                    block = A[s];
                    if(S != start) block *= dag(comp);
                    block *= Index(sites.si(s))(n);
                    if(is_mpo) block *= Index(sites.siP(s))(u);

                    //Initialize D Vector (D records which values of
                    //the right Link Index to keep for the current QN q)
                    auto count = qD.count(q);
                    Vector& D = qD[q];
                    if(count == 0)
                    {
                        resize(D, bond.m());
                        for(auto& el : D) el = 0;
                    }

                    if(s == show_s)
                    {
                        println("For n = ",n);
                        printfln("Got a block with norm %.10f",norm(block));
                        println("bond.m() = ",bond.m());
                        PrintData(block);
                        if(s != 1) PrintData(comp);
                    }

                    bool keep_block = false;
                    if(S == Send)
                    { keep_block = true; }
                    else
                    {
                        if(bond.m() == 1 && norm(block) != 0)
                        {
                            for(auto& el : D) el = 1;
                            keep_block = true;
                        }
                        else
                        {
                            ITensor summed_block;
                            if(S==start)
                            { summed_block = block; }
                            else
                            {
                                //Here we sum over the previous link index
                                //which is already ok, analyze the one to the right
                                assert(comp.r()==2);
                                auto ci = comp.inds().begin();
                                const Index& new_ind = (*ci==prev_bond ? *(ci+1) : *ci);
                                summed_block = delta(new_ind) * block;
                            }
                            //summed_block.print("summed_block");

                            Real rel_cut = -1;
                            const ITensor& sb = summed_block;
                            for(int j = 1; j <= bond.m(); ++j)
                            { rel_cut = std::max(std::fabs(sb.real(bond(j))),rel_cut); }
                            assert(rel_cut >= 0);
                            //Real rel_cut = summed_block.norm()/summed_block.vecSize();
                            rel_cut *= cut;
                            //cout << "rel_cut == " << rel_cut << "\n";

                            if(rel_cut > 0)
                                for(int j = 1; j <= bond.m(); ++j)
                                {
                                    if(std::fabs(sb.real(bond(j))) > rel_cut)
                                    {
                                        D(j-1) = 1;
                                        keep_block = true;
                                    }
                                }
                        }
                    } //else (S != Send)

                    if(keep_block)
                    {
                        qD[q] = D;

                        auto newindbuilder = IndexSetBuilder(block.r());
                        for(auto& I: block.inds())
                            newindbuilder.nextIndex(std::move(I));

                        if(is_mpo)
                        {
                            newindbuilder.resize(block.r()+2);
                            newindbuilder.nextIndex(Index(dag(sites.si(s)(n).indexqn())));
                            newindbuilder.nextIndex(Index(sites.siP(s)(u).indexqn()));
                        }
                        else
                        {
                            newindbuilder.resize(block.r()+1);
                            newindbuilder.nextIndex(Index(sites.si(s)(n).indexqn()));
                        }

                        IndexSet newinds = newindbuilder.build();

                        if(s==show_s)
                        {
                            PrintData(block);
                            cout << "D = " << D << "\n";
                        }

                        qt[q].push_back(ITensor(std::move(newinds),std::move(block.store()),std::move(block.scale())));

                    }
                }

        qC.clear();

        for(const qt_vt& x : qt)
        {
            const vector<ITensor>& blks = x.second;

            if(blks.size() != 0)
            {
                q = x.first;
                if(S == Send)
                { for(const ITensor& t : blks) nblock.push_back(t); }
                else
                {
                    Matrix M;
                    auto mm = collapseCols(qD[q],M);
                    if(s==show_s)
                    {
                        println("Adding block, mm = ",mm);
                        Print(q);
                        cout << "qD[q] = " << qD[q] << "\n";
                        cout << "M = \n" << M << "\n";
                        int count = 0;
                        for(const ITensor& t : blks)
                        {
                            printfln("t%02d",++count," ",t);
                        }
                    }
                    string qname = format("ql%d(%+d:%d)",s,Sz(q),Nf(q));
                    Index qbond(qname,mm);
                    auto compressor = matrixTensor(std::move(M),bond,qbond);
                    for(const ITensor& t : blks) nblock.push_back(t * compressor);
                    iq.push_back(IndexQN(qbond,q));
                    qC[q] = compressor;
                }
            }
        }

        if(S != Send)
        {
            if(iq.empty())
            {
                cout << "At site " << s << "\n";
                Error("convertToIQ: no compatible QNs to put into Link.");
            }
            linkind[s] = IQIndex(nameint("qL",s),std::move(iq));
        }
        if(S == start)
        {
            qA.at(s) = (is_mpo ? IQTensor(dag(sites.si(s)),sites.siP(s),linkind.at(s))
                               : IQTensor(sites.si(s),linkind[s]));
        }
        else
        if(S == Send)
        {
            qA.at(s) = (is_mpo ? IQTensor(dag(linkind[sprev]),dag(sites.si(s)),sites.siP(s))
                               : IQTensor(dag(linkind[sprev]),sites.si(s)));
        }
        else
        {
            qA.at(s) = (is_mpo ? IQTensor(dag(linkind[sprev]),dag(sites.si(s)),sites.siP(s),linkind[s])
                               : IQTensor(dag(linkind[sprev]),sites.si(s),linkind[s]));
        }

        for(const ITensor& nb : nblock)
        { qA.at(s) += nb; }
        nblock.clear();

        if(s==show_s)
        {
            printfln("qA[%d]",s,qA[s]);
            Error("Stopping");
        }

    } //for loop over s

} //void convertToIQ

/*
template <class Tensor> 
template <class IQMPSType> 
void MPSt<Tensor>::convertToIQ(IQMPSType& iqpsi, QN totalq, Real cut) const
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
findCenter(const IQMPS& psi)
    {
    for(int j = 1; j <= psi.N(); ++j) 
        {
        const IQTensor& A = psi.A(j);
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

bool 
checkQNs(const IQMPS& psi)
    {
    const int N = psi.N();

    QN Zero;

    int center = findCenter(psi);
    if(center == -1)
        {
        cout << "Did not find an ortho. center\n";
        return false;
        }

    //Check that all IQTensors have zero div
    //except possibly the ortho. center
    for(int i = 1; i <= N; ++i) 
        {
        if(i == center) continue;
        if(!psi.A(i))
            {
            println("A(",i,") null, QNs not well defined");
            return false;
            }
        if(div(psi.A(i)) != Zero)
            {
            cout << "At i = " << i << "\n";
            Print(psi.A(i));
            cout << "IQTensor other than the ortho center had non-zero divergence\n";
            return false;
            }
        }

    //Check arrows from left edge
    for(int i = 1; i < center; ++i)
        {
        if(rightLinkInd(psi,i).dir() != In) 
            {
            println("checkQNs: At site ",i," to the left of the OC, Right side Link not pointing In");
            return false;
            }
        if(i > 1)
            {
            if(leftLinkInd(psi,i).dir() != Out) 
                {
                println("checkQNs: At site ",i," to the left of the OC, Left side Link not pointing Out");
                return false;
                }
            }
        }

    //Check arrows from right edge
    for(int i = N; i > center; --i)
        {
        if(i < N)
        if(rightLinkInd(psi,i).dir() != Out) 
            {
            println("checkQNs: At site ",i," to the right of the OC, Right side Link not pointing Out");
            return false;
            }
        if(leftLinkInd(psi,i).dir() != In) 
            {
            println("checkQNs: At site ",i," to the right of the OC, Left side Link not pointing In");
            return false;
            }
        }

    //Done checking arrows
    return true;
    }

QN
totalQN(const IQMPS& psi)
    {
    const int center = findCenter(psi);
    if(center == -1)
        Error("Could not find ortho. center");
    return div(psi.A(center));
    }

template <class Tensor>
void 
fitWF(const MPSt<Tensor>& psi_basis, MPSt<Tensor>& psi_to_fit)
    {
    if(!itensor::isOrtho(psi_basis)) 
        Error("psi_basis must be orthogonolized.");
    if(orthoCenter(psi_basis) != 1) 
        Error("psi_basis must be orthogonolized to site 1.");

    auto N = psi_basis.N();
    if(psi_to_fit.N() != N) 
        Error("Wavefunctions must have same number of sites.");

    auto A = psi_to_fit.A(N) * dag(prime(psi_basis.A(N),Link));
    for(int n = N-1; n > 1; --n)
        {
        A *= dag(prime(psi_basis.A(n),Link));
        A *= psi_to_fit.A(n);
        }
    A = psi_to_fit.A(1) * A;
    A.noprime();

    auto nrm = norm(A);
    if(nrm == 0) Error("Zero overlap of psi_to_fit and psi_basis");
    A /= nrm;

    psi_to_fit = psi_basis;
    psi_to_fit.Anc(1) = A;
    }
template void fitWF(const MPSt<ITensor>& psi_basis, MPSt<ITensor>& psi_to_fit);
template void fitWF(const MPSt<IQTensor>& psi_basis, MPSt<IQTensor>& psi_to_fit);

template <class Tensor>
std::ostream& 
operator<<(std::ostream& s, const MPSt<Tensor>& M)
    {
    s << "\n";
    for(int i = 1; i <= M.N(); ++i) 
        s << M.A(i) << "\n";
    return s;
    }
template std::ostream& operator<<(std::ostream& s, const MPSt<ITensor>& M);
template std::ostream& operator<<(std::ostream& s, const MPSt<IQTensor>& M);

std::ostream& 
operator<<(std::ostream& s, const InitState& state)
    {
    s << "\n";
    for(int i = 1; i <= state.sites().N(); ++i) 
        s << state(i) << "\n";
    return s;
    }

} //namespace itensor
