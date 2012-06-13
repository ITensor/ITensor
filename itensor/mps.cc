//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#include "mps.h"
#include "localmpo.h"

using namespace std;
using boost::format;

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
    N(0), 
    model_(0),
    atb_(1),
    writedir_("psi"),
    do_write_(false)
    { }
template MPSt<ITensor>::
MPSt();
template MPSt<IQTensor>::
MPSt();

template <class Tensor>
MPSt<Tensor>::
MPSt(const Model& mod_,int maxmm, Real cut) 
    : 
    N(mod_.NN()), 
    A(mod_.NN()+1),
    l_orth_lim_(0),
    r_orth_lim_(mod_.NN()),
    model_(&mod_), 
    svd_(N,cut,1,maxmm,false,LogNumber(1)),
    atb_(1),
    writedir_("psi"),
    do_write_(false)
    { 
    random_tensors(A);
    }
template MPSt<ITensor>::
MPSt(const Model& mod_,int maxmm, Real cut);
template MPSt<IQTensor>::
MPSt(const Model& mod_,int maxmm, Real cut);

template <class Tensor>
MPSt<Tensor>::
MPSt(const Model& mod_,const InitState& initState,int maxmm, Real cut)
    : 
    N(mod_.NN()),
    A(mod_.NN()+1),
    l_orth_lim_(0),
    r_orth_lim_(2),
    model_(&mod_), 
    svd_(N,cut,1,maxmm,false,LogNumber(1)),
    atb_(1),
    writedir_("psi"),
    do_write_(false)
    { 
    init_tensors(A,initState);
    }
template MPSt<ITensor>::
MPSt(const Model& mod_,const InitState& initState,int maxmm, Real cut);
template MPSt<IQTensor>::
MPSt(const Model& mod_,const InitState& initState,int maxmm, Real cut);

template <class Tensor>
MPSt<Tensor>::
MPSt(const Model& model, std::istream& s)
    : 
    N(model.NN()), 
    A(model.NN()+1), 
    model_(&model),
    atb_(1),
    writedir_("psi"),
    do_write_(false)
    { 
    read(s); 
    }
template MPSt<ITensor>::
MPSt(const Model& model, std::istream& s);
template MPSt<IQTensor>::
MPSt(const Model& model, std::istream& s);

template <class Tensor>
Tensor& MPSt<Tensor>::
AAnc(int i) //nc means 'non const'
    { 
    setSite(i);
    if(i <= l_orth_lim_) l_orth_lim_ = i-1;
    if(i >= r_orth_lim_) r_orth_lim_ = i+1;
    return A.at(i); 
    }
template
ITensor& MPSt<ITensor>::AAnc(int i);
template
IQTensor& MPSt<IQTensor>::AAnc(int i);

template <class Tensor>
Tensor MPSt<Tensor>::
bondTensor(int b) const 
    { 
    setBond(b);
    Tensor res = A.at(b) * A.at(b+1); 
    return res; 
    }
template
ITensor MPSt<ITensor>::bondTensor(int b) const;
template
IQTensor MPSt<IQTensor>::bondTensor(int b) const;


template <class Tensor>
void MPSt<Tensor>::
read(std::istream& s)
    {
    if(model_ == 0)
        Error("Can't read to default constructed MPS");
    for(int j = 1; j <= N; ++j) 
        A.at(j).read(s);
    s.read((char*) &l_orth_lim_,sizeof(l_orth_lim_));
    s.read((char*) &r_orth_lim_,sizeof(r_orth_lim_));
    svd_.read(s);
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

    for(int j = 1; j <= N; ++j) 
        {
        A.at(j).write(s);
        }
    s.write((char*) &l_orth_lim_,sizeof(l_orth_lim_));
    s.write((char*) &r_orth_lim_,sizeof(r_orth_lim_));
    svd_.write(s);
    }
template
void MPSt<ITensor>::write(std::ostream& s) const;
template
void MPSt<IQTensor>::write(std::ostream& s) const;

template <class Tensor>
string MPSt<Tensor>::
AFName(int j) const
    { 
    return (format("%s/A_%03d")%writedir_%j).str();
    }
template
string MPSt<ITensor>::AFName(int j) const;
template
string MPSt<IQTensor>::AFName(int j) const;

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
    //
    //Shift atb_ (location of bond that is loaded into RAM)
    //to requested value b, writing any non-Null tensors to
    //disk along the way
    //
    while(b > atb_)
        {
        if(A.at(atb_).isNotNull())
            {
            //std::cerr << boost::format("Writing A(%d) to %s\n")%atb_%writedir_;
            writeToFile(AFName(atb_),A.at(atb_));
            A.at(atb_) = Tensor();
            }
        if(A.at(atb_+1).isNotNull())
            {
            //std::cerr << boost::format("Writing A(%d) to %s\n")%(atb_+1)%writedir_;
            writeToFile(AFName(atb_+1),A.at(atb_+1));
            if(atb_+1 != b) A.at(atb_+1) = Tensor();
            }
        ++atb_;
        }
    while(b < atb_)
        {
        if(A.at(atb_).isNotNull())
            {
            //std::cerr << boost::format("Writing A(%d) to %s\n")%atb_%writedir_;
            writeToFile(AFName(atb_),A.at(atb_));
            if(atb_ != b+1) A.at(atb_) = Tensor();
            }
        if(A.at(atb_+1).isNotNull())
            {
            //std::cerr << boost::format("Writing A(%d) to %s\n")%(atb_+1)%writedir_;
            writeToFile(AFName(atb_+1),A.at(atb_+1));
            A.at(atb_+1) = Tensor();
            }
        --atb_;
        }
    assert(atb_ == b);
    //
    //Load tensors at bond b into RAM if
    //they aren't loaded already
    //
    if(A.at(b).isNull())
        {
        std::string fname = AFName(b);
        std::ifstream s(fname.c_str());
        if(s.good())
            {
            A.at(b).read(s);
            s.close();
            }
        else
            {
            std::cerr << boost::format("Tried to read file %s\n")%fname;
            Error("Missing file");
            }
        }
    if(A.at(b+1).isNull())
        {
        std::string fname = AFName(b+1);
        std::ifstream s(fname.c_str());
        if(s.good())
            {
            A.at(b+1).read(s);
            s.close();
            }
        else
            {
            std::cerr << boost::format("Tried read A[%d]\n")%(b+1);
            Error("Missing file");
            }
        }
    if(b == 1)
        {
        writeToFile(writedir_+"/model",*model_);
        std::ofstream inf((boost::format("%s/info")%writedir_).str().c_str());
            inf.write((char*) &l_orth_lim_,sizeof(l_orth_lim_));
            inf.write((char*) &r_orth_lim_,sizeof(r_orth_lim_));
            svd_.write(inf);
        inf.close();
        }
    }
template
void MPSt<ITensor>::setBond(int b) const;
template
void MPSt<IQTensor>::setBond(int b) const;


template <class Tensor>
void MPSt<Tensor>::
new_tensors(std::vector<ITensor>& A_)
    {
    std::vector<Index> a(N+1);
    for(int i = 1; i <= N; ++i)
        { a[i] = Index(nameint("a",i)); }
    A_[1] = ITensor(si(1),a[1]);
    for(int i = 2; i < N; i++)
        { A_[i] = ITensor(conj(a[i-1]),si(i),a[i]); }
    A_[N] = ITensor(conj(a[N-1]),si(N));
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
    for(int i = 1; i <= N; ++i)
        A_[i].Randomize(); 
    }
template
void MPSt<ITensor>::random_tensors(std::vector<ITensor>& A_);
template
void MPSt<IQTensor>::random_tensors(std::vector<ITensor>& A_);

template <class Tensor>
void MPSt<Tensor>::
init_tensors(std::vector<ITensor>& A_, const InitState& initState)
    { 
    new_tensors(A_); 
    for(int i = 1; i <= N; ++i) 
        {
        A_[i] = ITensor(initState(i)); 
        }

    std::vector<Index> a(N+1);
    for(int i = 1; i <= N; ++i)
        { a[i] = Index(nameint("l",i)); }

    A_[1].addindex1(a[1]);
    for(int i = 2; i < N; ++i)
        {
        A_[i].addindex1(a[i-1]);
        A_[i].addindex1(a[i]);
        }
    A_[N].addindex1(a[N-1]);
    }
template
void MPSt<ITensor>::
init_tensors(std::vector<ITensor>& A_, const InitState& initState);


template <class Tensor>
void MPSt<Tensor>::
init_tensors(std::vector<IQTensor>& A_, const InitState& initState)
    {
    std::vector<QN> qa(N+1); //qn[i] = qn on i^th bond
    for(int i = 1; i <= N; ++i) { qa[0] -= initState(i).qn()*In; }

    //Taking OC to be at the leftmost site,
    //compute the QuantumNumbers of all the Links.
    for(int i = 1; i <= N; ++i)
        {
        //Taking the divergence to be zero,solve for qa[i]
        qa[i] = Out*(-qa[i-1]*In - initState(i).qn());
        }

    std::vector<IQIndex> a(N+1);
    for(int i = 1; i <= N; ++i)
        { a[i] = IQIndex(nameint("L",i),Index(nameint("l",i)),qa[i]); }

    A_[1] = IQTensor(initState(1),a[1](1));
    for(int i = 2; i < N; ++i)
        A_[i] = IQTensor(conj(a[i-1])(1),initState(i),a[i](1)); 
    A_[N] = IQTensor(conj(a[N-1])(1),initState(N));
    }
template
void MPSt<IQTensor>::
init_tensors(std::vector<IQTensor>& A_, const InitState& initState);


template <class Tensor>
int MPSt<Tensor>::
averageM() const
    {
    Real avgm = 0;
    for(int b = 1; b < NN(); ++b)
        avgm += LinkInd(b).m();
    avgm /= (NN()-1);
    return (int)avgm;
    }
template
int MPSt<ITensor>::averageM() const;
template
int MPSt<IQTensor>::averageM() const;


void 
plussers(const Index& l1, const Index& l2, Index& sumind, 
          ITensor& first, ITensor& second)
    {
    sumind = Index(sumind.rawname(),l1.m()+l2.m(),sumind.type());
    first = ITensor(l1,sumind,1);
    second = ITensor(l2,sumind);
    for(int i = 1; i <= l2.m(); ++i) second(l2(i),sumind(l1.m()+i)) = 1;
    }

void 
plussers(const IQIndex& l1, const IQIndex& l2, IQIndex& sumind, 
         IQTensor& first, IQTensor& second)
    {
    map<Index,Index> l1map, l2map;
    vector<inqn> iq;
    Foreach(const inqn& x, l1.iq())
        {
        Index ii = x.index;
        Index jj(ii.name(),ii.m(),ii.type());
        l1map[ii] = jj;
        iq.push_back(inqn(jj,x.qn));
        }
    Foreach(const inqn& x, l2.iq())
        {
        Index ii = x.index;
        Index jj(ii.name(),ii.m(),ii.type());
        l2map[ii] = jj;
        iq.push_back(inqn(jj,x.qn));
        }
    sumind = IQIndex(sumind,iq);
    first = IQTensor(conj(l1),sumind);
    Foreach(const inqn& x, l1.iq())
        {
        Index il1 = x.index;
        Index s1 = l1map[il1];
        ITensor t(il1,s1,1.0);
        first += t;
        }
    second = IQTensor(conj(l2),sumind);
    Foreach(const inqn& x, l2.iq())
        {
        Index il2 = x.index;
        Index s2 = l2map[il2];
        ITensor t(il2,s2,1.0);
        second += t;
        }
    }

//#define NEW_MPS_ADDITION

#ifdef NEW_MPS_ADDITION

template <>
MPSt<IQTensor>& MPSt<IQTensor>::operator+=(const MPSt<IQTensor>& other)
    {
    if(do_write_)
        Error("operator+= not supported if doWrite(true)");

    primelinks(0,4);

    //Create new link indices
    vector<IQIndex> nlinks(N);
    for(int b = 1; b < N; ++b)
        {
        IQIndex l1 = this->LinkInd(b);
        IQIndex l2 = other.LinkInd(b);
        vector<inqn> iq(l1.iq());
        iq.insert(iq.begin(),l2.iq().begin(),l2.iq().end());
        nlinks.at(b) = IQIndex(l2,iq);
        }
    //Create new A tensors
    vector<IQTensor> nA(N+1);
    nA[1] = IQTensor(si(1),nlinks[1]);
    for(int j = 2; j < N; ++j)
        nA[j] = IQTensor(conj(nlinks[j-1]),si(j),nlinks[j]);
    nA[N] = IQTensor(conj(nlinks[N-1]),si(N));

    for(int j = 1; j <= N; ++j)
        {
        Foreach(const ITensor& t, AA(j).itensors())
            { nA[j].insert(t); }
        Foreach(const ITensor& t, other.AA(j).itensors())
            { nA[j].insert(t); }
        }

    A.swap(nA);

    orthogonalize();

    return *this;
    }

template <class Tensor>
MPSt<Tensor>& MPSt<Tensor>::operator+=(const MPSt<Tensor>& other)
    {
    if(do_write_)
        Error("operator+= not supported if doWrite(true)");

    primelinks(0,4);

    vector<Tensor> first(N), second(N);
    for(int i = 1; i < N; ++i)
        {
        IndexT l1 = this->RightLinkInd(i);
        IndexT l2 = other.RightLinkInd(i);
        IndexT r(l1.rawname());
        plussers(l1,l2,r,first[i],second[i]);
        }

    AAnc(1) = AA(1) * first[1] + other.AA(1) * second[1];
    for(int i = 2; i < N; ++i)
        {
        AAnc(i) = conj(first[i-1]) * AA(i) * first[i] 
                  + conj(second[i-1]) * other.AA(i) * second[i];
        }
    AAnc(N) = conj(first[N-1]) * AA(N) + conj(second[N-1]) * other.AA(N);

    noprimelink();

    //cerr << "WARNING: skipping orthogonalize in operator+=\n";
    orthogonalize();

    return *this;
    }
template
MPSt<ITensor>& MPSt<ITensor>::operator+=(const MPSt<ITensor>& other);

#else
template <class Tensor>
MPSt<Tensor>& MPSt<Tensor>::operator+=(const MPSt<Tensor>& other)
    {
    if(do_write_)
        Error("operator+= not supported if doWrite(true)");

    primelinks(0,4);

    vector<Tensor> first(N), second(N);
    for(int i = 1; i < N; ++i)
        {
        IndexT l1 = this->RightLinkInd(i);
        IndexT l2 = other.RightLinkInd(i);
        IndexT r(l1.rawname());
        plussers(l1,l2,r,first[i],second[i]);
        }

    AAnc(1) = AA(1) * first[1] + other.AA(1) * second[1];
    for(int i = 2; i < N; ++i)
        {
        AAnc(i) = conj(first[i-1]) * AA(i) * first[i] 
                  + conj(second[i-1]) * other.AA(i) * second[i];
        }
    AAnc(N) = conj(first[N-1]) * AA(N) + conj(second[N-1]) * other.AA(N);

    noprimelink();

    //cerr << "WARNING: skipping orthogonalize in operator+=\n";
    orthogonalize();

    return *this;
    }
template
MPSt<ITensor>& MPSt<ITensor>::operator+=(const MPSt<ITensor>& other);
template
MPSt<IQTensor>& MPSt<IQTensor>::operator+=(const MPSt<IQTensor>& other);
#endif


//
//MPSt Index Methods
//

template <class Tensor>
void MPSt<Tensor>::
mapprime(int oldp, int newp, PrimeType pt)
    { 
    if(do_write_)
        Error("mapprime not supported if doWrite(true)");
    for(int i = 1; i <= N; ++i) 
        A[i].mapprime(oldp,newp,pt); 
    }
template
void MPSt<ITensor>::mapprime(int oldp, int newp, PrimeType pt);
template
void MPSt<IQTensor>::mapprime(int oldp, int newp, PrimeType pt);

template <class Tensor>
void MPSt<Tensor>::
primelinks(int oldp, int newp)
    { 
    if(do_write_)
        Error("primelinks not supported if doWrite(true)");
    for(int i = 1; i <= N; ++i) 
        A[i].mapprime(oldp,newp,primeLink); 
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
    for(int i = 1; i <= N; ++i) 
        A[i].noprime(primeLink); 
    }
template
void MPSt<ITensor>::noprimelink();
template
void MPSt<IQTensor>::noprimelink();

template<class Tensor> void
MPSt<Tensor>::
svdBond(int b, const Tensor& AA, Direction dir, bool preserve_shape)
    {
    svdBond(b,AA,dir,LocalMPO<Tensor>::Null(),preserve_shape);
    }
template void MPSt<ITensor>::
svdBond(int b, const ITensor& AA, Direction dir, bool);
template void MPSt<IQTensor>::
svdBond(int b, const IQTensor& AA, Direction dir, bool);


template<class Tensor> void
MPSt<Tensor>::
position(int i, bool preserve_shape)
    {
    if(isNull()) Error("position: MPS is null");
    while(l_orth_lim_ < i-1)
        {
        if(l_orth_lim_ < 0) l_orth_lim_ = 0;
        setBond(l_orth_lim_+1);
        Tensor WF = AA(l_orth_lim_+1) * AA(l_orth_lim_+2);
        svdBond(l_orth_lim_+1,WF,Fromleft,preserve_shape);
        }
    while(r_orth_lim_ > i+1)
        {
        if(r_orth_lim_ > N+1) r_orth_lim_ = N+1;
        setBond(r_orth_lim_-2);
        Tensor WF = AA(r_orth_lim_-2) * AA(r_orth_lim_-1);
        svdBond(r_orth_lim_-2,WF,Fromright,preserve_shape);
        }
    }
template void MPSt<ITensor>::
position(int b, bool preserve_shape);
template void MPSt<IQTensor>::
position(int b, bool preserve_shape);

template <class Tensor>
void MPSt<Tensor>::
orthogonalize(bool verbose)
    {
    //Do a half-sweep to the right, orthogonalizing each bond
    //but do not truncate since the basis to the right might not
    //be ortho (i.e. use the current m).
    svd_.useOrigM(true);
    position(N);
    if(verbose)
        {
        std::cout << "Done orthogonalizing, starting truncation." 
                  << std::endl;
        }
    //Now basis is ortho, ok to truncate
    svd_.useOrigM(false);
    position(1);
    }
template
void MPSt<ITensor>::orthogonalize(bool verbose);
template
void MPSt<IQTensor>::orthogonalize(bool verbose);

//Methods for use internally by checkOrtho
ITensor
makeKroneckerDelta(const Index& i, int plev)
    {
    return ITensor(i,primed(i,plev),1);
    }
IQTensor
makeKroneckerDelta(const IQIndex& I, int plev)
    {
    IQTensor D(I,primed(I,plev));

    for(int j = 1; j <= I.nindex(); ++j)
        {
        D += makeKroneckerDelta(I.index(j),plev);
        }
    return D;
    }

template <class Tensor>
bool MPSt<Tensor>::
checkOrtho(int i, bool left) const
    {
    setSite(i);
    IndexT link = (left ? RightLinkInd(i) : LeftLinkInd(i));
    Tensor A = AA(i);
    Tensor Ac = conj(A); 
    Ac.primeind(link,4);

    Tensor rho = A * Ac;

    Tensor Delta = makeKroneckerDelta(link,4);

    Tensor Diff = rho - Delta;

    Vector diff(Diff.vecSize());
    Diff.assignToVec(diff);

    Real threshold = 1E-13;
    if(Norm(diff) < threshold) return true;

    //Print any helpful debugging info here:
    std::cerr << "checkOrtho: on line " << __LINE__ 
              << " of mps.h," << std::endl;
    std::cerr << "checkOrtho: Tensor at position " << i 
              << " failed to be " << (left ? "left" : "right") 
              << " ortho." << std::endl;
    std::cerr << "checkOrtho: Norm(diff) = " << boost::format("%E") 
              % Norm(diff) << std::endl;
    std::cerr << "checkOrtho: Error threshold set to " 
              << boost::format("%E") % threshold << std::endl;
    //-----------------------------

    return false;
    }
template
bool MPSt<ITensor>::checkOrtho(int i, bool left) const;
template
bool MPSt<IQTensor>::checkOrtho(int i, bool left) const;

template <class Tensor>
bool MPSt<Tensor>::
checkOrtho() const
    {
    for(int i = 1; i <= l_orth_lim_; ++i)
    if(!checkLeftOrtho(i))
        {
        std::cerr << "checkOrtho: A[i] not left orthogonal at site i=" 
                  << i << std::endl;
        return false;
        }

    for(int i = NN(); i >= r_orth_lim_; --i)
    if(!checkRightOrtho(i))
        {
        std::cerr << "checkOrtho: A[i] not right orthogonal at site i=" 
                  << i << std::endl;
        return false;
        }
    return true;
    }
template
bool MPSt<ITensor>::checkOrtho() const;
template
bool MPSt<IQTensor>::checkOrtho() const;

template <class Tensor>
void MPSt<Tensor>::
projectOp(int j, Direction dir, 
          const Tensor& E, const Tensor& X, Tensor& nE) const
    {
    if(dir == Fromleft)
        setBond(j);
    else
        setBond(j-1);

    if(dir==Fromleft && j > l_orth_lim_) 
        { 
        std::cerr << boost::format("projectOp: from left j > l_orth_lim_ (j=%d,l_orth_lim_=%d)\n")%j%l_orth_lim_; 
        Error("Projecting operator at j > l_orth_lim_"); 
        }
    if(dir==Fromright && j < r_orth_lim_) 
        { 
        std::cerr << boost::format("projectOp: from left j < r_orth_lim_ (j=%d,r_orth_lim_=%d)\n")%j%r_orth_lim_; 
        Error("Projecting operator at j < r_orth_lim_"); 
        }
    nE = (E.isNull() ? AA(j) : E * AA(j));
    nE *= X; 
    nE *= conj(primed(AA(j)));
    }
template
void MPSt<ITensor>::projectOp(int j, Direction dir, 
                              const ITensor& E, const ITensor& X, 
                              ITensor& nE) const;
template
void MPSt<IQTensor>::projectOp(int j, Direction dir, 
                              const IQTensor& E, const IQTensor& X, 
                              IQTensor& nE) const;

template <class Tensor>
void MPSt<Tensor>::
applygate(const Tensor& gate)
    {
    setBond(l_orth_lim_+1);
    Tensor AA = A[l_orth_lim_+1] * A[l_orth_lim_+2] * gate;
    AA.noprime();
    doSVD(l_orth_lim_+1,AA,Fromleft);
    }
template
void MPSt<ITensor>::applygate(const ITensor& gate);
template
void MPSt<IQTensor>::applygate(const IQTensor& gate);


//Auxilary method for convertToIQ
int 
collapseCols(const Vector& Diag, Matrix& M)
    {
    int nr = Diag.Length(), nc = int(Diag.sumels());
    assert(nr != 0);
    if(nc == 0) return nc;
    M = Matrix(nr,nc); M = 0;
    int c = 0;
    for(int r = 1; r <= nr; ++r)
    if(Diag(r) == 1) { M(r,++c) = 1; }
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
convertToIQ(const Model& model, const vector<ITensor>& A, 
            vector<IQTensor>& qA, QN totalq, Real cut)
    {
    const int N = A.size()-1;
    qA.resize(A.size());
    const bool is_mpo = A[1].hasindex(model.siP(1));
    const int Dim = model.si(1).m();
    if(model.si(2).m() != Dim)
        Error("convertToIQ assumes uniform site dimension");
    const int PDim = (is_mpo ? Dim : 1);

    // If MPO, set all tensors to identity ops initially
    if(is_mpo)
        {
        for(int j = 1; j <= N; ++j)
            qA.at(j) = model.id(j);
        }

    const int fullrank = (is_mpo ? 4 : 3);
    int start = 1, end = N;

    for(int j = 1; j <= N; ++j)
        if(A[j].r() == fullrank)
            if(A.at(periodicWrap(j-1,N)).r() < fullrank) 
                {
                start = periodicWrap(j-1,N);
                //cerr << "Got start at " << start << "\n";
                break;
                }

    for(int j = 1; j <= N; ++j)
        if(A[j].r() == fullrank)
            if(A.at(periodicWrap(j+1,N)).r() < fullrank) 
                {
                end = periodicWrap(j+1,N);
                //cerr << "Got end at " << end << "\n";
                break;
                }

    //cout << "Converting to IQ with (start, end) = " << start SP end << endl;

    vector<IQIndex> linkind(N+1);

    typedef map<QN,Vector>::value_type qD_vt;
    map<QN,Vector> qD; //Diags of compressor matrices by QN

    typedef map<QN,vector<ITensor> >::value_type qt_vt;
    map<QN,vector<ITensor> > qt; //ITensor blocks by QN

    typedef map<QN,ITensor>::value_type qC_vt;
    map<QN,ITensor> qC; //Compressor ITensors by QN

    ITensor block;
    vector<ITensor> nblock;
    vector<inqn> iq;

    QN q;

    qC[totalq] = ITensor(); //Represents Virtual index
    //First value of prev_q below set to totalq

    const int show_s = 0;

    Index bond, prev_bond;
    int Send = (end < start ? N+end : end); 
    for(int S = start; S <= Send; ++S)
        {
        int s = periodicWrap(S,N);
        int sprev = periodicWrap(S-1,N);
        int snext = periodicWrap(S+1,N);
        //cerr << format("S = %d, s = %d, sprev = %d, snext = %d\n")%S%s%sprev%snext;

        qD.clear(); 
        qt.clear();

        if(S > start) prev_bond = index_in_common(A[sprev],A[s],Link);
        if(S < Send) bond = index_in_common(A[s],A[snext],Link);

        if(s == show_s) { PrintDat(A[s]); }

        Foreach(const qC_vt& x, qC) 
        for(int n = 1; n <= Dim;  ++n)
        for(int u = 1; u <= PDim; ++u)
            {
            //Each compressor represents a particular
            //QN channel given by prev_q
            const QN& prev_q = x.first; 
            //Alias previous compressor ITensor to comp
            const ITensor& comp = x.second; 

            q = (is_mpo ? prev_q+model.si(s).qn(n)-model.si(s).qn(u) 
                        : prev_q-model.si(s).qn(n));

            //For the last site, only keep blocks 
            //compatible with specified totalq i.e. q=0 here
            if(S == Send && q != QN()) continue;

            //Set Site indices of A[s] and its previous Link Index
            block = A[s];
            if(S != start) block *= conj(comp);
            block *= Index(model.si(s))(n);
            if(is_mpo) block *= Index(model.siP(s))(u);

            //Initialize D Vector (D records which values of
            //the right Link Index to keep for the current QN q)
            int count = qD.count(q);
            Vector& D = qD[q];
            if(count == 0) { D.ReDimension(bond.m()); D = 0; }

            if(s == show_s)
                {
                cerr << boost::format("For n = %d\n")%n;
                cerr << boost::format("Got a block with norm %.10f\n")%block.norm();
                cerr << boost::format("bond.m() = %d\n")%bond.m();
                PrintDat(block);
                if(s != 1) PrintDat(comp);
                }

            bool keep_block = false;
            if(S == Send) 
                { keep_block = true; }
            else
                {
                if(bond.m() == 1 && block.norm() != 0) 
                    { D = 1; keep_block = true; }
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
                        Index new_ind = (comp.index(1)==prev_bond ? comp.index(2) : comp.index(1));
                        summed_block = ITensor(new_ind,1) * block;
                        }
                    //cerr << boost::format("s = %d, bond=")%s << bond << "\n";
                    //summed_block.print("summed_block");

                    Real rel_cut = -1;
                    for(int j = 1; j <= bond.m(); ++j)
                        { rel_cut = max(fabs(summed_block.val1(j)),rel_cut); }
                    assert(rel_cut >= 0);
                    //Real rel_cut = summed_block.norm()/summed_block.vecSize();
                    rel_cut *= cut;
                    //cerr << "rel_cut == " << rel_cut << "\n";

                    if(rel_cut > 0)
                    for(int j = 1; j <= bond.m(); ++j)
                        {
                        if(fabs(summed_block.val1(j)) > rel_cut) 
                            { D(j) = 1; keep_block = true; }
                        }
                    }
                } //else (S != Send)

            if(keep_block)
                {
                qD[q] = D;

                if(is_mpo) 
                    {
                    block.addindex1(conj(model.si(s)(n).index()));
                    block.addindex1(model.siP(s)(u).index());
                    }
                else 
                    { block.addindex1(model.si(s)(n).index()); }

                qt[q].push_back(block);

                if(s==show_s)
                    {
                    block.print("Kept block",ShowData);
                    cerr << "D = " << D << "\n";
                    }
                }
            }

        qC.clear();

        Foreach(const qt_vt& x, qt)
            {
            const vector<ITensor>& blks = x.second;
            if(blks.size() != 0)
                {
                q = x.first; 
                if(S == Send) 
                    { Foreach(const ITensor& t, blks) nblock.push_back(t); }
                else
                    {
                    Matrix M; 
                    int mm = collapseCols(qD[q],M);
                    if(s==show_s)
                        {
                        cerr << boost::format("Adding block, mm = %d\n")%mm;
                        q.print("q");
                        cerr << "qD[q] = " << qD[q] << "\n";
                        cerr << "M = \n" << M << "\n";
                        int count = 0;
                        Foreach(const ITensor& t, blks) 
                        t.print((boost::format("t%02d")%(++count)).str(),ShowData);
                        }
                    //string qname = (boost::format("ql%d(%+d:%d:%s)")%s%q.sz()%q.Nf()%(q.Nfp() == 0 ? "+" : "-")).str();
                    string qname = (boost::format("ql%d(%+d:%d)")%s%q.sz()%q.Nf()).str();
                    Index qbond(qname,mm);
                    ITensor compressor(bond,qbond,M);
                    Foreach(const ITensor& t, blks) nblock.push_back(t * compressor);
                    iq.push_back(inqn(qbond,q));
                    qC[q] = compressor;
                    }
                }
            }

        if(S != Send) 
            { 
            if(iq.empty()) 
                {
                cerr << "At site " << s << "\n";
                Error("convertToIQ: no compatible QNs to put into Link.");
                }
            linkind[s] = IQIndex(nameint("qL",s),iq); iq.clear(); 
            }
        if(S == start)
            {
            qA.at(s) = (is_mpo ? IQTensor(conj(model.si(s)),model.siP(s),linkind.at(s)) 
                            : IQTensor(model.si(s),linkind[s]));
            }
        else 
        if(S == Send)
            {
            qA.at(s) = (is_mpo ? IQTensor(conj(linkind[sprev]),conj(model.si(s)),model.siP(s)) 
                            : IQTensor(conj(linkind[sprev]),model.si(s)));
            }
        else
            {
            qA.at(s) = (is_mpo ? IQTensor(conj(linkind[sprev]),conj(model.si(s)),model.siP(s),linkind[s]) 
                            : IQTensor(conj(linkind[sprev]),model.si(s),linkind[s]));
            }

        Foreach(const ITensor& nb, nblock) 
            { qA.at(s) += nb; } 
        nblock.clear();

        if(s==show_s)
            {
            qA[s].print((boost::format("qA[%d]")%s).str(),ShowData);
            Error("Stopping");
            }

        } //for loop over s

    } //void convertToIQ

/*
template <class Tensor> 
template <class IQMPSType> 
void MPSt<Tensor>::convertToIQ(IQMPSType& iqpsi, QN totalq, Real cut) const
{
    assert(model_ != 0);
    const Model& sst = *model_;

    iqpsi = IQMPSType(sst,maxm,cutoff);

    if(!A[1].hasindex(si(1))) Error("convertToIQ: incorrect primelevel for conversion");
    bool is_mpo = A[1].hasindex(primed(si(1)));
    const int Dim = si(1).m();
    const int PDim = (is_mpo ? Dim : 1);

    vector<IQIndex> linkind(N);

    typedef map<QN,Vector>::value_type qD_vt;
    map<QN,Vector> qD; //Diags of compressor matrices by QN
    typedef map<QN,vector<ITensor> >::value_type qt_vt;
    map<QN,vector<ITensor> > qt; //ITensor blocks by QN
    typedef map<QN,ITensor>::value_type qC_vt;
    map<QN,ITensor> qC; //Compressor ITensors by QN
    ITensor block;
    vector<ITensor> nblock;
    vector<inqn> iq;

    QN q;

    qC[totalq] = ITensor(); //Represents Virtual index
    //First value of prev_q below set to totalq

    const int show_s = 0;

    Index bond, prev_bond;
    for(int s = 1; s <= N; ++s)
    {

        qD.clear(); qt.clear();
        if(s > 1) prev_bond = LinkInd(s-1); 
        if(s < N) bond = LinkInd(s);

        if(s == show_s) 
        {
            PrintDat(A[s]);
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

            //Set Site indices of A[s] and its previous Link Index
            block = A[s];
            if(s != 1) block *= conj(comp);
            block *= si(s)(n);
            if(is_mpo) block *= siP(s)(u);

            //Initialize D Vector (D records which values of
            //the right Link Index to keep for the current QN q)
            int count = qD.count(q);
            Vector& D = qD[q];
            if(count == 0) { D.ReDimension(bond.m()); D = 0; }

            if(s == show_s)
            {
                cerr << boost::format("For n = %d\n")%n;
                cerr << boost::format("Got a block with norm %.10f\n")%block.norm();
                cerr << boost::format("bond.m() = %d\n")%bond.m();
                PrintDat(block);
                if(s != 1) PrintDat(comp);
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
                    //cerr << boost::format("s = %d, bond=")%s << bond << "\n";
                    //summed_block.print("summed_block");

                    Real rel_cut = -1;
                    for(int j = 1; j <= bond.m(); ++j)
                    { rel_cut = max(fabs(summed_block.val1(j)),rel_cut); }
                    assert(rel_cut >= 0);
                    //Real rel_cut = summed_block.norm()/summed_block.vecSize();
                    rel_cut *= cut;
                    //cerr << "rel_cut == " << rel_cut << "\n";

                    if(rel_cut > 0)
                    for(int j = 1; j <= bond.m(); ++j)
                    if(fabs(summed_block.val1(j)) > rel_cut) 
                    { D(j) = 1; keep_block = true; }
                }
            } //else (s != N)

            //if(!keep_block && q == totalq) { D(1) = 1; keep_block = true; }

            if(keep_block)
            {
                qD[q] = D;

                if(is_mpo) 
                {
                block.addindex1(conj(si(s)(n).index()));
                block.addindex1(siP(s)(u).index());
                }
                else { block.addindex1(si(s)(n).index()); }

                qt[q].push_back(block);

                if(s==show_s)
                {
                block.print("Kept block",ShowData);
                cerr << "D = " << D << "\n";
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
                    Matrix M; int mm = collapseCols(qD[q],M);
                    if(s==show_s)
                    {
                        cerr << boost::format("Adding block, mm = %d\n")%mm;
                        q.print("q");
                        cerr << "qD[q] = " << qD[q] << "\n";
                        cerr << "M = \n" << M << "\n";
                        int count = 0;
                        Foreach(const ITensor& t, blks) 
                        t.print((boost::format("t%02d")%(++count)).str(),ShowData);
                    }
                    //string qname = (boost::format("ql%d(%+d:%d:%s)")%s%q.sz()%q.Nf()%(q.Nfp() == 0 ? "+" : "-")).str();
                    string qname = (boost::format("ql%d(%+d:%d)")%s%q.sz()%q.Nf()).str();
                    Index qbond(qname,mm);
                    ITensor compressor(bond,qbond,M);
                    Foreach(const ITensor& t, blks) nblock.push_back(t * compressor);
                    iq.push_back(inqn(qbond,q));
                    qC[q] = compressor;
                }
            }
        }

        if(s != N) 
        { 
            if(iq.empty()) 
            {
                cerr << "At site " << s << "\n";
                Error("convertToIQ: no compatible QNs to put into Link.");
            }
            linkind[s] = IQIndex(nameint("qL",s),iq); iq.clear(); 
        }
        if(s == 1)
        {
            iqpsi.AAnc(s) = (is_mpo ? IQTensor(conj(si(s)),siP(s),linkind[s]) : IQTensor(si(s),linkind[s]));
        }
        else if(s == N)
        {
            iqpsi.AAnc(s) = (is_mpo ? IQTensor(conj(linkind[s-1]),conj(si(s)),siP(s)) 
                                    : IQTensor(conj(linkind[s-1]),si(s)));
        }
        else
        {
            iqpsi.AAnc(s) = (is_mpo ? IQTensor(conj(linkind[s-1]),conj(si(s)),siP(s),linkind[s]) 
                                    : IQTensor(conj(linkind[s-1]),si(s),linkind[s]));
        }

        Foreach(const ITensor& nb, nblock) { iqpsi.AAnc(s) += nb; } nblock.clear();

        if(0) //try to get this working ideally
        if(!is_mpo && s > 1) 
        {
            IQTensor AA = iqpsi.bondTensor(s-1);
            iqpsi.doSVD(s-1,AA,Fromleft);
        }

        if(s==show_s)
        {
        iqpsi.AA(s).print((boost::format("qA[%d]")%s).str(),ShowData);
        Error("Stopping");
        }

    } //for loop over s

    assert(checkQNs(iqpsi));

} //void convertToIQ(IQMPSType& iqpsi) const
*/

int 
findCenter(const IQMPS& psi)
    {
    for(int j = 1; j <= psi.NN(); ++j) 
        {
        const IQTensor& A = psi.AA(j);
        if(A.r() == 0) Error("Zero rank tensor in MPS");
        bool allSameDir = true;
        Arrow dir = A.index(1).dir();
        for(int i = 2; i <= A.r(); ++i)
            {
            if(A.index(i).type() == ReIm) continue;
            if(A.index(i).dir() != dir)
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
    const int N = psi.NN();

    QN Zero;

    int center = findCenter(psi);
    if(center == -1)
        {
        std::cerr << "Did not find an ortho. center\n";
        return false;
        }

    //Check that all IQTensors have zero div
    //except possibly the ortho. center
    for(int i = 1; i <= N; ++i) 
        {
        if(i == center) continue;
        if(psi.AA(i).isNull())
            {
            std::cerr << boost::format("AA(%d) null, QNs not well defined\n")%i;
            return false;
            }
        if(psi.AA(i).div() != Zero)
            {
            std::cerr << "At i = " << i << "\n";
            Print(psi.AA(i));
            std::cerr << "IQTensor other than the ortho center had non-zero divergence\n";
            return false;
            }
        }

    //Check arrows from left edge
    for(int i = 1; i < center; ++i)
        {
        if(psi.RightLinkInd(i).dir() != In) 
            {
            std::cerr << boost::format("checkQNs: At site %d to the left of the OC, Right side Link not pointing In\n")%i;
            return false;
            }
        if(i > 1)
            {
            if(psi.LeftLinkInd(i).dir() != Out) 
                {
                std::cerr << boost::format("checkQNs: At site %d to the left of the OC, Left side Link not pointing Out\n")%i;
                return false;
                }
            }
        }

    //Check arrows from right edge
    for(int i = N; i > center; --i)
        {
        if(i < N)
        if(psi.RightLinkInd(i).dir() != Out) 
            {
            std::cerr << boost::format("checkQNs: At site %d to the right of the OC, Right side Link not pointing Out\n")%i;
            return false;
            }
        if(psi.LeftLinkInd(i).dir() != In) 
            {
            std::cerr << boost::format("checkQNs: At site %d to the right of the OC, Left side Link not pointing In\n")%i;
            return false;
            }
        }

    //Done checking arrows
    return true;
    }

void 
fitWF(const IQMPS& psi_basis, IQMPS& psi_to_fit)
    {
    if(!psi_basis.is_ortho()) Error("psi_basis must be orthogonolized.");
    if(psi_basis.ortho_center() != 1) Error("psi_basis must be orthogonolized to site 1.");


    int N = psi_basis.NN();
    if(psi_to_fit.NN() != N) Error("Wavefunctions must have same number of sites.");

    IQTensor A = psi_to_fit.AA(N) * conj(primelink(psi_basis.AA(N)));
    for(int n = N-1; n > 1; --n)
        {
        A *= conj(primelink(psi_basis.AA(n)));
        A *= psi_to_fit.AA(n);
        }
    A = psi_to_fit.AA(1) * A;
    A.noprime();

    psi_to_fit = psi_basis;
    psi_to_fit.AAnc(1) = A;
    }
