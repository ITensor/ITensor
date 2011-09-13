#ifndef __MPS_H
#define __MPS_H
#include "iqcombiner.h"
#include "model.h"

enum Direction { Fromright, Fromleft, Both, None };

static const LogNumber DefaultRefScale(7.58273202392352185);
extern bool showeigs;

void convertToIQ(const BaseModel& model, const vector<ITensor>& A, vector<IQTensor>& qA, QN totalq = QN(), Real cut = 1E-12);

template<class Tensor, class TensorSet>
Real doDavidson(Tensor& phi, const TensorSet& mpoh, const TensorSet& LH, const TensorSet& RH, int niter, int debuglevel, Real errgoal);

extern Real truncerror, svdtruncerr;

template<class Tensor, class IndexT>
IndexT index_in_common(const Tensor& A, const Tensor& B, IndexType t)
{
    for(int j = 1; j <= A.r(); ++j)
    {
        const IndexT& I = A.index(j);
        if(I.type() == t && B.hasindex(I)) { return I; }
    }
    return IndexT();
}
inline Index index_in_common(const ITensor& A, const ITensor& B, IndexType t)
{ return index_in_common<ITensor,Index>(A,B,t); }
inline IQIndex index_in_common(const IQTensor& A, const IQTensor& B, IndexType t)
{ return index_in_common<IQTensor,IQIndex>(A,B,t); }

namespace {
int collapseCols(const Vector& Diag, Matrix& M)
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
}


class InitState
{
    int N;
    vector<IQIndexVal> state;
    typedef IQIndexVal (*SetFuncPtr)(int);
public:
    int NN() const { return N; }

    InitState(int nsite) : N(nsite), state(N+1) { }
    InitState(int nsite,SetFuncPtr setter) : N(nsite), state(N+1) 
    { set_all(setter); }

    void set_all(SetFuncPtr setter)
    { for(int j = 1; j <= N; ++j) GET(state,j-1) = (*setter)(j); }

    IQIndexVal& operator()(int i) { return GET(state,i-1); }
    const IQIndexVal& operator()(int i) const { return GET(state,i-1); }
    operator vector<IQIndexVal>() const { return state; }
};

template <class Tensor>
class MPSt //the lowercase t stands for "type" or "template"
{
public:
    typedef Tensor TensorT;
    typedef typename Tensor::IndexT IndexT;
    typedef typename Tensor::IndexValT IndexValT;
    typedef BaseModel ModelT;
protected:
    int N;
    vector<Tensor> A;
    int left_orth_lim,right_orth_lim;
    const ModelT* model_;
    bool doRelCutoff_;
    LogNumber refNorm_;

    void new_tensors(vector<ITensor>& A_)
    {
        vector<Index> a(N+1);
        for(int i = 1; i <= N; ++i)
        { a[i] = Index(nameint("a",i)); }
        A_[1] = ITensor(si(1),a[1]);
        for(int i = 2; i < N; i++)
        { A_[i] = ITensor(conj(a[i-1]),si(i),a[i]); }
        A_[N] = ITensor(conj(a[N-1]),si(N));
    }

    void random_tensors(vector<ITensor>& A_)
    { new_tensors(A_); for(int i = 1; i <= N; ++i) A_[i].Randomize(); }

    void random_tensors(vector<IQTensor>& A_) { }

    void init_tensors(vector<ITensor>& A_, const InitState& initState)
    { new_tensors(A_); for(int i = 1; i <= N; ++i) A_[i](initState(i))=1; }

    void init_tensors(vector<IQTensor>& A_, const InitState& initState)
    {
        vector<QN> qa(N+1); //qn[i] = qn on i^th bond
        for(int i = 1; i <= N; ++i) { qa[0] -= initState(i).qn()*In; }
        IQIndex Center("Center",Index("center",1,Virtual),qa[0],In);

        //Taking OC to be at the leftmost site,
        //compute the QuantumNumbers of all the Links.
        for(int i = 1; i <= N; ++i)
        {
            //Taking the divergence to be zero,solve for qa[i]
            qa[i] = Out*(-qa[i-1]*In - initState(i).qn());
        }

        vector<IQIndex> a(N+1);
        for(int i = 1; i <= N; ++i)
        { a[i] = IQIndex(nameint("L",i),Index(nameint("l",i)),qa[i]); }
        A_[1] = IQTensor(si(1),a[1]); A_[1].addindex1(Center); A_[1](initState(1))=1;
        for(int i = 2; i < N; ++i)
        { 
        A_[i] = IQTensor(conj(a[i-1]),si(i),a[i]); 
        A_[i](initState(i))=1;
        }
        A_[N] = IQTensor(conj(a[N-1]),si(N)); A_[N](initState(N))=1;
    }

    typedef pair<typename vector<Tensor>::const_iterator,typename vector<Tensor>::const_iterator> const_range_type;
public:
    int minm,maxm;
    Real cutoff;

    //Accessor Methods ------------------------------

    int NN() const { return N;}
    int right_lim() const { return right_orth_lim; }
    int left_lim() const { return left_orth_lim; }
    IQIndex si(int i) const { return model_->si(i); }
    IQIndex siP(int i) const { return model_->siP(i); }
    typedef typename vector<Tensor>::const_iterator AA_it;
    const pair<AA_it,AA_it> AA() const { return make_pair(A.begin()+1,A.end()); }
    const Tensor& AA(int i) const { return GET(A,i); }
    const ModelT& model() const { return *model_; }
    Tensor& AAnc(int i) //nc means 'non const'
    { 
        if(i <= left_orth_lim) left_orth_lim = i-1;
        if(i >= right_orth_lim) right_orth_lim = i+1;
        return GET(A,i); 
    }
    Tensor& setU(int i, Direction dir) //set unitary
    {
        if(dir == Fromleft)
        {
        if(i == left_orth_lim+1) left_orth_lim = i;
        else Error("left_orth_lim not at i-1");
        }
        else if(dir == Fromright)
        {
        if(i == right_orth_lim-1) right_orth_lim = i;
        else Error("right_orth_lim not at i-1");
        }
        return GET(A,i);
    }
    bool is_null() const { return (model_==0); }
    bool is_not_null() const { return (model_!=0); }

    bool doRelCutoff() const { return doRelCutoff_; }
    void doRelCutoff(bool val) 
        { doRelCutoff_ = val; }
    LogNumber refNorm() const { return refNorm_; }
    void refNorm(LogNumber val) 
        { if(val == 0) { Error("zero refNorm"); } refNorm_ = val; }

    Tensor bondTensor(int b) const 
        { Tensor res = A.at(b) * A.at(b+1); return res; }

    //MPSt: Constructors --------------------------------------------

    MPSt() 
    : N(0), model_(0), doRelCutoff_(false), refNorm_(1), 
    minm(1), maxm(MAX_M), cutoff(MAX_CUT)
    { }

    MPSt(const ModelT& mod_,int maxmm = MAX_M, Real cut = MAX_CUT) 
    : N(mod_.NN()), A(mod_.NN()+1),left_orth_lim(0),right_orth_lim(mod_.NN()),
    model_(&mod_), doRelCutoff_(false), refNorm_(1), 
    minm(1), maxm(maxmm), cutoff(cut)
	{ random_tensors(A); }

    MPSt(const ModelT& mod_,const InitState& initState,int maxmm = MAX_M, Real cut = MAX_CUT) 
    : N(mod_.NN()),A(mod_.NN()+1),left_orth_lim(0),right_orth_lim(2),
    model_(&mod_), doRelCutoff_(false), refNorm_(1), 
    minm(1), maxm(maxmm), cutoff(cut)
	{ init_tensors(A,initState); }

    MPSt(const ModelT& model, istream& s) { read(model,s); }

    virtual ~MPSt() { }

    void read(const ModelT& model, istream& s)
    {
        model_ = &model;
        N = model_->NN();
        A.resize(N);
        for(int j = 1; j <= N; ++j) A[j].read(s);
        s.read((char*) &left_orth_lim,sizeof(left_orth_lim));
        s.read((char*) &right_orth_lim,sizeof(right_orth_lim));
        s.read((char*) &minm,sizeof(minm));
        s.read((char*) &maxm,sizeof(maxm));
        s.read((char*) &cutoff,sizeof(cutoff));
        s.read((char*) &doRelCutoff_,sizeof(doRelCutoff_));
        s.read((char*) &refNorm_,sizeof(refNorm_));
    }

    void write(ostream& s) const
    {
        for(int j = 1; j <= N; ++j) A[j].write(s);
        s.write((char*) &left_orth_lim,sizeof(left_orth_lim));
        s.write((char*) &right_orth_lim,sizeof(right_orth_lim));
        s.write((char*) &minm,sizeof(minm));
        s.write((char*) &maxm,sizeof(maxm));
        s.write((char*) &cutoff,sizeof(cutoff));
        s.write((char*) &doRelCutoff_,sizeof(doRelCutoff_));
        s.write((char*) &refNorm_,sizeof(refNorm_));
    }


    //MPSt: operators ------------------------------------------------------

    MPSt& operator*=(Real a)
    {
        AAnc(left_orth_lim+1) *= a;
        return *this;
    }
    inline MPSt operator*(Real r) const { MPSt res(*this); res *= r; return res; }
    friend inline MPSt operator*(Real r, MPSt res) { res *= r; return res; }

    MPSt& operator+=(const MPSt& oth);
    inline MPSt operator+(MPSt res) const { res += *this; return res; }
    //friend inline MPSt operator+(MPSt A, const MPSt& B) { A += B; return A; }
    inline MPSt operator-(MPSt res) const { res *= -1; res += *this; return res; }

    //MPSt: index methods --------------------------------------------------

    void mapprime(int oldp, int newp, PrimeType pt = primeBoth)
	{ for(int i = 1; i <= N; ++i) A[i].mapprime(oldp,newp,pt); }

    void primelinks(int oldp, int newp)
	{ for(int i = 1; i <= N; ++i) A[i].mapprime(oldp,newp,primeLink); }

    void noprimelink()
	{ for(int i = 1; i <= N; ++i) A[i].noprime(primeLink); }

    IndexT LinkInd(int i) const { return index_in_common(A[i],A[i+1],Link); }
    IndexT RightLinkInd(int i) const { assert(i<NN()); return index_in_common(AA(i),AA(i+1),Link); }
    IndexT LeftLinkInd(int i)  const { assert(i>1); return index_in_common(AA(i),AA(i-1),Link); }

    //MPSt: orthogonalization methods -------------------------------------

    virtual void doSVD(int b, const Tensor& AA, Direction dir, bool preserve_shape = false)
	{
        assert(b > 0);
        assert(b < N);

        if(dir == Fromleft && b-1 > left_orth_lim)
        {
            cerr << boost::format("b=%d, left_orth_lim=%d\n")%b%left_orth_lim;
            Error("b-1 > left_orth_lim");
        }
        if(dir == Fromright && b+2 < right_orth_lim)
        {
            cerr << boost::format("b=%d, right_orth_lim=%d\n")%b%right_orth_lim;
            Error("b+2 < right_orth_lim");
        }

        tensorSVD(AA,GET(A,b),GET(A,b+1),
                  cutoff,minm,maxm,dir,doRelCutoff_,refNorm_);
        truncerror = svdtruncerr;

        if(dir == Fromleft)
        {
            //if(left_orth_lim >= b-1 || b == 1) 
            left_orth_lim = b;
            if(right_orth_lim < b+2) right_orth_lim = b+2;
        }
        else //dir == Fromright
        {
            if(left_orth_lim > b-1) left_orth_lim = b-1;
            //if(right_orth_lim <= b+2 || b == N-1) 
            right_orth_lim = b+1;
        }
	}

    //Move the orthogonality center to site i (left_orth_lim = i-1, right_orth_lim = i+1)
    void position(int i, bool preserve_shape = false)
	{
        if(is_null()) Error("position: MPS is null");
        while(left_orth_lim < i-1)
        {
            if(left_orth_lim < 0) left_orth_lim = 0;
            Tensor WF = AA(left_orth_lim+1) * AA(left_orth_lim+2);
            doSVD(left_orth_lim+1,WF,Fromleft,preserve_shape);
        }
        while(right_orth_lim > i+1)
        {
            if(right_orth_lim > N+1) right_orth_lim = N+1;
            Tensor WF = AA(right_orth_lim-2) * AA(right_orth_lim-1);
            doSVD(right_orth_lim-2,WF,Fromright,preserve_shape);
        }
	}

    bool is_ortho() const { return (left_orth_lim + 1 == right_orth_lim - 1); }

    int ortho_center() const 
    { 
        if(!is_ortho()) Error("orthogonality center not well defined.");
        return (left_orth_lim + 1);
    }

    //Checks if A[i] is left (left == true) or right (left == false) orthogonalized
    bool checkOrtho(int i, bool left) const
    {
        IndexT link = (left ? RightLinkInd(i) : LeftLinkInd(i));
        Tensor A = AA(i);
        Tensor Ac = conj(A); Ac.primeind(link,4);

        Tensor rho = A * Ac;

        Tensor Delta(link,link.primed(4),1);

        Tensor Diff = rho - Delta;

        Vector diff(Diff.Length());
        Diff.AssignToVec(diff);

        Real threshold = 1E-13;
        if(Norm(diff) < threshold) return true;

        //Print any helpful debugging info here:
        cerr << "checkOrtho: on line " << __LINE__ << " of mps.h," << endl;
        cerr << "checkOrtho: Tensor at position " << i << " failed to be " << (left ? "left" : "right") << " ortho." << endl;
        cerr << "checkOrtho: Norm(diff) = " << boost::format("%E") % Norm(diff) << endl;
        cerr << "checkOrtho: Error threshold set to " << boost::format("%E") % threshold << endl;
        //-----------------------------

        return false;
    }
    bool checkRightOrtho(int i) const { return checkOrtho(i,false); }
    bool checkLeftOrtho(int i) const { return checkOrtho(i,true); }
    
    bool checkOrtho() const
    {
        for(int i = 1; i <= left_orth_lim; ++i)
        if(!checkLeftOrtho(i))
        {
            cerr << "checkOrtho: A[i] not left orthogonal at site i=" << i << endl;
            return false;
        }

        for(int i = NN(); i >= right_orth_lim; --i)
        if(!checkRightOrtho(i))
        {
            cerr << "checkOrtho: A[i] not right orthogonal at site i=" << i << endl;
            return false;
        }
        return true;
    }

    void getCenter(int j, Direction dir, Tensor& lambda, bool do_signfix = false)
    {
        getCenterMatrix(AAnc(j),(dir == Fromleft ? RightLinkInd(j) : LeftLinkInd(j)),cutoff,minm,maxm,lambda,nameint("c",j));

        if(dir == Fromleft) if(left_orth_lim == j-1 || j == 1) left_orth_lim = j;
        else if(right_orth_lim == j+1 || j == N) right_orth_lim = j;

        if(do_signfix) Error("do_signfix not implemented.");
    }

    template<class TensorSet>
    Real bondDavidson(int b, const TensorSet& mpoh, const TensorSet& LH, const TensorSet& RH, 
    int niter, int debuglevel, Direction dir, Real errgoal=1E-4)
    {
        if(b-1 > left_orth_lim)
        {
            cerr << boost::format("b=%d, Lb=%d\n")%b%left_orth_lim;
            Error("b-1 > left_orth_lim");
        }
        if(b+2 < right_orth_lim)
        {
            cerr << boost::format("b+1=%d, Rb=%d\n")%(b+1)%right_orth_lim;
            Error("b+1 < right_orth_lim");
        }
        Tensor phi = GET(A,b); phi *= GET(A,b+1);
        Real En = doDavidson(phi,mpoh,LH,RH,niter,debuglevel,errgoal);
        doSVD(b,phi,dir);
        return En;
    }

    template<class TensorSet, class OpTensorSet>
    void projectOp(int j, Direction dir, const TensorSet& P, const OpTensorSet& Op, TensorSet& res) const
    {
        if(res.size() != Op.size()) res.resize(Op.size());
        const TensorSet& nP = (P.size() == Op.size() ? P : TensorSet(Op.size()));
        for(unsigned int n = 0; n < Op.size(); ++n) projectOp(j,dir,GET(nP,n),Op[n],GET(res,n));
    }

    template<class OpTensor>
    void projectOp(int j, Direction dir, const Tensor& P, const OpTensor& Op, Tensor& res) const
    {
        if(dir==Fromleft && j > left_orth_lim) 
        { cerr << boost::format("projectOp: from left j > left_orth_lim (j=%d,left_orth_lim=%d)\n")%j%left_orth_lim, Error(""); }
        if(dir==Fromright && j < right_orth_lim) 
        { cerr << boost::format("projectOp: from left j < right_orth_lim (j=%d,right_orth_lim=%d)\n")%j%right_orth_lim, Error(""); }

        res = (P.is_null() ? AA(j) : P * AA(j));
        res *= Op; res *= conj(primed(AA(j)));
    }


    void applygate(const Tensor& gate)
	{
        Tensor AA = A[left_orth_lim+1] * A[left_orth_lim+2] * gate;
        AA.noprime();
        doSVD(left_orth_lim+1,AA,Fromleft);
	}

    Real norm() const { return sqrt(psiphi(*this,*this)); }

    Real normalize()
    {
        Real norm_ = norm();
        if(fabs(norm_) < 1E-20) Error("Zero norm");
        *this *= 1.0/norm_;
        return norm_;
    }

    bool is_complex() const
    { return A[left_orth_lim+1].is_complex(); }

    friend inline ostream& operator<<(ostream& s, const MPSt& M)
    {
        s << "\n";
        for(int i = 1; i <= M.NN(); ++i) s << M.AA(i) << "\n";
        return s;
    }

    void print(string name = "",Printdat pdat = HideData) const 
    { printdat = (pdat==ShowData); cerr << "\n" << name << " =\n" << *this << "\n"; printdat = false; }

    //-----------------------------------------------------------------
    //IQMPS specific methods

    template <class IQMPSType> 
    void convertToIQ(IQMPSType& iqpsi, QN totalq = QN(), Real cut = 1E-12) const;

}; //class MPSt<Tensor>

typedef MPSt<ITensor> MPS;
typedef MPSt<IQTensor> IQMPS;

template<class Tensor>
extern Vector tensorSVD(const Tensor& AA, Tensor& A, Tensor& B, 
                        Real cutoff, int minm, int maxm, 
                        Direction dir, bool doRelCutoff,
                        LogNumber refNorm);

void diag_denmat(const ITensor& rho, Real cutoff, int minm, int maxm, 
                 ITensor& nU, Vector& D, bool doRelCutoff, LogNumber refNorm);

void diag_denmat(const IQTensor& rho, Real cutoff, int minm, int maxm, 
                 IQTensor& nU, Vector& eigs_kept, 
                 bool doRelCutoff, LogNumber refNorm);

inline bool check_QNs(const MPS& psi) { return true; }

inline bool check_QNs(const IQMPS& psi)
{
    const int N = psi.NN();
    //Check Link arrows

    //Get the orthogonality center (based on location of Virtual IQIndex)
    int center = -1;
    for(int i = 1; i <= N; ++i) 
    {
        if(psi.AA(i).hastype(Virtual)) { center = i; break; } 
        //else if(i == N) { cerr << "check_QNs: couldn't find a Virtual IQIndex.\n"; return false; }
    }

    if(center > 0)
    {
        //Check arrows from left edge
        for(int i = 1; i < center; ++i)
        {
            if(psi.RightLinkInd(i).dir() != In) 
            {
                cerr << boost::format("check_QNs: At site %d to the left of the OC, Right side Link not pointing In\n")%i;
                return false;
            }
            if(i > 1)
            {
                if(psi.LeftLinkInd(i).dir() != Out) 
                {
                    cerr << boost::format("check_QNs: At site %d to the left of the OC, Left side Link not pointing Out\n")%i;
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
                cerr << boost::format("check_QNs: At site %d to the right of the OC, Right side Link not pointing Out\n")%i;
                return false;
            }
            if(psi.LeftLinkInd(i).dir() != In) 
            {
                cerr << boost::format("check_QNs: At site %d to the right of the OC, Left side Link not pointing In\n")%i;
                return false;
            }
        }
    }
    //Done checking arrows

    //Check IQTensors
    for(int i = 1; i <= N; ++i)
    {
        if(!check_QNs(psi.AA(i)))
        {
            cerr << "check_QNs: IQTensor AA(" << i << ") had non-zero divergence.\n";
            return false;
        }
    }
    return true;
}

inline QN total_QN(const IQMPS& psi)
{
    assert(psi.AA(psi.ortho_center()).has_virtual());
    return psi.AA(psi.ortho_center()).virtualQN();
}

template <class MPSType>
void psiphi(const MPSType& psi, const MPSType& phi, Real& re, Real& im)  // <psi | phi>
{
    const int N = psi.NN();
    if(N != phi.NN()) Error("psiphi: mismatched N");

    typename MPSType::TensorT L = phi.AA(1) * conj(primeind(psi.AA(1),psi.LinkInd(1))); 
    for(int i = 2; i < psi.NN(); ++i) 
        { L = L * phi.AA(i) * conj(primelink(psi.AA(i))); }
    L = L * phi.AA(N);

    Dot(primeind(psi.AA(N),psi.LinkInd(N-1)),L,re,im);
}
template <class MPSType>
Real psiphi(const MPSType& psi, const MPSType& phi) //Re[<psi|phi>]
{
    Real re, im;
    psiphi(psi,phi,re,im);
    if(im != 0) cerr << "Real psiphi: WARNING, dropping non-zero imaginary part of expectation value.\n";
    return re;
}

//Computes an MPS which has the same overlap with psi_basis as psi_to_fit,
//but which differs from psi_basis only on the first site, and has same index
//structure as psi_basis. Result is stored to psi_to_fit on return.
inline void fitWF(const IQMPS& psi_basis, IQMPS& psi_to_fit)
{
    if(!psi_basis.is_ortho()) Error("fitWF: psi_basis must be orthogonolized.");
    if(psi_basis.ortho_center() != 1) Error("fitWF: psi_basis must be orthogonolized to site 1.");
    psi_to_fit.position(1);

    const IQMPS& psib = psi_basis;
    IQMPS& psif = psi_to_fit;
    IQMPS res = psib;

    int N = psib.NN();
    if(psif.NN() != N) Error("fitWF: sizes of wavefunctions must match.");

    IQTensor A = psif.AA(N) * conj(psib.AA(N));
    for(int n = N-1; n > 1; --n)
    {
        A = conj(psib.AA(n)) * A;
        A = psif.AA(n) * A;
    }
    A = psif.AA(1) * A;

    res.AAnc(1) = A;

    assert(check_QNs(res));

    psi_to_fit = res;
}

//Template method for efficiently summing a set of MPS's or MPO's (or any class supporting operator+=)
template <typename MPSType>
void sum(const vector<MPSType>& terms, MPSType& res, Real cut = MAX_CUT, int maxm = MAX_M)
{
    int Nt = terms.size();
    if(Nt == 1) 
    {
        res = terms[0];
        res.cutoff = cut; res.maxm = maxm;
        return;
    }
    else if(Nt == 2)
	{ 
        res = terms[0];
        res.cutoff = cut; res.maxm = maxm;
        //cerr << boost::format("Before +=, cutoff = %.1E, maxm = %d\n")%(res.cutoff)%(res.maxm);
        res += terms[1];
        return;
    }
    else if(Nt > 2)
	{
        //Add all MPS's in pairs
        vector<MPSType> terms2(2), nterms; nterms.reserve(Nt/2);
        for(int n = 0; n < Nt-1; n += 2)
        {
            terms2[0] = terms[n]; terms2[1] = terms[n+1];
            sum(terms2,res,cut,maxm);
            nterms.push_back(res);
        }
        if(Nt%2 == 1) nterms.push_back(terms.back());
        //Recursively call sum again
        sum(nterms,res,cut,maxm);
        return;
	}
    return;
} // void sum(const vector<MPSType>& terms, Real cut, int maxm, MPSType& res)


#ifdef THIS_IS_MAIN
Real truncerror = 0.0;
Real svdtruncerr = 0.0;
bool showeigs = false;
#endif

#endif
