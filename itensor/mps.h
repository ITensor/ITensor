#ifndef __MPS_H
#define __MPS_H
#include "svdworker.h"
#include "model.h"

static const LogNumber DefaultRefScale(7.58273202392352185);

void convertToIQ(const BaseModel& model, const std::vector<ITensor>& A, std::vector<IQTensor>& qA, QN totalq = QN(), Real cut = 1E-12);

template<class Tensor, class TensorSet>
Real doDavidson(Tensor& phi, const TensorSet& mpoh, const TensorSet& LH, const TensorSet& RH, int niter, int debuglevel, Real errgoal);

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

class InitState
{
    int N;
    std::vector<IQIndexVal> state;
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
    operator std::vector<IQIndexVal>() const { return state; }
};

template <class Tensor>
class MPSt //the lowercase t stands for "type" or "template"
{
public:

    typedef Tensor TensorT;
    typedef typename Tensor::IndexT IndexT;
    typedef typename Tensor::IndexValT IndexValT;
    typedef BaseModel ModelT;

    //Accessor Methods ------------------------------

    int 
    NN() const { return N;}

    int 
    right_lim() const { return right_orth_lim; }

    int 
    left_lim() const { return left_orth_lim; }

    IQIndex 
    si(int i) const { return model_->si(i); }

    IQIndex 
    siP(int i) const { return model_->siP(i); }

    typedef typename std::vector<Tensor>::const_iterator AA_it;
    const std::pair<AA_it,AA_it> 
    AA() const { return std::make_pair(A.begin()+1,A.end()); }

    const Tensor& 
    AA(int i) const { return GET(A,i); }

    Tensor& 
    AAnc(int i) //nc means 'non const'
        { 
        if(i <= left_orth_lim) left_orth_lim = i-1;
        if(i >= right_orth_lim) right_orth_lim = i+1;
        return GET(A,i); 
        }

    const ModelT& 
    model() const { return *model_; }

    const SVDWorker& 
    svd() const { return svd_; }


    Tensor& 
    setU(int i, Direction dir); //set unitary

    bool 
    is_null() const { return (model_==0); }
    bool 
    is_not_null() const { return (model_!=0); }

    bool 
    doRelCutoff() const { return svd_.doRelCutoff(); }
    void 
    doRelCutoff(bool val) { svd_.doRelCutoff(val); }

    bool 
    absoluteCutoff() const { return svd_.absoluteCutoff(); }
    void 
    absoluteCutoff(bool val) { svd_.absoluteCutoff(val); }

    LogNumber 
    refNorm() const { return svd_.refNorm(); }
    void 
    refNorm(LogNumber val) { svd_.refNorm(val); }

    Real 
    cutoff() const { return svd_.cutoff(); }
    void 
    cutoff(Real val) { svd_.cutoff(val); }

    int 
    minm() const { return svd_.minm(); }
    void 
    minm(int val) { svd_.minm(val); }

    int 
    maxm() const { return svd_.maxm(); }
    void 
    maxm(int val) { svd_.maxm(val); }

    Real 
    truncerr(int b) const { return svd_.truncerr(b); }

    const Vector& 
    eigsKept(int b) const { return svd_.eigsKept(b); }

    bool 
    showeigs() const { return svd_.showeigs(); }
    void 
    showeigs(bool val) { svd_.showeigs(val); }

    Tensor 
    bondTensor(int b) const 
        { Tensor res = A.at(b) * A.at(b+1); return res; }

    //MPSt: Constructors --------------------------------------------

    MPSt() 
        : N(0), model_(0)
        { }

    MPSt(const ModelT& mod_,int maxmm = MAX_M, Real cut = MIN_CUT) 
    : N(mod_.NN()), A(mod_.NN()+1),left_orth_lim(0),right_orth_lim(mod_.NN()),
    model_(&mod_), svd_(N,cut,1,maxmm,false,LogNumber(1))
        { 
        random_tensors(A);
        }

    MPSt(const ModelT& mod_,const InitState& initState,int maxmm = MAX_M, Real cut = MIN_CUT) 
    : N(mod_.NN()),A(mod_.NN()+1),left_orth_lim(0),right_orth_lim(2),
    model_(&mod_), svd_(N,cut,1,maxmm,false,LogNumber(1))
        { 
        init_tensors(A,initState);
        }

    MPSt(const ModelT& model, std::istream& s)
        : N(model.NN()), A(model.NN()+1), model_(&model)
        { 
        read(s); 
        }

    virtual ~MPSt() { }

    void 
    read(std::istream& s);

    void 
    write(std::ostream& s) const;

    //MPSt: operators ------------------------------------------------------

    inline MPSt& 
    operator*=(Real a) { AAnc(left_orth_lim+1) *= a; return *this; }

    inline MPSt 
    operator*(Real r) const { MPSt res(*this); res *= r; return res; }

    friend inline MPSt 
    operator*(Real r, MPSt res) { res *= r; return res; }

    MPSt& 
    operator+=(const MPSt& oth);

    inline MPSt 
    operator+(MPSt res) const { res += *this; return res; }

    inline MPSt 
    operator-(MPSt res) const { res *= -1; res += *this; return res; }

    //MPSt: index methods --------------------------------------------------

    void mapprime(int oldp, int newp, PrimeType pt = primeBoth)
	{ for(int i = 1; i <= N; ++i) A[i].mapprime(oldp,newp,pt); }

    void primelinks(int oldp, int newp)
	{ for(int i = 1; i <= N; ++i) A[i].mapprime(oldp,newp,primeLink); }

    void noprimelink()
	{ for(int i = 1; i <= N; ++i) A[i].noprime(primeLink); }

    IndexT LinkInd(int b) const 
        { return index_in_common(AA(b),AA(b+1),Link); }
    IndexT RightLinkInd(int i) const 
        { assert(i < N); return index_in_common(AA(i),AA(i+1),Link); }
    IndexT LeftLinkInd(int i)  const 
        { assert(i > 1); return index_in_common(AA(i),AA(i-1),Link); }

    //MPSt: orthogonalization methods -------------------------------------

    void doSVD(int b, const Tensor& AA, Direction dir, 
               bool preserve_shape = false)
	{
        assert(b > 0);
        assert(b < N);

        if(dir == Fromleft && b-1 > left_orth_lim)
        {
            std::cerr << boost::format("b=%d, left_orth_lim=%d\n")
                    %b%left_orth_lim;
            Error("b-1 > left_orth_lim");
        }
        if(dir == Fromright && b+2 < right_orth_lim)
        {
            std::cerr << boost::format("b=%d, right_orth_lim=%d\n")
                    %b%right_orth_lim;
            Error("b+2 < right_orth_lim");
        }

        svd_(b,AA,A[b],A[b+1],dir);
                 
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

    void orthogonalize(bool verbose = false)
    {
        //Do a half-sweep to the right, orthogonalizing each bond
        //but do not truncate since the basis to the right might not
        //be ortho (i.e. use the current m).
        svd_.truncate(false);
        position(N);
        if(verbose)
            std::cout << "Done orthogonalizing, starting truncation." << std::endl;
        //Now basis is ortho, ok to truncate
        svd_.truncate(true);
        position(1);
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
        std::cerr << "checkOrtho: on line " << __LINE__ << " of mps.h," << std::endl;
        std::cerr << "checkOrtho: Tensor at position " << i << " failed to be " << (left ? "left" : "right") << " ortho." << std::endl;
        std::cerr << "checkOrtho: Norm(diff) = " << boost::format("%E") % Norm(diff) << std::endl;
        std::cerr << "checkOrtho: Error threshold set to " << boost::format("%E") % threshold << std::endl;
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
            std::cerr << "checkOrtho: A[i] not left orthogonal at site i=" << i << std::endl;
            return false;
        }

        for(int i = NN(); i >= right_orth_lim; --i)
        if(!checkRightOrtho(i))
        {
            std::cerr << "checkOrtho: A[i] not right orthogonal at site i=" << i << std::endl;
            return false;
        }
        return true;
    }

    void getCenter(int j, Direction dir, Tensor& lambda, bool do_signfix = false)
    {
        getCenterMatrix(AAnc(j),(dir == Fromleft ? RightLinkInd(j) : LeftLinkInd(j)),cutoff,minm,maxm,lambda,nameint("c",j));

        if(dir == Fromleft) {
		if(left_orth_lim == j-1 || j == 1) {
		       	left_orth_lim = j;
		}
	}
        else if(right_orth_lim == j+1 || j == N) right_orth_lim = j;

        if(do_signfix) Error("do_signfix not implemented.");
    }

    template<class TensorSet>
    Real bondDavidson(int b, const TensorSet& mpoh, const TensorSet& LH, const TensorSet& RH, 
    int niter, int debuglevel, Direction dir, Real errgoal=1E-4)
    {
        if(b-1 > left_orth_lim)
            {
            std::cerr << boost::format("b=%d, Lb=%d\n")%b%left_orth_lim;
            Error("b-1 > left_orth_lim");
            }
        if(b+2 < right_orth_lim)
            {
            std::cerr << boost::format("b+1=%d, Rb=%d\n")%(b+1)%right_orth_lim;
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
            { 
            std::cerr << boost::format("projectOp: from left j > left_orth_lim (j=%d,left_orth_lim=%d)\n")%j%left_orth_lim; 
            Error("Projecting operator at j > left_orth_lim"); 
            }
        if(dir==Fromright && j < right_orth_lim) 
            { 
            std::cerr << boost::format("projectOp: from left j < right_orth_lim (j=%d,right_orth_lim=%d)\n")%j%right_orth_lim; 
            Error("Projecting operator at j < right_orth_lim"); 
            }
        res = (P.is_null() ? AA(j) : P * AA(j));
        res *= Op; res *= conj(primed(AA(j)));
    }


    void applygate(const Tensor& gate)
	{
        Tensor AA = A[left_orth_lim+1] * A[left_orth_lim+2] * gate;
        AA.noprime();
        doSVD(left_orth_lim+1,AA,Fromleft);
	}

    Real 
    norm() const { return sqrt(psiphi(*this,*this)); }

    int
    averageM() const;

    Real normalize()
    {
        Real norm_ = norm();
        if(fabs(norm_) < 1E-20) Error("Zero norm");
        *this *= 1.0/norm_;
        return norm_;
    }

    bool is_complex() const
    { return A[left_orth_lim+1].is_complex(); }

    friend inline std::ostream& 
    operator<<(std::ostream& s, const MPSt& M)
        {
        s << "\n";
        for(int i = 1; i <= M.NN(); ++i) s << M.AA(i) << "\n";
        return s;
        }

    void print(std::string name = "",Printdat pdat = HideData) const 
        { 
        bool savep = Globals::printdat();
        Globals::printdat() = (pdat==ShowData); 
        std::cerr << "\n" << name << " =\n" << *this << "\n"; 
        Globals::printdat() = savep;
        }

    void 
    printIndices(const std::string& name = "") const
        {
        std::cout << name << "=\n";
        for(int i = 1; i <= NN(); ++i) 
            AA(i).printIndices(boost::format("AA(%d)")%i);
        }

    void 
    printIndices(const boost::format& fname) const
        { printIndices(fname.str()); }

    void toIQ(QN totalq, MPSt<IQTensor>& iqpsi, Real cut = 1E-12) const
    {
        iqpsi = MPSt<IQTensor>(*model_,maxm(),cutoff());
        iqpsi.svd_ = svd_;
        convertToIQ(*model_,A,iqpsi.A,totalq,cut);
    }

protected:

    //Data Members ----------

    int N;

    std::vector<Tensor> A;

    int left_orth_lim,
        right_orth_lim;

    const ModelT* model_;

    SVDWorker svd_;

    //Constructor Helpers ----------

    void 
    new_tensors(std::vector<ITensor>& A_);

    void 
    random_tensors(std::vector<ITensor>& A_);

    void 
    random_tensors(std::vector<IQTensor>& A_) { }

    void 
    init_tensors(std::vector<ITensor>& A_, const InitState& initState)
        { 
        new_tensors(A_); 
        for(int i = 1; i <= N; ++i) 
            {
            A_[i] = ITensor(initState(i)); 
            }
        }

    void 
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

private:
    friend class MPSt<ITensor>;
    friend class MPSt<IQTensor>;
}; //class MPSt<Tensor>
typedef MPSt<ITensor> MPS;
typedef MPSt<IQTensor> IQMPS;

int 
findCenter(const IQMPS& psi);

inline bool 
checkQNs(const MPS& psi) { return true; }

bool 
checkQNs(const IQMPS& psi);

inline QN 
totalQN(const IQMPS& psi)
    {
    int center = findCenter(psi);
    if(center == -1)
        Error("Could not find ortho. center");
    return psi.AA(center).div();
    }

template <class MPSType>
void 
psiphi(const MPSType& psi, const MPSType& phi, Real& re, Real& im)  // <psi | phi>
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
Real 
psiphi(const MPSType& psi, const MPSType& phi) //Re[<psi|phi>]
    {
    Real re, im;
    psiphi(psi,phi,re,im);
    if(im != 0) 
	if(fabs(im) > 1.0e-12 * fabs(re))
	    std::cerr << "Real psiphi: WARNING, dropping non-zero imaginary part of expectation value.\n";
    return re;
    }

//Computes an MPS which has the same overlap with psi_basis as psi_to_fit,
//but which differs from psi_basis only on the first site, and has same index
//structure as psi_basis. Result is stored to psi_to_fit on return.
void 
fitWF(const IQMPS& psi_basis, IQMPS& psi_to_fit);

//Template method for efficiently summing a set of MPS's or MPO's (or any class supporting operator+=)
template <typename MPSType>
void 
sum(const std::vector<MPSType>& terms, MPSType& res, Real cut = MIN_CUT, int maxm = MAX_M)
{
    int Nt = terms.size();
    if(Nt == 1) 
        {
        res = terms[0];
        res.cutoff(cut); res.maxm(maxm);
        return;
        }
    else if(Nt == 2)
        { 
        res = terms[0];
        res.cutoff(cut); res.maxm(maxm);
        //std::cerr << boost::format("Before +=, cutoff = %.1E, maxm = %d\n")%(res.cutoff)%(res.maxm);
        res += terms[1];
        return;
        }
    else if(Nt > 2)
        {
        //Add all MPS's in pairs
        std::vector<MPSType> terms2(2), nterms; nterms.reserve(Nt/2);
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
} // void sum(const std::vector<MPSType>& terms, Real cut, int maxm, MPSType& res)


#endif
