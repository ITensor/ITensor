#ifndef __MPO_H
#define __MPO_H
#include "mps.h"

namespace Internal {

template<class Tensor>
class MPO
{
public:
    typedef Tensor TensorT;
    typedef typename Tensor::IndexT IndexT;
    typedef typename Tensor::IndexValT IndexValT;
    typedef typename Tensor::CombinerT CombinerT;
private:
    int N;
    vector<Tensor> A;
    int Lb,Rb;
    const BaseModel* model_;
    Real lref_;
public:
    int minm,maxm;
    Real cutoff;

    operator MPO<IQTensor>()
    { 
        MPO<IQTensor> res(*model_,maxm,cutoff,lref_); 
        res.minm = minm;
        convertToIQ(*model_,A,res.A);
        return res; 
    }

    //Accessor Methods ------------------------------

    int NN() const { return N;}

    // Norm of psi^2 = 1 = norm = sum of denmat evals. 
    // This translates to Tr{Adag A} = norm.  
    // Ref. norm is Tr{1} = d^N, d = 2 S=1/2, d = 4 for Hubbard, etc
    Real lref() const { return lref_; }
    void lref(Real val) 
	{  if(val == 0) { Error("bad lref"); } lref_ = val; }

    const BaseModel& model() const { return *model_; }
    bool is_null() const { return (model_==0); }
    bool is_not_null() const { return (model_!=0); }

    IQIndex si(int i) const { return model_->si(i); }
    IQIndex siP(int i) const { return model_->siP(i); }

    int right_lim() const { return Rb; }
    int left_lim() const { return Lb; }

    const Tensor& AA(int i) const { return GET(A,i); }
    Tensor& AAnc(int i) //nc means 'non const'
    { 
        if(i <= Lb) Lb = i-1;
        if(i >= Rb) Rb = i+1;
        return GET(A,i); 
    }
    Tensor bondTensor(int b) const { Tensor res = A.at(b) * A.at(b+1); return res; }

    //MPO: Constructors -----------------------------------------

    MPO() : N(0), lref_(DefaultLogRef) { }

    MPO(const BaseModel& model, int maxm_ = MAX_M, Real cutoff_ = MAX_CUT, Real _lref = DefaultLogRef) 
    : N(model.NN()), A(N+1), Lb(0), Rb(N), model_(&model), lref_(_lref), minm(1), maxm(maxm_), cutoff(cutoff_)
	{ 
        if(_lref == 0) Error("MPO<Tensor>: Setting lref_ to zero");
        if(_lref == DefaultLogRef) lref_ = model.NN() * log(2.0); 
	}

    MPO(BaseModel& model, istream& s) : N(model.NN()), A(N+1), model_(&model)
    {
        for(int j = 1; j <= N; ++j) A[j].read(s);
        s.read((char*) &Lb,sizeof(Lb));
        s.read((char*) &Rb,sizeof(Rb));
        s.read((char*) &lref_,sizeof(lref_));
        s.read((char*) &minm,sizeof(minm));
        s.read((char*) &maxm,sizeof(maxm));
        s.read((char*) &cutoff,sizeof(cutoff));
    }

    void write(ostream& s) const
    {
        for(int j = 1; j <= N; ++j) A[j].write(s);
        s.write((char*) &Lb,sizeof(Lb));
        s.write((char*) &Rb,sizeof(Rb));
        s.write((char*) &lref_,sizeof(lref_));
        s.write((char*) &minm,sizeof(minm));
        s.write((char*) &maxm,sizeof(maxm));
        s.write((char*) &cutoff,sizeof(cutoff));
    }

    //MPO: index methods --------------------------------------------------

    void mapprime(int oldp, int newp, PrimeType pt = primeBoth)
	{ for(int i = 1; i <= N; ++i) A[i].mapprime(oldp,newp,pt); }

    void primelinks(int oldp, int newp)
	{ for(int i = 1; i <= N; ++i) A[i].mapprime(oldp,newp,primeLink); }

    void noprimelink()
	{ for(int i = 1; i <= N; ++i) A[i].noprime(primeLink); }

    IndexT LinkInd(int i) const { return index_in_common(A[i],A[i+1],Link); }
    IndexT RightLinkInd(int i) const { assert(i<NN()); return index_in_common(AA(i),AA(i+1),Link); }
    IndexT LeftLinkInd(int i)  const { assert(i>1); return index_in_common(AA(i),AA(i-1),Link); }

    void primeall()	// sites i,i' -> i',i'';  link:  l -> l'
	{
        for(int i = 1; i <= this->NN(); i++)
        {
            AAnc(i).mapprime(0,1,primeLink);
            AAnc(i).mapprime(1,2,primeSite);
            AAnc(i).mapprime(0,1,primeSite);
        }
	}

    void doSVD(int i, const Tensor& AA, Direction dir, bool preserve_shape = false)
	{
        tensorSVD(AA,A[i],A[i+1],cutoff,minm,maxm,dir,lref_);
        truncerror = svdtruncerr;

        if(dir == Fromleft)
        {
            if(Lb == i-1 || i == 1) Lb = i;
            if(Rb < i+2) Rb = i+2;
        }
        else
        {
            if(Lb > i-1) Lb = i-1;
            if(Rb == i+2 || i == N-1) Rb = i+1;
        }
	}

    //Move the orthogonality center to site i (left_orth_lim = i-1, right_orth_lim = i+1)
    void position(int i, bool preserve_shape = false)
	{
        if(is_null()) Error("position: MPS is null");
        while(Lb < i-1)
        {
            if(Lb < 0) Lb = 0;
            Tensor WF = AA(Lb+1) * AA(Lb+2);
            doSVD(Lb+1,WF,Fromleft,preserve_shape);
        }
        while(Rb > i+1)
        {
            if(Rb > N+1) Rb = N+1;
            Tensor WF = AA(Rb-2) * AA(Rb-1);
            doSVD(Rb-2,WF,Fromright,preserve_shape);
        }
	}

    bool is_ortho() const { return (Lb + 1 == Rb - 1); }

    int ortho_center() const 
    { 
        if(!is_ortho()) Error("MPS: orthogonality center not well defined.");
        return (Lb + 1);
    }

    bool is_complex() const
    { return A[Lb+1].is_complex(); }

    friend inline ostream& operator<<(ostream& s, const MPO& M)
    {
        s << "\n";
        for(int i = 1; i <= M.NN(); ++i) s << M.AA(i) << "\n";
        return s;
    }

    void print(string name = "",Printdat pdat = HideData) const 
    { printdat = (pdat==ShowData); cerr << "\n" << name << " =\n" << *this << "\n"; printdat = false; }

private:
    friend class MPO<ITensor>;
    friend class MPO<IQTensor>;
}; //class MPO<Tensor>
} //namespace Internal
typedef Internal::MPO<ITensor> MPO;
typedef Internal::MPO<IQTensor> IQMPO;


namespace Internal {

template<class Tensor>
class MPOSet
{
    int N, size_;
    vector<vector<const Tensor*> > A;
public:
    typedef vector<Tensor> TensorT;

    MPOSet() : N(-1), size_(0) { }

    MPOSet(const MPS<Tensor>& Op1) 
    : N(-1), size_(0) 
    { include(Op1); }

    MPOSet(const MPS<Tensor>& Op1, 
           const MPS<Tensor>& Op2) 
    : N(-1), size_(0) 
    { include(Op1); include(Op2); }

    MPOSet(const MPS<Tensor>& Op1, 
           const MPS<Tensor>& Op2,
           const MPS<Tensor>& Op3) 
    : N(-1), size_(0) 
    { include(Op1); include(Op2); include(Op3); }

    MPOSet(const MPS<Tensor>& Op1, 
           const MPS<Tensor>& Op2,
           const MPS<Tensor>& Op3, 
           const MPS<Tensor>& Op4) 
    : N(-1), size_(0) 
    { include(Op1); include(Op2); include(Op3); include(Op4); }

    void include(const MPS<Tensor>& Op)
    {
        if(N < 0) { N = Op.NN(); A.resize(N+1); }
        for(int n = 1; n <= N; ++n) GET(A,n).push_back(&(Op.AA(n))); 
        ++size_;
    }

    int NN() const { return N; }
    int size() const { return size_; }
    const vector<const Tensor*>& AA(int j) const { return GET(A,j); }
    const vector<Tensor> bondTensor(int b) const
    { vector<Tensor> res = A[b] * A[b+1]; return res; }

}; //class Internal::MPOSet

} //namespace Internal
typedef Internal::MPOSet<ITensor> MPOSet;
typedef Internal::MPOSet<IQTensor> IQMPOSet;

template <class MPSType, class MPOType>
void psiHphi(const MPSType& psi, const MPOType& H, const MPSType& phi, Real& re, Real& im) //<psi|H|phi>
{
    typedef typename MPSType::TensorT Tensor;
    const int N = H.NN();

    if(psi.NN() != phi.NN() || psi.NN() != N) Error("psiHphi: mismatched N");

    Tensor L = phi.AA(1) * H.AA(1) * conj(primed(psi.AA(1)));
    for(int i = 2; i < N; ++i) { L = L * phi.AA(i) * H.AA(i) * conj(primed(psi.AA(i))); }
    L = L * phi.AA(N) * H.AA(N);

    Dot(primed(psi.AA(N)),L,re,im);
}
template <class MPSType, class MPOType>
Real psiHphi(const MPSType& psi, const MPOType& H, const MPSType& phi) //Re[<psi|H|phi>]
{
    Real re, im;
    psiHphi(psi,H,phi,re,im);
    if(im != 0) cerr << format("\nReal psiHphi: WARNING, dropping non-zero (im = %.5f) imaginary part of expectation value.\n")%im;
    return re;
}

inline void psiHphi(const MPS& psi, const MPO& H, const ITensor& LB, const ITensor& RB, const MPS& phi, Real& re, Real& im) //<psi|H|phi>
{
    int N = psi.NN();
    if(N != phi.NN() || H.NN() < N) Error("mismatched N in psiHphi");
    MPS psiconj(psi);
    for(int i = 1; i <= N; ++i) psiconj.AAnc(i) = conj(primed(psi.AA(i)));
    ITensor L = (LB.is_null() ? phi.AA(1) : LB * phi.AA(1));
    L *= H.AA(1); L *= psiconj.AA(1);
    for(int i = 2; i <= N; ++i)
	{ L *= phi.AA(i); L *= H.AA(i); L *= psiconj.AA(i); }
    if(!RB.is_null()) L *= RB;
    if(L.is_complex())
    {
        if(L.Length() != 1) Error("Non-scalar result in psiHphi.");
        const int sign = (L.neg() ? -1 : 1);
        re = sign * L(IndReIm(1)) * exp(L.logfac());
        im = sign * L(IndReIm(2)) * exp(L.logfac());
    }
    else 
    {
        if(L.Length() != 1) Error("Non-scalar result in psiHphi.");
        re = L.val0()*(L.neg() ? -1 : 1)*exp(L.logfac());
        im = 0;
    }
}
inline Real psiHphi(const MPS& psi, const MPO& H, const ITensor& LB, const ITensor& RB, const MPS& phi) //Re[<psi|H|phi>]
{
    Real re,im; psiHphi(psi,H,LB,RB,phi,re,im);
    if(im != 0) cerr << "Real psiHphi: WARNING, dropping non-zero imaginary part of expectation value.\n";
    return re;
}

#endif
