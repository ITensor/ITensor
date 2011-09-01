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

    IndexT LinkInd(int i) const { IndexT res; index_in_common(A[i],A[i+1],Link,res); return res; }
    IndexT RightLinkInd(int i) const { assert(i<NN()); IndexT res; index_in_common(AA(i),AA(i+1),Link,res); return res; }
    IndexT LeftLinkInd(int i)  const { assert(i>1); IndexT res; index_in_common(AA(i),AA(i-1),Link,res); return res; }

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
    if(phi.NN() != N || psi.NN() != N) Error("psiHphi: mismatched N");

    Tensor L = phi.AA(1); L *= H.AA(1); L *= conj(primed(psi.AA(1)));
    for(int i = 2; i < N; ++i) 
    { L *= phi.AA(i); L *= H.AA(i); L *= conj(primed(psi.AA(i))); }
    L *= phi.AA(N); L *= H.AA(N);

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

inline void psiHKphi(const IQMPS& psi, const IQMPO& H, const IQMPO& K,const IQMPS& phi, Real& re, Real& im) //<psi|H K|phi>
{
    if(psi.NN() != phi.NN() || psi.NN() != H.NN() || psi.NN() != K.NN()) Error("Mismatched N in psiHKphi");
    int N = psi.NN();
    IQMPS psiconj(psi);
    for(int i = 1; i <= N; i++)
	{
        psiconj.AAnc(i) = conj(psi.AA(i));
        psiconj.AAnc(i).mapprime(0,2);
	}
    IQMPO Kp(K);
    Kp.mapprime(1,2);
    Kp.mapprime(0,1);

    //scales as m^2 k^2 d
    IQTensor L = (((phi.AA(1) * H.AA(1)) * Kp.AA(1)) * psiconj.AA(1));
    for(int i = 2; i < N; i++)
    {
        //scales as m^3 k^2 d + m^2 k^3 d^2
        L = ((((L * phi.AA(i)) * H.AA(i)) * Kp.AA(i)) * psiconj.AA(i));
    }
    //scales as m^2 k^2 d
    L = ((((L * phi.AA(N)) * H.AA(N)) * Kp.AA(N)) * psiconj.AA(N)) * IQTSing;
    //cout << "in psiHKpsi, L is "; PrintDat(L);
    L.GetSingComplex(re,im);
}
inline Real psiHKphi(const IQMPS& psi, const IQMPO& H, const IQMPO& K,const IQMPS& phi) //<psi|H K|phi>
{
    Real re,im;
    psiHKphi(psi,H,K,phi,re,im);
    if(fabs(im) > 1E-12) Error("Non-zero imaginary part in psiHKphi");
    return re;
}

inline void nmultMPO(const IQMPO& Aorig, const IQMPO& Borig, IQMPO& res,Real cut, int maxm)
{
    if(Aorig.NN() != Borig.NN()) Error("nmultMPO(IQMPO): Mismatched N");
    int N = Borig.NN();
    IQMPO A(Aorig), B(Borig);

    A.position(1);
    B.position(1);
    B.primeall();

    res=A;
    res.primelinks(0,4);
    res.mapprime(1,2,primeSite);

    IQTensor clust,nfork;
    vector<int> midsize(N);
    for(int i = 1; i < N; ++i)
	{
        if(i == 1) { clust = A.AA(i) * B.AA(i); }
        else       { clust = nfork * A.AA(i) * B.AA(i); }
        if(i == N-1) break;

        IQIndex oldmid = res.RightLinkInd(i);
        nfork = IQTensor(A.RightLinkInd(i),B.RightLinkInd(i),oldmid);
        if(clust.iten_size() == 0)	// this product gives 0 !!
        { cerr << format("WARNING: clust.iten_size()==0 in nmultMPO (i=%d).\n")%i; res = IQMPO(); return; }
        tensorSVD(clust, res.AAnc(i), nfork,cut,1,maxm,Fromleft,A.lref());
        IQIndex mid; index_in_common(res.AA(i),nfork,Link,mid);
        assert(mid.dir() == In);
        mid.conj();
        midsize[i] = mid.m();
        assert(res.RightLinkInd(i+1).dir() == Out);
        assert(res.si(i+1).dir() == Out);
        res.AAnc(i+1) = IQTensor(mid,conj(res.si(i+1)),res.si(i+1).primed().primed(),res.RightLinkInd(i+1));
	}

    nfork = clust * A.AA(N) * B.AA(N);
    if(nfork.iten_size() == 0)	// this product gives 0 !!
    { cerr << "WARNING: nfork.iten_size()==0 in nmultMPO\n"; res = IQMPO(); return; }

    res.doSVD(N-1,nfork,Fromright,false);
    res.noprimelink();
    res.mapprime(2,1,primeSite);
    res.cutoff = cut;
    res.position(N);
    res.position(1);

}//void nmultMPO(const IQMPO& Aorig, const IQMPO& Borig, IQMPO& res,Real cut, int maxm)

inline void napplyMPO(const IQMPS& x, const IQMPO& K, IQMPS& res, Real cutoff, int maxm)
{
    if(cutoff < 0) cutoff = x.cutoff;
    if(maxm < 0) maxm = x.maxm;
    int N = x.NN();
    if(K.NN() != N) Error("Mismatched N in napplyMPO");
    if(x.right_lim() > 3)
    {
        cerr << "x is " << endl << x << endl;
        Error("bad right_lim for x");
    }
    if(K.right_lim() > 3)
    {
        //cerr << "K is " << endl << K << endl;
        Error("bad right_lim for K");
    }

    res = x; res.maxm = maxm; res.cutoff = cutoff;
    res.primelinks(0,4);
    res.mapprime(0,1,primeSite);

    IQTensor clust,nfork;
    vector<int> midsize(N);
    int maxdim = 1;
    for(int i = 1; i < N; i++)
	{
        if(i == 1) { clust = x.AA(i) * K.AA(i); }
        else { clust = nfork * (x.AA(i) * K.AA(i)); }
        if(i == N-1) break; //No need to SVD for i == N-1

        IQIndex oldmid = res.RightLinkInd(i); assert(oldmid.dir() == Out);
        nfork = IQTensor(x.RightLinkInd(i),K.RightLinkInd(i),oldmid);
        if(clust.iten_size() == 0)	// this product gives 0 !!
        { res = IQMPS(); return; }
        tensorSVD(clust, res.AAnc(i), nfork,cutoff,1,maxm,Fromleft);
        IQIndex mid; index_in_common(res.AA(i),nfork,Link,mid);
        assert(mid.dir() == In);
        mid.conj();
        midsize[i] = mid.m();
        maxdim = max(midsize[i],maxdim);
        assert(res.RightLinkInd(i+1).dir() == Out);
        res.AAnc(i+1) = IQTensor(mid,res.si(i+1).primed(),res.RightLinkInd(i+1));
	}
    nfork = clust * x.AA(N) * K.AA(N);
    if(nfork.iten_size() == 0)	// this product gives 0 !!
	{ res = IQMPS(); return; }

    res.doSVD(N-1,nfork,Fromright,false);
    res.noprimelink();
    res.mapprime(1,0,primeSite);
    res.position(1);
    res.maxm = x.maxm; res.cutoff = x.cutoff;

} //void napplyMPO

//Expensive: scales as m^3 k^3!
inline void exact_applyMPO(const IQMPS& x, const IQMPO& K, IQMPS& res)
{
    int N = x.NN();
    if(K.NN() != N) Error("Mismatched N in exact_applyMPO");

    res = x;
    res.position(1);

    res.AAnc(1) = x.AA(1) * K.AA(1);
    for(int j = 1; j < N; ++j)
	{
        //cerr << format("exact_applyMPO: step %d\n") % j;
        //Compute product of MPS tensor and MPO tensor
        res.AAnc(j+1) = x.AA(j+1) * K.AA(j+1); //m^2 k^2 d^2

        //Add common IQIndices to IQCombiner
        IQCombiner comb; comb.doCondense(false);
        foreach(const IQIndex& I, res.AA(j).iqinds())
        if(res.AA(j+1).hasindex(I) && I != IQIndReIm && I.type() != Virtual)
        { assert(I.dir() == Out); comb.addleft(I);}
        comb.init(nameint("a",j));

        //Apply combiner to product tensors
        res.AAnc(j) = res.AA(j) * comb; //m^3 k^3 d
        res.AAnc(j+1) = conj(comb) * res.AA(j+1); //m^3 k^3 d
	}
    res.mapprime(1,0,primeSite);
    //res.position(1);
} //void exact_applyMPO

#endif
