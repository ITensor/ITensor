#ifndef __MPO_H
#define __MPO_H
#include "mps.h"


template<class Tensor>
class MPOt : private MPSt<Tensor>
{
public:
    typedef Tensor TensorT;
    typedef typename Tensor::IndexT IndexT;
    typedef typename Tensor::IndexValT IndexValT;
    typedef typename Tensor::CombinerT CombinerT;
private:
    typedef MPSt<Tensor> Parent;
    using Parent::N;
    using Parent::A;
    using Parent::left_orth_lim;
    using Parent::right_orth_lim;
    using Parent::model_;
    using Parent::svd_;
public:

    operator MPOt<IQTensor>()
    { 
        MPOt<IQTensor> res(*model_,maxm(),cutoff(),doRelCutoff(),refNorm()); 
        res.svd_ = svd_;
        convertToIQ(*model_,A,res.A);
        return res; 
    }

    //Accessor Methods ------------------------------

    using Parent::NN;

    using Parent::model;
    using Parent::is_null;
    using Parent::is_not_null;

    using Parent::si;
    using Parent::siP;

    using Parent::right_lim;
    using Parent::left_lim;

    using Parent::AA;
    using Parent::AAnc;
    using Parent::bondTensor;

    using Parent::doRelCutoff;
    using Parent::refNorm;
    using Parent::cutoff;
    using Parent::minm;
    using Parent::maxm;
    using Parent::truncerr;
    using Parent::eigs_kept;
    using Parent::showeigs;
    using Parent::svd;

    //MPOt: Constructors -----------------------------------------

    MPOt() : Parent() { doRelCutoff(true); }

    MPOt(const BaseModel& model, int maxm_ = MAX_M, Real cutoff_ = MAX_CUT, 
    bool _doRelCutoff = true, LogNumber _refNorm = DefaultRefScale) 
    : Parent(model,maxm_,cutoff_)
	{ 
        doRelCutoff(_doRelCutoff);
        refNorm(_refNorm);
        // Norm of psi^2 = 1 = norm = sum of denmat evals. 
        // This translates to Tr{Adag A} = norm.  
        // Ref. norm is Tr{1} = d^N, d = 2 S=1/2, d = 4 for Hubbard, etc
        if(_refNorm == DefaultRefScale) refNorm(exp(model.NN()));
	}

    MPOt(BaseModel& model, istream& s) { read(model,s); }

    virtual ~MPOt() { }

    using Parent::read;
    using Parent::write;

    //MPOt: operators ------------------------------------------------------

    MPOt& operator*=(Real a) { Parent::operator*=(a); return *this; }
    inline MPOt operator*(Real r) const { MPOt res(*this); res *= r; return res; }
    friend inline MPOt operator*(Real r, MPOt res) { res *= r; return res; }

    MPOt& operator+=(const MPOt& oth) { Parent::operator+=(oth); return *this; }
    inline MPOt operator+(MPOt res) const { res += *this; return res; }
    inline MPOt operator-(MPOt res) const { res *= -1; res += *this; return res; }

    //MPOt: index methods --------------------------------------------------

    using Parent::mapprime;
    using Parent::primelinks;
    using Parent::noprimelink;

    using Parent::LinkInd;
    using Parent::RightLinkInd;
    using Parent::LeftLinkInd;

    void primeall()	// sites i,i' -> i',i'';  link:  l -> l'
	{
        for(int i = 1; i <= this->NN(); i++)
        {
            AAnc(i).mapprime(0,1,primeLink);
            AAnc(i).mapprime(1,2,primeSite);
            AAnc(i).mapprime(0,1,primeSite);
        }
	}

    using Parent::doSVD;
    using Parent::position;

    using Parent::is_ortho;
    using Parent::ortho_center;
    using Parent::orthogonalize;

    using Parent::is_complex;

    using Parent::applygate;

    friend inline ostream& operator<<(ostream& s, const MPOt& M)
    {
        s << "\n";
        for(int i = 1; i <= M.NN(); ++i) s << M.AA(i) << "\n";
        return s;
    }

    using Parent::print;

private:
    friend class MPOt<ITensor>;
    friend class MPOt<IQTensor>;
}; //class MPOt<Tensor>
typedef MPOt<ITensor> MPO;
typedef MPOt<IQTensor> IQMPO;


namespace Internal {

template<class Tensor>
class MPOSet
{
    int N, size_;
    vector<vector<const Tensor*> > A;
public:
    typedef vector<Tensor> TensorT;

    MPOSet() : N(-1), size_(0) { }

    MPOSet(const MPOt<Tensor>& Op1) 
    : N(-1), size_(0) 
    { include(Op1); }

    MPOSet(const MPOt<Tensor>& Op1, 
           const MPOt<Tensor>& Op2) 
    : N(-1), size_(0) 
    { include(Op1); include(Op2); }

    MPOSet(const MPOt<Tensor>& Op1, 
           const MPOt<Tensor>& Op2,
           const MPOt<Tensor>& Op3) 
    : N(-1), size_(0) 
    { include(Op1); include(Op2); include(Op3); }

    MPOSet(const MPOt<Tensor>& Op1, 
           const MPOt<Tensor>& Op2,
           const MPOt<Tensor>& Op3, 
           const MPOt<Tensor>& Op4) 
    : N(-1), size_(0) 
    { include(Op1); include(Op2); include(Op3); include(Op4); }

    void include(const MPOt<Tensor>& Op)
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
    if(im != 0) cerr << boost::format("\nReal psiHphi: WARNING, dropping non-zero (im = %.5f) imaginary part of expectation value.\n")%im;
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
        if(L.vec_size() != 2) Error("Non-scalar result in psiHphi.");
        re = L(IndReIm(1));
        im = L(IndReIm(2));
    }
    else 
    {
        if(L.vec_size() != 1) Error("Non-scalar result in psiHphi.");
        re = L.val0();
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

template <class MPOType>
void nmultMPO(const MPOType& Aorig, const MPOType& Borig, MPOType& res,Real cut, int maxm);
void napplyMPO(const IQMPS& x, const IQMPO& K, IQMPS& res, Real cutoff, int maxm);
void exact_applyMPO(const IQMPS& x, const IQMPO& K, IQMPS& res);

#endif
