#ifndef __MPS_H
#define __MPS_H
#include "tensor.h"
#include "iq.h"

Vector do_denmat_Real(const ITensor& nA, ITensor& A, ITensor& B, Real cutoff,int minm, int maxm, Direction dir);
Vector do_denmat_Real(const IQTensor& nA, IQTensor& A, IQTensor& B, Real cutoff, int minm,int maxm, Direction dir);
Vector do_denmat_Real(const vector<IQTensor>& nA, const IQTensor& A, const IQTensor& B, IQTensor& U,
	Real cutoff,Real lognormfac,int minm,int maxm, Direction dir, bool donormalize, bool do_relative_cutoff);
void getCenterMatrix(ITensor& A, const Index& bond, Real cutoff,int minm, int maxm, ITensor& Lambda, string newbondname = "");
void plussers(const Index& l1, const Index& l2, Index& sumind, ITensor& first, ITensor& second);
void plussers(const IQIndex& l1, const IQIndex& l2, IQIndex& sumind, IQTensor& first, IQTensor& second);

template<class Tensor, class TensorSet>
Real doDavidson(Tensor& phi, const TensorSet& mpoh, const TensorSet& LH, const TensorSet& RH, int niter, int debuglevel, Real errgoal);

//class MPS; void psiphi(const MPS& psi, const MPS& phi, Real& re, Real& im); //Allows class MPS to call this later

extern Real truncerror, svdtruncerr;

namespace Internal {
template <class IndexT>
class SiteSet
{
protected:
    int N;
    vector<IndexT> site;
    virtual void initSites(vector<IndexT>& site) = 0;
public:
    virtual int NN() const { return N; }
    virtual IndexT si(int i) const { return site.at(i); }
    virtual IndexT siP(int i) const { return site.at(i).primed(); }
    virtual ~SiteSet() { }
    SiteSet() : N(-1) { }
    SiteSet(int nsite) : N(nsite), site(nsite+1) { }
    SiteSet(istream& s) { read(s); }

    virtual void read(istream& s)
    {
        s.read((char*) &N,sizeof(N));
        site.resize(N+1);
        for(int j = 1; j <= N; ++j) site.at(j) = IndexT(s);
    }
    virtual void write(ostream& s) const
    {
        s.write((char*) &N,sizeof(N));
        for(int j = 1; j <= N; ++j) site.at(j).write(s);
    }
};
} //namespace Internal
typedef Internal::SiteSet<IQIndex> SiteSet;

template <class IndexT>
ostream& operator<<(ostream& s, const Internal::SiteSet<IndexT>& st)
{
    s << "SiteSet:\n";
    for(int j = 1; j <= st.NN(); ++j) s << format("si(%d) = ")%j << st.si(j) << "\n";
    return s;
}

//---------------------------------------------------------
//Definition of Model Types
//---------------------------------------------------------

namespace SpinOne {

const int Dim = 3;

class Model : public SiteSet
{
    typedef SiteSet Parent;

    //void initSites(vector<Index>& site)
    //{
    //    for(int i = 1; i <= this->N; ++i) site.at(i) = Index(nameint("S=1, site=",i),Dim,Site);
    //}

    void initSites(vector<IQIndex>& site)
    {
        for(int i = 1; i <= this->N; ++i)
        {
            site.at(i) = IQIndex(nameint("S=1, site=",i),
            Index(nameint("Up for site",i),1,Site),QN(+2,0),
            Index(nameint("Z0 for site",i),1,Site),QN( 0,0),
            Index(nameint("Dn for site",i),1,Site),QN(-2,0));
        }
    }

public:
    Model() : Parent() { }
    Model(int nsite) : Parent(nsite) { initSites(this->site); }
    Model(istream& s) : Parent(s) { }
    virtual ~Model() { }

    IQIndex si(int i) const { return Parent::si(i); }
    IQIndex siP(int i) const { return Parent::siP(i); }

    IQIndexVal Up(int i) const { return si(i)(1); }
    IQIndexVal Z0(int i) const { return si(i)(2); }
    IQIndexVal Dn(int i) const { return si(i)(3); }

    IQIndexVal UpP(int i) const { return siP(i)(1); }
    IQIndexVal Z0P(int i) const { return siP(i)(2); }
    IQIndexVal DnP(int i) const { return siP(i)(3); }

    SiteOp id(int i) const
    {
        SiteOp res(si(i));
        res(Up(i),UpP(i)) = 1; res(Z0(i),Z0P(i)) = 1; res(Dn(i),DnP(i)) = 1;
        return res;
    }

    SiteOp sz(int i) const
    {
        SiteOp res(si(i));
        res(Up(i),UpP(i)) = 1; res(Dn(i),DnP(i)) = -1;
        return res;
    }

    SiteOp sx(int i) const
    {
        SiteOp res(si(i));
        res(Up(i),Z0P(i)) = ISqrt2; res(Z0(i),UpP(i)) = ISqrt2;
        res(Z0(i),DnP(i)) = ISqrt2; res(Dn(i),Z0P(i)) = ISqrt2;
        return res;
    }

    SiteOp isy(int i) const
    {
        SiteOp res(si(i));
        res(Up(i),Z0P(i)) = +ISqrt2; res(Z0(i),UpP(i)) = -ISqrt2;
        res(Z0(i),DnP(i)) = +ISqrt2; res(Dn(i),Z0P(i)) = -ISqrt2;
        return res;
    }

    SiteOp sp(int i) const
    {
        SiteOp res(si(i));
        res(Dn(i),Z0P(i)) = Sqrt2; res(Z0(i),UpP(i)) = Sqrt2;
        return res;
    }

    SiteOp sm(int i) const
    {
        SiteOp res(si(i));
        res(Up(i),Z0P(i)) = Sqrt2; res(Z0(i),DnP(i)) = Sqrt2;
        return res;
    }

    SiteOp sz2(int i) const
    {
        SiteOp res(si(i));
        res(Up(i),UpP(i)) = 1; res(Dn(i),DnP(i)) = 1;
        return res;
    }

    SiteOp sx2(int i) const
    {
        SiteOp res(si(i));
        res(Up(i),UpP(i)) = 0.5; res(Up(i),DnP(i)) = 0.5;
        res(Z0(i),Z0P(i)) = 1;
        res(Dn(i),DnP(i)) = 0.5; res(Dn(i),UpP(i)) = 0.5;
        return res;
    }

    SiteOp sy2(int i) const
    {
        SiteOp res(si(i));
        res(Up(i),UpP(i)) = 0.5; res(Up(i),DnP(i)) = -0.5;
        res(Z0(i),Z0P(i)) = 1;
        res(Dn(i),DnP(i)) = 0.5; res(Dn(i),UpP(i)) = -0.5;
        return res;
    }

};

} //end namespace SpinOne

namespace SpinHalf {

const int Dim = 2;

class Model : public SiteSet
{
    typedef SiteSet Parent;

    //void initSites(vector<Index>& site)
    //{ for(int i = 1; i <= this->N; ++i) site.at(i) = Index(nameint("S=1/2, site=",i),Dim,Site); }

    void initSites(vector<IQIndex>& site)
    {
        for(int i = 1; i <= this->N; ++i)
        {
            site.at(i) = IQIndex(nameint("S=1/2, site=",i),
            Index(nameint("Up for site",i),1,Site),QN(+1,0),
            Index(nameint("Dn for site",i),1,Site),QN(-1,0));
        }
    }
public:
    Model() : Parent() { }
    Model(int nsite) : Parent(nsite) { initSites(this->site); }
    Model(istream& s) : Parent(s) { }
    virtual ~Model() { }

    IQIndex si(int i) const { return Parent::si(i); }
    IQIndex siP(int i) const { return Parent::siP(i); }

    IQIndexVal Up(int i) const { return si(i)(1); }
    IQIndexVal Dn(int i) const { return si(i)(2); }

    IQIndexVal UpP(int i) const { return siP(i)(1); }
    IQIndexVal DnP(int i) const { return siP(i)(2); }

    SiteOp id(int i) const
    {
        SiteOp res(si(i));
        res(Up(i),UpP(i)) = 1; res(Dn(i),DnP(i)) = 1;
        return res;
    }

    SiteOp sz(int i) const
    {
        SiteOp res(si(i));
        res(Up(i),UpP(i)) = 0.5; res(Dn(i),DnP(i)) = -0.5;
        return res;
    }

    SiteOp sx(int i) const
    {
        SiteOp res(si(i));
        res(Up(i),DnP(i)) = 0.5; res(Dn(i),UpP(i)) = 0.5;
        return res;
    }

    SiteOp isy(int i) const
    {
        SiteOp res(si(i));
        res(Up(i),DnP(i)) = 0.5; res(Dn(i),UpP(i)) = -0.5;
        return res;
    }

    //S^+
    SiteOp sp(int i) const
    {
        SiteOp res(si(i));
        res(Dn(i),UpP(i)) = 1;
        return res;
    }

    //S^-
    SiteOp sm(int i) const
    {
        SiteOp res(si(i));
        res(Up(i),DnP(i)) = 1;
        return res;
    }

};

} //end namespace SpinHalf

namespace Fermion {

const int Dim = 2;

class Model : public SiteSet
{
public:
    bool odd_even_up_down,conserve_Nf;
private:
    typedef SiteSet Parent;

    //void initSites(vector<Index>& site)
    //{
    //    if(odd_even_up_down)
    //    {
    //    for(int i = 1; i <= this->N; ++i) 
    //    {
    //        if(i%2==1) site.at(i) = Index(nameint("Fermion, Up site=",i),Dim,Site);
    //        else site.at(i) = Index(nameint("Fermion, Dn site=",i),Dim,Site);
    //    }
    //    }
    //    else
    //    { for(int i = 1; i <= this->N; ++i) site.at(i) = Index(nameint("Fermion, site=",i),Dim,Site); }
    //}

    void initSites(vector<IQIndex>& site)
    {
        int occ = (conserve_Nf ? 1 : 0);
        if(odd_even_up_down)
        {
        for(int i = 1; i <= this->N; ++i)
        {
            if(i%2==1)
            {
                //site.at(i) = IQIndex(nameint("Fermion, Up site=",i),Site);
                //site.at(i).insert(Index(nameint("Emp for Up site",i),1,Site),QN(0,0,0));
                //site.at(i).insert(Index(nameint("Occ for Up site",i),1,Site),QN(+1,occ,1));
                site.at(i) = IQIndex(nameint("Fermion, Up site=",i),
                Index(nameint("Emp for Up site",i),1,Site),QN(0,0,0),
                Index(nameint("Occ for Up site",i),1,Site),QN(+1,occ,1));
            }
            else
            {
                site.at(i) = IQIndex(nameint("Fermion, Dn site=",i),
                Index(nameint("Emp for Dn site",i),1,Site),QN(0,0,0),
                Index(nameint("Occ for Dn site",i),1,Site),QN(-1,occ,1));
            }
        }
        }
        else
        {
        for(int i = 1; i <= this->N; ++i)
        {
            site.at(i) = IQIndex(nameint("Fermion, site=",i),
            Index(nameint("Emp for site",i),1,Site),QN(0,0,0),
            Index(nameint("Occ for site",i),1,Site),QN(0,occ,1));
        }
        }
    }
public:
    Model() : odd_even_up_down(false),conserve_Nf(true) { }
    Model(int nsite, bool odd_even_up_down_ = false, bool conserve_Nf_ = true) 
    : Parent(nsite), odd_even_up_down(odd_even_up_down_), conserve_Nf(conserve_Nf_) { initSites(this->site); }
    virtual ~Model() { }
    Model(istream& s) : Parent(s) 
    { 
        s.read((char*) &odd_even_up_down,sizeof(odd_even_up_down));
        s.read((char*) &conserve_Nf,sizeof(conserve_Nf));
    }
    void write(ostream& s) const
    {
        Parent::write(s);
        s.write((char*) &odd_even_up_down,sizeof(odd_even_up_down));
        s.write((char*) &conserve_Nf,sizeof(conserve_Nf));
    }

    IQIndex si(int i) const { return Parent::si(i); }
    IQIndex siP(int i) const { return Parent::siP(i); }

    IQIndexVal Emp(int i) const { return si(i)(1); }
    IQIndexVal Occ(int i) const { return si(i)(2); }

    IQIndexVal EmpP(int i) const { return siP(i)(1); }
    IQIndexVal OccP(int i) const { return siP(i)(2); }

    SiteOp id(int i) const
    {
        SiteOp res(si(i));
        res(Emp(i),EmpP(i)) = 1; res(Occ(i),OccP(i)) = 1;
        return res;
    }

    SiteOp C(int i) const
    {
        SiteOp res(si(i));
        res(Occ(i),EmpP(i)) = 1;
        return res;
    }

    SiteOp Cdag(int i) const
    {
        SiteOp res(si(i));
        res(Emp(i),OccP(i)) = 1;
        return res;
    }

    SiteOp n(int i) const
    {
        SiteOp res(si(i));
        res(Occ(i),OccP(i)) = 1;
        return res;
    }

    //String operator F_i = (-1)^{n_i} = (1-2*n_i)
    SiteOp FermiPhase(int i) const
    {
        SiteOp res(si(i));
        res(Emp(i),EmpP(i)) = 1; res(Occ(i),OccP(i)) = -1;
        return res;
    }
};

} //end namespace Spinless

namespace Hubbard {	// Full Hubbard sites, srw 8/10/11

const int Dim = 4;

class Model : public SiteSet
{
private:
    typedef SiteSet Parent;

    void initSites(vector<IQIndex>& site)
	{
        for(int i = 1; i <= this->N; ++i)
	    {
	    site.at(i) = IQIndex(nameint("Hubbard, site=",i),
		    Index(nameint("Emp for site ",i),1,Site),  QN( 0,0,0),
		    Index(nameint("Up for site ",i),1,Site),   QN(+1,1,1),
		    Index(nameint("Dn for site ",i),1,Site),   QN(-1,1,1),
		    Index(nameint("Up-Dn for site ",i),1,Site),QN( 0,2,0));
	    }
	}
public:
    Model() { }
    Model(int nsite) 
    : Parent(nsite)  { initSites(this->site); }
    virtual ~Model() { }

    IQIndex si(int i) const { return Parent::si(i); }
    IQIndex siP(int i) const { return Parent::siP(i); }

    IQIndexVal Emp(int i) const { return si(i)(1); }
    IQIndexVal UpState(int i) const { return si(i)(2); }
    IQIndexVal DnState(int i) const { return si(i)(3); }
    IQIndexVal UpDnState(int i) const { return si(i)(4); } // cdag_dn cdag_up | vac >

    IQIndexVal EmpP(int i) const { return siP(i)(1); }
    IQIndexVal UpStateP(int i) const { return siP(i)(2); }
    IQIndexVal DnStateP(int i) const { return siP(i)(3); }
    IQIndexVal UpDnStateP(int i) const { return siP(i)(4); }

    SiteOp id(int i) const
	{
	SiteOp res(si(i));
	res(Emp(i),EmpP(i)) = 1; res(UpState(i),UpStateP(i)) = 1;
	res(DnState(i),DnStateP(i)) = 1; res(UpDnState(i),UpDnStateP(i)) = 1;
	return res;
	}

    SiteOp Cup(int i) const
	{
	SiteOp res(si(i));
	res(UpState(i),EmpP(i)) = 1;
	res(UpDnState(i),DnStateP(i)) = -1;
	return res;
	}

    SiteOp Cdagup(int i) const
	{
	SiteOp res(si(i));
	res(UpStateP(i),Emp(i)) = 1;
	res(UpDnStateP(i),DnState(i)) = -1;
	return res;
	}

    SiteOp Cdn(int i) const
	{
	SiteOp res(si(i));
	res(DnState(i),EmpP(i)) = 1;
	res(UpDnState(i),UpStateP(i)) = 1;
	return res;
	}

    SiteOp Cdagdn(int i) const
	{
	SiteOp res(si(i));
	res(DnStateP(i),Emp(i)) = 1;
	res(UpDnStateP(i),UpState(i)) = 1;
	return res;
	}

    SiteOp Nup(int i) const
	{
	SiteOp res(si(i));
	res(UpState(i),UpStateP(i)) = 1;
	res(UpDnState(i),UpDnStateP(i)) = 1;
	return res;
	}

    SiteOp Ndn(int i) const
	{
	SiteOp res(si(i));
	res(DnState(i),DnStateP(i)) = 1;
	res(UpDnState(i),UpDnStateP(i)) = 1;
	return res;
	}

    SiteOp Ntot(int i) const
	{
	SiteOp res(si(i));
	res(Emp(i),EmpP(i)) = 0; res(UpState(i),UpStateP(i)) = 1;
	res(DnState(i),DnStateP(i)) = 1; res(UpDnState(i),UpDnStateP(i)) = 2;
	return res;
	}

    SiteOp NupNdn(int i) const
	{
	SiteOp res(si(i));
	res(UpDnState(i),UpDnStateP(i)) = 1;
	return res;
	}

    SiteOp Sz(int i) const
	{
	SiteOp res(si(i));
	res(UpState(i),UpStateP(i)) = 0.5;
	res(DnState(i),DnStateP(i)) = -0.5;
	return res;
	}

    //String operator F_i = (-1)^{n_i} = (1-2*n_i)
    SiteOp FermiPhase(int i) const
	{
	SiteOp res(si(i));
	res(Emp(i),EmpP(i)) = 1; res(UpState(i),UpStateP(i)) = -1;
	res(DnState(i),DnStateP(i)) = -1; res(UpDnState(i),UpDnStateP(i)) = 1;
	return res;
	}
};

} //end namespace Hubbard

//---------------------------------------------------------
//---------------------------------------------------------
//---------------------------------------------------------

class InitState
{
    vector<IQIndexVal> state;
public:
    InitState(int nsite) : state(nsite+1) { }
    IQIndexVal& operator()(int i) { return GET(state,i-1); }
    const IQIndexVal& operator()(int i) const { return GET(state,i-1); }
    operator vector<IQIndexVal>() const { return state; }
};

namespace Internal {

template <class Tensor>
class MPS
{
public:
    typedef Tensor TensorT;
    typedef typename Tensor::IndexT IndexT;
    typedef typename Tensor::IndexValT IndexValT;
    typedef Internal::SiteSet<IQIndex> SiteSetT;
protected:
    int N;
    vector<Tensor> A;
    int left_orth_lim,right_orth_lim;


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
        A_[1] = IQTensor(si(1),a[1],Center); A_[1](initState(1))=1;
        for(int i = 2; i < N; ++i)
        { 
        A_[i] = IQTensor(conj(a[i-1]),si(i),a[i]); 
        A_[i](initState(i))=1;
        }
        A_[N] = IQTensor(conj(a[N-1]),si(N)); A_[N](initState(N))=1;
    }

    typedef pair<typename vector<Tensor>::const_iterator,typename vector<Tensor>::const_iterator> const_range_type;
public:
    const SiteSetT* sst_;
    int minm,maxm;
    Real cutoff;

    //Accessor Methods ------------------------------

    int NN() const { return N;}
    IQIndex si(int i) const { return sst_->si(i); }
    IQIndex siP(int i) const { return sst_->siP(i); }
    typedef typename vector<Tensor>::const_iterator AA_it;
    const pair<AA_it,AA_it> AA() const { return make_pair(A.begin()+1,A.end()); }
    const Tensor& AA(int i) const { return GET(A,i); }
    const SiteSetT& sst() const { return *sst_; }
    Tensor& AAnc(int i) //nc means 'non const'
    { 
        if(i <= left_orth_lim) left_orth_lim = i-1;
        if(i >= right_orth_lim) right_orth_lim = i+1;
        return GET(A,i); 
    }
    bool is_null() const { return A[1].is_null(); }
    bool is_not_null() const { return A[1].is_not_null(); }

    //Useful for iterating over A's in a foreach loop
    //Do foreach(const I[Q]Tensor& A, psi.all_As()) { ... }
    const_range_type all_As() const { return make_pair(A.begin(),A.end()-1); }

    //MPS: Constructors --------------------------------------------

    MPS() : sst_(0), minm(1), maxm(MAX_M), cutoff(MAX_CUT) {}

    MPS(const SiteSetT& model,int maxmm = MAX_M, Real cut = MAX_CUT) 
		: N(model.NN()),A(model.NN()+1),left_orth_lim(0),right_orth_lim(model.NN()),
        sst_(&model), minm(1), maxm(maxmm), cutoff(cut)
	{ random_tensors(A); }

    MPS(const SiteSetT& model,const InitState& initState,int maxmm = MAX_M, Real cut = MAX_CUT) 
		: N(model.NN()),A(model.NN()+1),left_orth_lim(0),right_orth_lim(2),
        sst_(&model), minm(1), maxm(maxmm), cutoff(cut)
	{ init_tensors(A,initState);}

    MPS(const SiteSetT& model, istream& s) : N(model.NN()), A(model.NN()+1), sst_(&model)
    { read(s); }

    void read(istream& s)
    {
        for(int j = 1; j <= N; ++j) A[j] = Tensor(s);
        s.read((char*) &left_orth_lim,sizeof(left_orth_lim));
        s.read((char*) &right_orth_lim,sizeof(right_orth_lim));
        s.read((char*) &minm,sizeof(minm));
        s.read((char*) &maxm,sizeof(maxm));
        s.read((char*) &cutoff,sizeof(cutoff));
    }

    void write(ostream& s) const
    {
        //s.write((char*) &N,sizeof(N));
        for(int j = 1; j <= N; ++j) A[j].write(s);
        s.write((char*) &left_orth_lim,sizeof(left_orth_lim));
        s.write((char*) &right_orth_lim,sizeof(right_orth_lim));
        s.write((char*) &minm,sizeof(minm));
        s.write((char*) &maxm,sizeof(maxm));
        s.write((char*) &cutoff,sizeof(cutoff));
    }


    //MPS: operators ------------------------------------------------------

    MPS& operator*=(Real a)
    {
        A[left_orth_lim+1] *= a;
        return *this;
    }
    inline MPS operator*(Real r) const { MPS res(*this); res *= r; return res; }
    friend inline MPS operator*(Real r, MPS res) { res *= r; return res; }

    MPS& operator+=(const MPS& oth);
    inline MPS operator+(MPS res) const { res += *this; return res; }
    //friend inline MPS operator+(MPS A, const MPS& B) { A += B; return A; }
    inline MPS operator-(MPS res) const { res *= -1; res += *this; return res; }

    //MPS: index methods --------------------------------------------------

    void mapprime(int oldp, int newp, PrimeType pt)
	{ for(int i = 1; i <= N; ++i) A[i].mapprime(oldp,newp,pt); }

    void primelinks(int oldp, int newp)
	{ for(int i = 1; i <= N; ++i) A[i].mapprime(oldp,newp,primeLink); }

    void noprimelink()
	{ for(int i = 1; i <= N; ++i) A[i].noprime(primeLink); }

    IndexT LinkInd(int i) const { return index_in_common(A[i],A[i+1],Link); }
    IndexT RightLinkInd(int i) const { assert(i<NN()); return index_in_common(AA(i),AA(i+1),Link); }
    IndexT LeftLinkInd(int i)  const { assert(i>1); return index_in_common(AA(i),AA(i-1),Link); }

    //MPS: orthogonalization methods -------------------------------------

    void doSVD(int i, const Tensor& AA, Direction dir, bool preserve_shape = false)
	{
        /*
        if(!preserve_shape) 
        {
            const int i1 = (dir == Fromleft ? i : i+1);
            const int i2 = (dir == Fromleft ? i-1 : i+2);
            const int boundary = (dir == Fromleft ? 1 : NN()-1);

            list<IndexT> keep_indices;
            IndexT site = A.at(i1).findtype(Site); site.noprime();
            keep_indices.push_back(site);

            int j = A.at(i1).findindexn(site.primed());
            if(j > -1) //part of an MPO
            { keep_indices.push_back(site.primed()); }

            if(i != boundary) 
            {
            foreach(const IndexT& I, AA.indexn())
            if(A.at(i2).hasindexn(I)) keep_indices.push_back(I);
            foreach(const IndexT& I, AA.index1())
            if(A.at(i2).hasindex1(I)) keep_indices.push_back(I);
            }
        }
        //otherwise the MPS will retain the same index structure
        */

        do_denmat_Real(AA,A[i],A[i+1],cutoff,minm,maxm,dir);
        truncerror = svdtruncerr;

        if(dir == Fromleft)
        {
            if(left_orth_lim == i-1 || i == 1) left_orth_lim = i;
            if(right_orth_lim < i+2) right_orth_lim = i+2;
        }
        else
        {
            if(left_orth_lim > i-1) left_orth_lim = i-1;
            if(right_orth_lim == i+2 || i == N-1) right_orth_lim = i+1;
        }
	}

    //Move the orthogonality center to site i (left_orth_lim = i-1, right_orth_lim = i+1)
    void position(int i, bool preserve_shape = false)
	{
        //if(is_null()) { left_orth_lim = i-1; right_orth_lim = i+1; return; }
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
        if(!is_ortho()) Error("MPS: orthogonality center not well defined.");
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
        /*
        Tensor Ajc = AA(j);
        list<IndexT> keep_inds;
        keep_inds.push_back(Ajc.findtype(Site));
        if(dir == Fromleft && j > 1)  keep_inds.push_back(index_in_common(Ajc,AA(j-1),Link));
        if(dir == Fromright && j < N) keep_inds.push_back(index_in_common(Ajc,AA(j+1),Link));
        const bool do_allocate = false;
        IndexT centerlink("c",1); keep_inds.push_back(centerlink);
        AAnc(j) = Tensor(keep_inds,do_allocate); lambda = Tensor(centerlink);
        do_denmat_Real(Ajc,AAnc(j),lambda,cutoff,minm,maxm,dir);
        */

        getCenterMatrix(AAnc(j),(dir == Fromleft ? RightLinkInd(j) : LeftLinkInd(j)),cutoff,minm,maxm,lambda,nameint("c",j));

        if(dir == Fromleft) if(left_orth_lim == j-1 || j == 1) left_orth_lim = j;
        else if(right_orth_lim == j+1 || j == N) right_orth_lim = j;

        if(do_signfix) Error("MPS::getCenter: do_signfix not implemented.");
    }

    Real bondDavidson(int b, Direction dir, const Tensor& mpoh, const Tensor& LH, const Tensor& RH, int niter, int debuglevel, Real errgoal=1E-4);

    void applygate(Tensor gate)
	{
        Tensor AA = A[left_orth_lim+1] * A[left_orth_lim+2] * gate;
        AA.noprime();
        doSVD(left_orth_lim+1,AA,Fromleft);
	}

    Real normalize()
    {
        Real norm2,im;
        psiphi(*this,*this,norm2,im);
        Real norm = sqrt(norm2);
        *this *= 1.0/norm;
        return norm;
    }

    bool is_complex() const
    {
        foreach(const Tensor& AA, A)
        if(AA.is_complex()) return true;
        return false;
    }

    friend inline ostream& operator<<(ostream& s, const MPS& M)
    {
        s << "\n";
        for(int i = 1; i <= M.NN(); ++i) s << M.AA(i) << "\n";
        return s;
    }

    void print(string name = "",Printdat pdat = HideData) const 
    { printdat = (pdat==ShowData); cerr << "\n" << name << " =\n" << *this << "\n"; printdat = false; }

    //-----------------------------------------------------------------
    //IQMPS specific methods


private:
    int collapseCols(const Vector& Diag, Matrix& M) const
    {
        int nr = Diag.Length(), nc = Diag.sumels();
        assert(nr != 0);
        if(nc == 0) return nc;
        M = Matrix(nr,nc); M = 0;
        int c = 0;
        for(int r = 1; r <= nr; ++r)
        if(Diag(r) == 1) { M(r,++c) = 1; }
        return nc;
    }
public:
    template <class IQMPSType> 
    void convertToIQ(IQMPSType& iqpsi) const
    {
        const Real cut = 1E-12;
        assert(sst_ != 0);
        const SiteSetT& sst = *sst_;

        iqpsi = IQMPSType(sst,maxm,cutoff);

        if(!A[1].hasindex(si(1))) Error("convertToIQ: incorrect primelevel for conversion");
        bool is_mpo = A[1].hasindex(si(1).primed());
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

        qC[QN()] = ITensor(); //So that A[1] case executes

        const int show_s = 0;

        Index bond, prev_bond;
        for(int s = 1; s <= N; ++s)
        {
            qD.clear(); qt.clear();
            if(s > 1) prev_bond = LinkInd(s-1); 
            if(s < N) bond = LinkInd(s);

            foreach(const qC_vt& x, qC) {
            const QN& prev_q = x.first; const ITensor& comp = x.second; 
            for(int n = 1; n <= Dim;  ++n)
            for(int u = 1; u <= PDim; ++u)
            {
                q = (is_mpo ? prev_q+si(s).qn(n)-si(s).qn(u) : prev_q-si(s).qn(n));

                block = (s == 1 ? A[s] : conj(comp) * A[s]);
                block *= si(s)(n);
                if(is_mpo) block *= siP(s)(u);

                int count = qD.count(q);
                Vector& D = qD[q];
                if(count == 0) { D.ReDimension(bond.m()); D = 0; }

                bool keep_block = false;
                if(s == N) keep_block = (block.norm() != 0);
                else
                {
                    if(bond.m() == 1) D = 1;
                    else
                    {
                        ITensor summed_block;
                        if(s==1) summed_block = block;
                        else
                        {
                            assert(comp.r()==2);
                            Index new_ind = (comp.index(1)==prev_bond ? comp.index(2) : comp.index(1));
                            summed_block = ITensor(new_ind,1) * block;
                        }
                        //cerr << format("s = %d, bond=")%s << bond << "\n";
                        //summed_block.print("summed_block");

                        Real rel_cut = -1;
                        for(int j = 1; j <= bond.m(); ++j)
                        { rel_cut = max(fabs(summed_block.val1(j)),rel_cut); }
                        assert(rel_cut >= 0);
                        rel_cut *= cut;

                        if(rel_cut > 0)
                        for(int j = 1; j <= bond.m(); ++j)
                        if(fabs(summed_block.val1(j)) > rel_cut) 
                        { D(j) = 1; keep_block = true; }
                    }
                } //if(s != N)

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

            foreach(const qt_vt& x, qt)
            {
                const vector<ITensor>& blks = x.second;
                if(blks.size() != 0)
                {
                    q = x.first; 
                    if(s == N) 
                    { foreach(const ITensor& t, blks) nblock.push_back(t); }
                    else
                    {
                        Matrix M; int mm = collapseCols(qD[q],M);
                        if(s==show_s)
                        {
                            cerr << format("Adding block, mm = %d\n")%mm;
                            q.print("q");
                            cerr << "qD[q] = " << qD[q] << "\n";
                            cerr << "M = \n" << M << "\n";
                            int count = 0;
                            foreach(const ITensor& t, blks) 
                            t.print((format("t%02d")%(++count)).str(),ShowData);
                        }
                        string qname = (format("ql%d (%+d:%d)")%s%q.sz()%q.Nf()).str();
                        Index qbond(qname,mm);
                        ITensor compressor(bond,qbond,M);
                        foreach(const ITensor& t, blks) nblock.push_back(t * compressor);
                        iq.push_back(inqn(qbond,q));
                        qC[q] = compressor;
                    }
                }
            }

            if(s != N) { linkind[s] = IQIndex(nameint("qL",s),iq); iq.clear(); }
            if(s == 1)
                iqpsi.AAnc(s) = (is_mpo ? IQTensor(conj(si(s)),siP(s),linkind[s]) : IQTensor(si(s),linkind[s]));
            else if(s == N)
                iqpsi.AAnc(s) = (is_mpo ? IQTensor(conj(linkind[s-1]),conj(si(s)),siP(s)) 
                                        : IQTensor(conj(linkind[s-1]),si(s)));
            else
                iqpsi.AAnc(s) = (is_mpo ? IQTensor(conj(linkind[s-1]),conj(si(s)),siP(s),linkind[s]) 
                                        : IQTensor(conj(linkind[s-1]),si(s),linkind[s]));

            foreach(const ITensor& nb, nblock) { iqpsi.AAnc(s) += nb; } nblock.clear();

            if(s==show_s)
            {
            iqpsi.AA(s).print((format("qA[%d]")%s).str(),ShowData);
            Error("Stopping");
            }

        } //for loop over s

        IQIndex Center("Center",Index("center",1,Virtual),q,In);
        iqpsi.AAnc(1).viqindex = Center;

    } //void convertToIQ(IQMPSType& iqpsi) const

    operator MPS<IQTensor>() const { MPS<IQTensor> res; convertToIQ(res); return res; }

};

template <class Tensor>
MPS<Tensor>& MPS<Tensor>::operator+=(const MPS<Tensor>& other)
{
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
        AAnc(i) = first[i-1] * AA(i) * first[i] + second[i-1] * other.AA(i) * second[i];
    }
    AAnc(N) = first[N-1] * AA(N) + second[N-1] * other.AA(N);

    noprimelink();

    position(N);
    position(1);

    return *this;
}


template <class Tensor>
Real MPS<Tensor>::bondDavidson(int b, Direction dir, const Tensor& mpoh, const Tensor& LH, const Tensor& RH, int niter, int debuglevel, Real errgoal)
{
    if(b-1 > left_orth_lim)
    {
        cerr << format("MPS::bondDavidson: b=%d, Lb=%d\n")%b%left_orth_lim;
        Error("MPS::bondDavidson: b > left_orth_lim");
    }
    if(b+2 < right_orth_lim)
    {
        cerr << format("MPS::bondDavidson: b+1=%d, Rb=%d\n")%(b+1)%right_orth_lim;
        Error("MPS::bondDavidson: b+1 < right_orth_lim");
    }
    Tensor phi = GET(A,b); phi *= GET(A,b+1);
    Real En = doDavidson(phi,mpoh,LH,RH,niter,debuglevel,errgoal);
    doSVD(b,phi,dir);
    return En;
}


} //namespace Internal
typedef Internal::MPS<ITensor> MPS;
typedef Internal::MPS<IQTensor> IQMPS;

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
                cerr << "check_QNs: Right side Link not pointing In to left of orthog. center at site " << i << endl;
                cerr << "AA(" << i << ") = " << psi.AA(i) << endl;
                return false;
            }
            if(i > 1)
            {
                if(psi.LeftLinkInd(i).dir() != Out) 
                {
                    cerr << "check_QNs: Left side Link not pointing Out to left of orthog. center at site " << i << endl;
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
                cerr << "check_QNs: Right side Link not pointing Out on right side of orthog. center at site " << i << endl;
                return false;
            }
            if(psi.LeftLinkInd(i).dir() != In) 
            {
                cerr << "check_QNs: Left side Link not pointing In on right side of orthog. center at site " << i << endl;
                return false;
            }
        }
    }
    //Done checking arrows

    //Check divergence
    for(int i = 1; i <= N; ++i)
    {
        if(!psi.AA(i).checkDivZero())
        {
            cerr << "check_QNs: IQTensor AA(" << i << ") had non-zero divergence." << endl;
            return false;
        }
    }
    return true;
}

inline QN total_QN(const IQMPS& psi)
{
    assert(psi.AA(psi.ortho_center()).viqindex != IQEmptyV);
    return psi.AA(psi.ortho_center()).viqindex.qn(1); 
}

class IQMPO : public IQMPS
{
    Real lref;
public:
    bool do_relative_cutoff;
    // Norm of psi^2 = 1 = norm = sum of denmat evals. 
    // This translates to Tr{Adag A} = norm.  
    // Ref. norm is Tr{1} = d^N, d = 2 S=1/2, d = 4 for Hubbard, etc
    Real LogReferenceNorm() const { return lref; }
    void setlref(Real _lref) 
	{ 
        //cout << "setting lref to " << setprecision(15) << _lref << endl;
        if(_lref == 0.0)
            Error("bad lref");
        lref = _lref; 
        do_relative_cutoff = false;
	}
    IQMPO() : lref(1.3e-14) { do_relative_cutoff = false; }
    IQMPO(const SiteSet& iss, int maxmm = MAX_M, Real cut = MAX_CUT, Real _lref = -1.0e-17) 
		: IQMPS(iss,maxmm,cut)
	{ 
        if(_lref == 0.0) Error("creating lref to zero");
        if(_lref == -1.0e-17) 
            lref = NN() * log(2.0); 
        else
            lref = _lref;
        do_relative_cutoff = false;
	}
    IQMPO& operator*=(Real a) 
	{
        IQTensor& OC = A[ortho_center()];
        if(!OC.is_null()) OC *= a;
        return *this;
	}
    void primeall()	// sites i,i' -> i',i'';  link:  l -> l'
	{
        for(int i = 1; i <= NN(); i++)
        {
            AAnc(i).mapprime(0,1,primeLink);
            AAnc(i).mapprime(1,2,primeSite);
            AAnc(i).mapprime(0,1,primeSite);
        }
	}
    /*
    bool recalculate_QNs()
	{
        bool did_not_change = true;
        for(int i = 1; i < NN(); i++)
        {
            IQTensor& a(AA(i));
            IQIndex link_prev, linknext(LinkInd(i)); 
            if(i > 1) link_prev = LinkInd(i-1);
            for(IQTensor::iten_it k = a.iten_begin(); k != a.iten_end(); ++k)
            {
                QuantumNumber q;
                Index ilink;
                for(int j = 1; j <= k->d; j++)
                {
                    Index jj(k->index(j));
                    QuantumNumber del = a.qn(jj);
                    if(link_prev.contains(jj))
                    q += del;
                    else if(jj.type() == Site && jj.primelevel == 1)
                    q += del;
                    else if(jj.type() == Site && jj.primelevel == 0)
                    q += -del;
                    else if(jj.type() == Site)
                    Error("bad primelevel on a site");
                    if(linknext.contains(jj))
                    ilink = jj;
                }
                if(q != a.qn(ilink)) 
                {
                    did_not_change = false;
                    a.set_qn(ilink,q);
                    AA(i+1).set_qn(ilink,q);
                }
            }
        }
        return did_not_change;
	}
    */
    void transpose()
	{
        if(A[1].is_null()) return;
        mapprime(1,2,primeSite);
        mapprime(0,1,primeSite);
        mapprime(2,0,primeSite);
        //recalculate_QNs();
        //assert(check_QNs());
	}
};

class MPO : public MPS
{
public:

    MPO() : MPS() {}

    MPO(const SiteSet& iss, int maxmm = MAX_M, Real cut = MAX_CUT)
    : MPS(iss,maxmm,cut)
    {
        maxm = maxmm;
        cutoff = cut;
        vector<Index> hind(N+1);
        for(int i = 1; i <= N; i++) hind[i] = Index(nameint("Hind",i),1);
        A[1] = ITensor(si(1),si(1).primed(),hind[1]);
        A[N] = ITensor(si(N),si(N).primed(),hind[N-1]);
        for(int i = 2; i < N; i++)
            A[i] = ITensor(hind[i-1],si(i),si(i).primed(),hind[i]);
    }

    MPO(SiteSet& iss, istream& s) : MPS(iss,s) { }

    operator IQMPO() const { IQMPO res; convertToIQ(res); return res; }

    MPO& operator*=(Real a)
    {
        A[left_orth_lim+1] *= a;
        return *this;
    }
    inline MPO operator*(Real r) const { MPO res(*this); res *= r; return res; }
    friend inline MPO operator*(Real r, MPO res) { res *= r; return res; }

    //MPO& operator+=(const MPO& oth);
    inline MPO operator+(MPO res) const { res += *this; return res; }
    inline MPO operator-(MPO res) const { res *= -1; res += *this; return res; }

    void newindices()
	{
        vector<Combiner> na(N+1);
        for(int i = 1; i < N; i++)
        {
            Index I = LinkInd(i);
            Index nI(nameint("a",i),I.m());
            na[i] = Combiner(nI,I);
        }
        A[1] = A[1] * na[1];
        for(int i = 2; i < N; i++)
	    { A[i] = na[i-1] * A[i] * na[i]; }
        A[N] = na[N-1] * A[N];
	}

};


namespace Internal {

template <class Tensor>
class HamBuilder
{
public:
    typedef Tensor TensorT;
    typedef typename Tensor::IndexT IndexT;
    typedef SiteSet<IQIndex> SiteSetT;
private:
    const SiteSetT& iss;
    const int N;
    vector<IndexT> currentlinks;
public:

    int NN() const { return N; }
    IndexT si(int i) const { return iss.si(i); }

    HamBuilder(const SiteSetT& iss_) : iss(iss_), N(iss.NN()), currentlinks(N) { }

    void newlinks(vector<Index>& currentlinks)
	{
        currentlinks.resize(N);
        static int ver = 0; ++ver;
        for(int i = 1; i < N; i++)
        {
            stringstream ss;
            ss << "h" << ver << "-" << i;
            currentlinks[i] = IndexT(ss.str(),1);
        }
	}

    void newlinks(vector<IQIndex>& currentlinks)
	{
        currentlinks.resize(N);
        static int ver = 0; ++ver;
        for(int i = 1; i < N; i++)
        {
            stringstream ss;
            ss << "h" << ver << "-" << i;
            currentlinks[i] = IndexT(ss.str());
        }
	}

    Tensor unit(int i) const { return Tensor(si(i),si(i).primed(),1); }

    void getidentity(Real factor, MPO& res)
	{
        newlinks(currentlinks);
        res = MPO(iss,res.maxm,res.cutoff);
        res.AAnc(1) = unit(1); res.AAnc(1).addindex1(GET(currentlinks,1));
        res.AAnc(1) *= factor;
        res.AAnc(N) = unit(N); res.AAnc(N).addindex1(GET(currentlinks,N-1));
        for(int i = 2; i < N; ++i)
        {
            res.AAnc(i) = unit(i); 
            res.AAnc(i).addindex1(GET(currentlinks,i-1));
            res.AAnc(i).addindex1(GET(currentlinks,i));
        }
	}

    void getMPO(Real factor, int i, Tensor op, MPO& res)
	{
        getidentity(1,res);
        res.AAnc(i) = op;
        if(i > 1) res.AAnc(i).addindex1(GET(currentlinks,i-1));
        if(i < N) res.AAnc(i).addindex1(GET(currentlinks,i));
        res *= factor;
	}

    void getMPO(Real factor, int i1, Tensor op1, int i2, Tensor op2, MPO& res)
	{
        if(i1 == i2) Error("HamBuilder::getMPO: i1 cannot equal i2.");
        getMPO(1,i2,op2,res);
        res.AAnc(i1) = op1;
        if(i1 > 1) res.AAnc(i1).addindex1(GET(currentlinks,i1-1));
        if(i1 < N) res.AAnc(i1).addindex1(GET(currentlinks,i1));
        res *= factor;
	}

    template <typename Iterable1, typename Iterable2>
    void getMPO(Real factor, Iterable1 sites, Iterable2 ops, MPO& res)
	{
        for(int i = 0; i < (int) sites.size(); ++i)
        for(int j = 0; j < (int) sites.size(); ++j)
        {
            if(i == j) continue;
            if(sites[i] == sites[j]) Error("HamBuilder::getMPO: all sites should be unique.");
        }

        if(sites.size() != ops.size()) Error("HamBuilder::getMPO: need same number of sites as ops.");

        getidentity(1,res);
        for(int i = 0; i < (int) sites.size(); ++i)
        {
            const int s = sites[i];
            res.AAnc(s) = GET(ops,i);
            assert(GET(ops,i).hasindex(si(sites[i])));
            if(s > 1) res.AAnc(s).addindex1(GET(currentlinks,s-1));
            if(s < N) res.AAnc(s).addindex1(GET(currentlinks,s));
        }
        res *= factor;
	}

};

} //namespace Internal
typedef Internal::HamBuilder<ITensor> HamBuilder;
//typedef Internal::HamBuilder<IQTensor> IQHamBuilder;

class MPOBuilder
{
protected:
    const SiteSet& sst;
    const int Ns;
public:

    MPOBuilder(const SiteSet& sst_) : sst(sst_), Ns(sst_.NN()) { }

    ITensor makeLedge(const Index& L) const
    {
        ITensor res(L); res(L(L.m())) = 1;
        return res;
    }

    ITensor makeRedge(const Index& R) const
    {
        ITensor res(R); res(R(1)) = 1;
        return res;
    }
};

template <class MPSType>
void psiphi(const MPSType& psi, const MPSType& phi, Real& re, Real& im)  // <psi | phi>
{
    const int N = psi.NN();
    if(N != phi.NN()) Error("psiphi: mismatched N");

    typename MPSType::TensorT L = phi.AA(1) * conj(primeind(psi.AA(1),psi.LinkInd(1))); 
    for(int i = 2; i < psi.NN(); ++i) L = L * phi.AA(i) * conj(primelink(psi.AA(i)));
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

template <class MPSType, class MPOType>
void psiHphi(const MPSType& psi, const MPOType& H, const MPSType& phi, Real& re, Real& im) //<psi|H|phi>
{
    typedef typename MPSType::TensorT Tensor;
    const int N = H.NN();

    if(psi.NN() != phi.NN() || psi.NN() != N) Error("psiHphi: mismatched N");

    Tensor L = phi.AA(1) * H.AA(1) * conj(primed(psi.AA(1)));
    for(int i = 2; i < N; ++i) L = L * phi.AA(i) * H.AA(i) * conj(primed(psi.AA(i)));
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

void diag_denmat(const ITensor& rho, Real cutoff, int minm, int maxm, Matrix& U, Vector& D);
void psiHphi(const MPS& psi, const MPO& H, const ITensor& LB, const ITensor& RB, const MPS& phi, Real& re, Real& im);
Real psiHphi(const MPS& psi, const MPO& H, const ITensor& LB, const ITensor& RB, const MPS& phi); //Re[<psi|H|phi>]

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
        //cerr << format("Before +=, cutoff = %.1E, maxm = %d\n")%(res.cutoff)%(res.maxm);
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
int specialdebug1 = 0, specialdebug2 = 0;
#endif

#endif
