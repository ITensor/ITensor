#ifndef __MODEL_H
#define __MODEL_H
#include "iqtensor.h"

class BaseModel
{
public:
    virtual int dim() const = 0;
    virtual int NN() const = 0;
    virtual IQIndex si(int i) const = 0;
    virtual IQIndex siP(int i) const = 0;

    virtual void read(std::istream& s) = 0;
    virtual void write(std::ostream& s) const = 0;

    virtual SiteOp id(int i) const = 0;
protected:
    virtual ~BaseModel() { }
};
typedef BaseModel SiteSet;

inline std::ostream& operator<<(std::ostream& s, const BaseModel& b)
{
    s << "Model:\n";
    for(int j = 1; j <= b.NN(); ++j) s << boost::format("si(%d) = ")%j << b.si(j) << "\n";
    return s;
}

//---------------------------------------------------------
//Definition of Model Types
//---------------------------------------------------------

namespace SpinOne {

const int Dim = 3;

class Model : public BaseModel
{
    typedef BaseModel Parent;

    int N;
    vector<IQIndex> site;
public:
    Model() : N(-1) { }
    Model(int N_) : N(N_), site(N_+1) 
    { 
        for(int i = 1; i <= N; ++i)
        {
        site.at(i) = IQIndex(nameint("S=1, site=",i),
        Index(nameint("Up for site",i),1,Site),QN(+2,0),
        Index(nameint("Z0 for site",i),1,Site),QN( 0,0),
        Index(nameint("Dn for site",i),1,Site),QN(-2,0));
        }
    }
    Model(std::istream& s) { read(s); }

    virtual ~Model() { }

    void read(std::istream& s)
    {
        s.read((char*) &N,sizeof(N));
        site.resize(N+1);
        for(int j = 1; j <= N; ++j) site.at(j).read(s);
    }
    void write(std::ostream& s) const
    {
        s.write((char*) &N,sizeof(N));
        for(int j = 1; j <= N; ++j) site.at(j).write(s);
    }

    int dim() const { return SpinOne::Dim; }
    inline int NN() const { return N; }
    inline IQIndex si(int i) const { return site.at(i); }
    inline IQIndex siP(int i) const { return site.at(i).primed(); }

    IQIndexVal Up(int i) const { return si(i)(1); }
    IQIndexVal Z0(int i) const { return si(i)(2); }
    IQIndexVal Dn(int i) const { return si(i)(3); }

    IQIndexVal UpP(int i) const { return siP(i)(1); }
    IQIndexVal Z0P(int i) const { return siP(i)(2); }
    IQIndexVal DnP(int i) const { return siP(i)(3); }

    virtual SiteOp id(int i) const
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

class Model : public BaseModel
{
    typedef BaseModel Parent;
    
    int N;
    vector<IQIndex> site;
public:
    Model() : Parent() { }
    Model(int N_) : N(N_), site(N_+1) 
    {
        for(int i = 1; i <= N; ++i)
        {
        site.at(i) = IQIndex(nameint("S=1/2, site=",i),
        Index(nameint("Up for site",i),1,Site),QN(+1,0),
        Index(nameint("Dn for site",i),1,Site),QN(-1,0));
        }
    }
    Model(std::istream& s) { read(s); }

    virtual ~Model() { }

    void read(std::istream& s)
    {
        s.read((char*) &N,sizeof(N));
        site.resize(N+1);
        for(int j = 1; j <= N; ++j) site.at(j).read(s);
    }
    void write(std::ostream& s) const
    {
        s.write((char*) &N,sizeof(N));
        for(int j = 1; j <= N; ++j) site.at(j).write(s);
    }

    int dim() const { return SpinHalf::Dim; }
    int NN() const { return N; }
    IQIndex si(int i) const { return site.at(i); }
    IQIndex siP(int i) const { return site.at(i).primed(); }

    IQIndexVal Up(int i) const { return si(i)(1); }
    IQIndexVal Dn(int i) const { return si(i)(2); }

    IQIndexVal UpP(int i) const { return siP(i)(1); }
    IQIndexVal DnP(int i) const { return siP(i)(2); }

    virtual SiteOp id(int i) const
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

namespace Spinless {

const int Dim = 2;

class Model : public BaseModel
{
public:
    bool odd_even_up_down,conserve_Nf;
private:
    typedef BaseModel Parent;

    int N;
    vector<IQIndex> site;

    void initSites()
    {
        int occ = (conserve_Nf ? 1 : 0);
        if(odd_even_up_down)
        {
        for(int i = 1; i <= N; ++i)
        {
            if(i%2==1)
            {
                site.at(i) = IQIndex(nameint("Spinless, Up site=",i),
                Index(nameint("Emp for Up site",i),1,Site),QN(0,0,0),
                Index(nameint("Occ for Up site",i),1,Site),QN(+1,occ,1));
            }
            else
            {
                site.at(i) = IQIndex(nameint("Spinless, Dn site=",i),
                Index(nameint("Emp for Dn site",i),1,Site),QN(0,0,0),
                Index(nameint("Occ for Dn site",i),1,Site),QN(-1,occ,1));
            }
        }
        }
        else
        {
        for(int i = 1; i <= N; ++i)
        {
            site.at(i) = IQIndex(nameint("Spinless, site=",i),
            Index(nameint("Emp for site",i),1,Site),QN(0,0,0),
            Index(nameint("Occ for site",i),1,Site),QN(0,occ,1));
        }
        }
    }
public:
    Model() : odd_even_up_down(false),conserve_Nf(true), N(-1) { }
    Model(int N_, bool odd_even_up_down_ = false, bool conserve_Nf_ = true) 
    : odd_even_up_down(odd_even_up_down_), conserve_Nf(conserve_Nf_), N(N_), site(N_+1)
    { initSites(); }
    Model(std::istream& s) { read(s); }

    virtual ~Model() { }

    void read(std::istream& s)
    { 
        s.read((char*) &odd_even_up_down,sizeof(odd_even_up_down));
        s.read((char*) &conserve_Nf,sizeof(conserve_Nf));
        s.read((char*) &N,sizeof(N));
        site.resize(N+1);
        for(int j = 1; j <= N; ++j) site.at(j).read(s);
    }
    void write(std::ostream& s) const
    {
        s.write((char*) &odd_even_up_down,sizeof(odd_even_up_down));
        s.write((char*) &conserve_Nf,sizeof(conserve_Nf));
        s.write((char*) &N,sizeof(N));
        for(int j = 1; j <= N; ++j) site.at(j).write(s);
    }

    int dim() const { return Spinless::Dim; }
    int NN() const { return N; }
    IQIndex si(int i) const { return GET(site,i); }
    IQIndex siP(int i) const { return GET(site,i).primed(); }

    IQIndexVal Emp(int i) const { return si(i)(1); }
    IQIndexVal Occ(int i) const { return si(i)(2); }

    IQIndexVal EmpP(int i) const { return siP(i)(1); }
    IQIndexVal OccP(int i) const { return siP(i)(2); }

    virtual SiteOp id(int i) const
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

class Model : public BaseModel
{
private:
    typedef BaseModel Parent;

    int N;
    vector<IQIndex> site;
public:
    Model() : N(-1) { }
    Model(int N_) : N(N_), site(N_+1)
    {
        for(int i = 1; i <= N; ++i)
	    {
	    site.at(i) = IQIndex(nameint("Hubbard, site=",i),
		    Index(nameint("Emp for site ",i),1,Site),  QN( 0,0,0),
		    Index(nameint("Up for site ",i),1,Site),   QN(+1,1,1),
		    Index(nameint("Dn for site ",i),1,Site),   QN(-1,1,1),
		    Index(nameint("Up-Dn for site ",i),1,Site),QN( 0,2,0));
	    }
    }
    Model(std::istream& s) { read(s); }

    virtual ~Model() { }

    void read(std::istream& s)
    {
        s.read((char*) &N,sizeof(N));
        site.resize(N+1);
        for(int j = 1; j <= N; ++j) site.at(j).read(s);
    }
    void write(std::ostream& s) const
    {
        s.write((char*) &N,sizeof(N));
        for(int j = 1; j <= N; ++j) site.at(j).write(s);
    }

    int dim() const { return Hubbard::Dim; }
    int NN() const { return N; }
    IQIndex si(int i) const { return GET(site,i); }
    IQIndex siP(int i) const { return GET(site,i).primed(); }

    IQIndexVal Emp(int i) const { return si(i)(1); }
    IQIndexVal UpState(int i) const { return si(i)(2); }
    IQIndexVal DnState(int i) const { return si(i)(3); }
    IQIndexVal UpDnState(int i) const { return si(i)(4); } // cdag_dn cdag_up | vac >

    IQIndexVal EmpP(int i) const { return siP(i)(1); }
    IQIndexVal UpStateP(int i) const { return siP(i)(2); }
    IQIndexVal DnStateP(int i) const { return siP(i)(3); }
    IQIndexVal UpDnStateP(int i) const { return siP(i)(4); }

    virtual SiteOp id(int i) const
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
	res(Emp(i),UpStateP(i)) = 1;
	res(DnState(i),UpDnStateP(i)) = -1;
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
	res(Emp(i),DnStateP(i)) = 1;
	res(UpState(i),UpDnStateP(i)) = 1;
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


#endif
