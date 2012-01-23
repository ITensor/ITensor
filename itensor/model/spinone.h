#ifndef __ITENSOR_SPINONE_H
#define __ITENSOR_SPINONE_H
#include "../model.h"

class SpinOne : public Model
    {
    public:

    SpinOne(int N);

    IQIndexVal
    Up(int i) const;

    IQIndexVal
    Z0(int i) const;

    IQIndexVal
    Dn(int i) const;

    private:

    virtual int
    getNN() const;

    virtual const IQIndex&
    getSi() const;

    virtual const IQIndex&
    getSiP() const;

    virtual IQTensor
    makeSz(int i) const;

    virtual IQTensor
    makeSx(int i) const;

    virtual IQTensor
    makeISy(int i) const;

    virtual IQTensor
    makeSp(int i) const;

    virtual IQTensor
    makeSm(int i) const;

    virtual void
    constructSites();
        
    //Data members -----------------

    int N_;

    std::vector<IQIndex> site_;

    };

SpinOne::
SpinOne(int N)
    : N_(N),
      site_(N_)
    { 
    constructSites();
    }

inline void SpinOne::
constructSites()
    {
    for(int j = 1; j <= N_; ++j)
        {
        site_.at(j) = IQIndex(nameint("S=1, site=",i),
            Index(nameint("Up for site",i),1,Site),QN(+2,0),
            Index(nameint("Z0 for site",i),1,Site),QN( 0,0),
            Index(nameint("Dn for site",i),1,Site),QN(-2,0));
        }
    }

inline void SpinOne::
doRead(std::istream& s)
    {
    s.read((char*) &N_,sizeof(N_));
    site.resize(N_+1);
    for(int j = 1; j <= N_; ++j) 
        site.at(j).read(s);
    }

inline void SpinOne::
doWrite(std::ostream& s) const
    {
    s.write((char*) &N,sizeof(N));
    for(int j = 1; j <= N_; ++j) 
        site.at(j).write(s);
    }

inline int SpinOne::
getNN() const
    { return N_; }

inline const IQIndex& SpinOne::
getSi(int i) const
    { return site_.at(i); }

inline const IQIndex& SpinOne::
getSiP(int i) const
    { return site_.at(i).primed(); }

inline IQTensor SpinOne::
makeSz(int i) const
    {
    IQTensor Sz(si(i),siP(i));
    Sz(Up(i),UpP(i)) = +1.;
    Sz(Dn(i),DnP(i)) = -1.;
    return Sz;
    }

inline IQTensor SpinOne::
makeSx(int i) const
    {
    IQTensor Sx(si(i),siP(i));
    Sx(Up(i),Z0P(i)) = ISqrt2; 
    Sx(Z0(i),UpP(i)) = ISqrt2;
    Sx(Z0(i),DnP(i)) = ISqrt2; 
    Sx(Dn(i),Z0P(i)) = ISqrt2;
    return Sx;
    }

inline IQTensor SpinOne::
makeISy(int i) const
    {
    IQTensor ISy(si(i),siP(i));
    ISy(Up(i),Z0P(i)) = +ISqrt2; 
    ISy(Z0(i),UpP(i)) = -ISqrt2;
    ISy(Z0(i),DnP(i)) = +ISqrt2; 
    ISy(Dn(i),Z0P(i)) = -ISqrt2;
    return ISy;
    }

inline IQTensor SpinOne::
makeSp(int i) const
    {
    IQTensor Sp(si(i),siP(i));
    Sp(Dn(i),Z0P(i)) = Sqrt2; 
    Sp(Z0(i),UpP(i)) = Sqrt2;
    return Sp;
    }

inline IQTensor SpinOne::
makeSm(int i) const
    {
    IQTensor Sm(si(i),siP(i));
    Sm(Up(i),Z0P(i)) = Sqrt2; 
    Sm(Z0(i),DnP(i)) = Sqrt2;
    return Sm;
    }

#endif
