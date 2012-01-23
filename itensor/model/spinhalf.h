#ifndef __ITENSOR_SPINHALF_H
#define __ITENSOR_SPINHALF_H
#include "../model.h"

class SpinHalf : public Model
    {
    public:

    SpinHalf();

    SpinHalf(int N);

    SpinHalf(std::ifstream& s) { doRead(s); }

    IQIndexVal
    Up(int i) const;

    IQIndexVal
    Dn(int i) const;

    IQIndexVal
    UpP(int i) const;

    IQIndexVal
    DnP(int i) const;

    private:

    virtual int
    getNN() const;

    virtual const IQIndex&
    getSi(int i) const;

    virtual IQIndex
    getSiP(int i) const;

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
    doRead(std::istream& s);

    virtual void
    doWrite(std::ostream& s) const;

    virtual void
    constructSites();

    //Data members -----------------

    int N_;

    std::vector<IQIndex> site_;
        
    };

inline SpinHalf::
SpinHalf()
    : N_(-1)
    { }

inline SpinHalf::
SpinHalf(int N)
    : N_(N),
      site_(N_+1)
    { 
    constructSites();
    }

inline void SpinHalf::
constructSites()
    {
    for(int j = 1; j <= N_; ++j)
        {
        site_.at(j) = IQIndex(nameint("S=1/2, site=",j),
            Index(nameint("Up for site",j),1,Site),QN(+1,0),
            Index(nameint("Dn for site",j),1,Site),QN(-1,0));
        }
    }

inline void SpinHalf::
doRead(std::istream& s)
    {
    s.read((char*) &N_,sizeof(N_));
    site_.resize(N_+1);
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).read(s);
    }

inline void SpinHalf::
doWrite(std::ostream& s) const
    {
    s.write((char*) &N_,sizeof(N_));
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).write(s);
    }

inline int SpinHalf::
getNN() const
    { return N_; }

inline const IQIndex& SpinHalf::
getSi(int i) const
    { return site_.at(i); }

inline IQIndex SpinHalf::
getSiP(int i) const
    { return site_.at(i).primed(); }

inline IQIndexVal SpinHalf::
Up(int i) const
    {
    return getSi(i)(1);
    }

inline IQIndexVal SpinHalf::
Dn(int i) const
    {
    return getSi(i)(2);
    }

inline IQIndexVal SpinHalf::
UpP(int i) const
    {
    return getSiP(i)(1);
    }

inline IQIndexVal SpinHalf::
DnP(int i) const
    {
    return getSiP(i)(2);
    }

inline IQTensor SpinHalf::
makeSz(int i) const
    {
    IQTensor Sz(si(i),siP(i));
    Sz(Up(i),UpP(i)) = +0.5;
    Sz(Dn(i),DnP(i)) = -0.5;
    return Sz;
    }

inline IQTensor SpinHalf::
makeSx(int i) const
    {
    IQTensor Sx(si(i),siP(i));
    Sx(Up(i),DnP(i)) = +0.5;
    Sx(Dn(i),UpP(i)) = +0.5;
    return Sx;
    }

inline IQTensor SpinHalf::
makeISy(int i) const
    {
    IQTensor ISy(si(i),siP(i));
    ISy(Up(i),DnP(i)) = -0.5;
    ISy(Dn(i),UpP(i)) = +0.5;
    return ISy;
    }

inline IQTensor SpinHalf::
makeSp(int i) const
    {
    IQTensor Sp(si(i),siP(i));
    Sp(Dn(i),UpP(i)) = 1;
    return Sp;
    }

inline IQTensor SpinHalf::
makeSm(int i) const
    {
    IQTensor Sm(si(i),siP(i));
    Sm(Up(i),DnP(i)) = 1;
    return Sm;
    }

#endif
