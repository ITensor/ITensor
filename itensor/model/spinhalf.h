#ifndef __ITENSOR_SPINHALF_H
#define __ITENSOR_SPINHALF_H
#include "../model.h"

class SpinHalf : public Model
    {
    public:

    SpinHalf(int N);

    IQIndexVal
    Up(int i) const;

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

    int N_;
    std::vector<IQIndex> site_;
        
    };

SpinHalf::
SpinHalf(int N)
    : N_(N),
      site_(N_)
    { }

inline int SpinHalf::
getNN() const
    { return N_; }

inline const IQIndex& SpinHalf::
getSi(int i) const
    { return site_.at(i); }

inline const IQIndex& SpinHalf::
getSiP(int i) const
    { return site_.at(i).primed(); }

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
