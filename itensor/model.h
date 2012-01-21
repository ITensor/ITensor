#ifndef __ITENSOR_MODEL_H
#define __ITENSOR_MODEL_H
#include "iqtensor.h"

class Model
    {
    public:

    Model() { }

    //Number of Sites
    int 
    NN() const { return getNN(); }

    //Index at Site i
    const IQIndex& 
    si(int i) const { return getSi(i); }

    //Primed Index at Site i
    IQIndex 
    siP(int i) const { return getSiP(i); }

    //Operators -----------------------

    IQTensor 
    id(int i) const { return makeId(i); }

    //Spin Operators -----------------------

    IQTensor 
    sx(int i) const { makeSx(i); }

    IQTensor 
    isy(int i) const { makeISy(i); }

    IQTensor 
    sz(int i) const { makeSz(i); }

    IQTensor 
    sp(int i) const { makeSp(i); }

    IQTensor 
    sm(int i) const { makeSm(i); }

    //Particle Operators -----------------------

    IQTensor
    n(int i) const { makeN(i); }

    IQTensor
    C(int i) const { makeC(i); }

    IQTensor
    Cdag(int i) const { makeCdag(i); }

    IQTensor
    fermiPhase(int i) const { makeFermiPhase(i); }

    //Hubbard Model Operators -----------

    IQTensor
    Nup(int i) const { makeNUp(i); }

    IQTensor
    Ndn(int i) const { makeNDn(i); }

    IQTensor
    Nupdn(int i) const { makeNDn(i); }

    IQTensor
    Ntot(int i) const { makeNTot(i); }

    IQTensor
    Cup(int i) const { makeCup(i); }

    IQTensor
    Cdagup(int i) const { makeCdagup(i); }

    IQTensor
    Cdn(int i) const { makeCdn(i); }

    IQTensor
    Cdagdn(int i) const { makeCdagdn(i); }

    //Other Methods -----------------------

    void 
    read(std::istream& s) { doRead(s); }

    void 
    write(std::ostream& s) const { doWrite(s); }

    virtual 
    ~Model() { }



    //Implementations (To Be Overridden by Derived Classes) 

    private:

    virtual int
    getNN() const = 0;

    virtual const IQIndex&
    getSi(int i) const = 0;

    virtual IQIndex
    getSiP(int i) const = 0;

    virtual IQTensor 
    makeId(int i) const;

    virtual IQTensor 
    makeSx(int i) const { Error("sx not implemented"); }

    virtual IQTensor 
    makeISy(int i) const { Error("isy not implemented"); }

    virtual IQTensor 
    makeSz(int i) const { Error("sz not implemented"); }

    virtual IQTensor 
    makeSp(int i) const { Error("sp not implemented"); }

    virtual IQTensor 
    makeSm(int i) const { Error("sm not implemented"); }

    virtual IQTensor 
    makeN(int i) const { Error("n not implemented"); }

    virtual IQTensor 
    makeC(int i) const { Error("C not implemented"); }

    virtual IQTensor 
    makeCdag(int i) const { Error("Cdag not implemented"); }

    virtual IQTensor 
    makeFermiPhase(int i) const { Error("fermiPhase not implemented"); }

    virtual IQTensor 
    makeNup(int i) const { Error("Nup not implemented"); }

    virtual IQTensor 
    makeNdn(int i) const { Error("Ndn not implemented"); }

    virtual IQTensor 
    makeNupdn(int i) const { Error("Nupdn not implemented"); }

    virtual IQTensor 
    makeNtot(int i) const { Error("Ntot not implemented"); }

    virtual IQTensor 
    makeCup(int i) const { Error("Cup not implemented"); }

    virtual IQTensor 
    makeCdagup(int i) const { Error("Cdagup not implemented"); }

    virtual IQTensor 
    makeCdn(int i) const { Error("Cdn not implemented"); }

    virtual IQTensor 
    makeCdagdn(int i) const { Error("Cdagdn not implemented"); }

    protected:

    virtual void
    doRead(std::istream& s) = 0;

    virtual void
    doWrite(std::ostream& s) const = 0;

    };

inline IQTensor Model::
makeId(int i) const
    { 
    IQTensor id_(si(i),siP(i));
    for(int j = 1; j <= si(i).m(); ++j)
        id_(si(i)(j),siP(i)(j)) = 1;
    return id_;
    }

inline std::ostream& 
operator<<(std::ostream& s, const Model& M)
    {
    s << "Model:\n";
    for(int j = 1; j <= M.NN(); ++j) 
        s << boost::format("si(%d) = ")%j << M.si(j) << "\n";
    return s;
    }

#endif
