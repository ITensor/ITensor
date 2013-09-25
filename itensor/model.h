//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MODEL_H
#define __ITENSOR_MODEL_H
#include "iqtensor.h"

//
// Classes derived from Model 
// represent the abstract lattice of a 
// system as a set of Site indices.
//
// Classes derived from Model are
// responsible for implementing 
// site operators such as Sz for 
// spin models, Cdag for particle
// models, etc. whereas the Model
// base class is reponsible for
// enforcing a consistent interface.
//
// The convention for operators is
// that they are 2-index IQTensors
// with the Site IQIndex pointing
// In and the Site' IQIndex pointing
// Out. This is so we can compute expectation
// values by doing conj(primesite(A)) * Op * A.
// (assuming the tensor A is an ortho center 
// of our MPS)
//

class Model
    {
    public:

    typedef std::string
    String;

    Model() { }

    Model(std::ifstream& s) { }

    //Number of Sites
    int 
    N() const { return getN(); }

    //Index at Site i
    const IQIndex&
    operator()(int i) const { return getSi(i); }

    //Index at Site i, alternate version
    const IQIndex& 
    si(int i) const { return getSi(i); }

    //Primed Index at Site i
    IQIndex 
    siP(int i) const { return primed(getSi(i)); }

    //Index at site i set to a certain state
    //indicated by the string "state"
    //e.g. model("Up",5) returns the IQIndexVal
    //representing the spin up state on site 5
    //(assuming a spin type model such as SpinHalf)
    IQIndexVal
    operator()(int i, const String& state) const
        { return getState(i,state); }

    IQIndexVal
    st(int i, const String& state) const
        { return getState(i,state); }

    IQIndexVal
    stP(int i, const String& state) const
        { return primed(getState(i,state)); }


    //Get the operator indicated by
    //"opname" located at site i
    IQTensor
    op(const String& opname, int i,
       const OptSet& opts = Global::opts()) const;

    void 
    read(std::istream& s) { doRead(s); }

    void 
    write(std::ostream& s) const { doWrite(s); }

    virtual 
    ~Model() { }

    //Implementations (To Be Overridden by Derived Classes) 

    private:

    virtual int
    getN() const = 0;

    virtual const IQIndex&
    getSi(int i) const = 0;

    virtual IQIndexVal
    getState(int i, const String& state) const;

    virtual IQTensor
    getOp(int i, const String& opname, 
          const OptSet& opts = Global::opts()) const;

    protected:

    virtual void
    doRead(std::istream& s) = 0;

    virtual void
    doWrite(std::ostream& s) const = 0;

    //
    //
    // Older interface (deprecated)
    // Can safely ignore everything below;
    // present only for backwards compatibility.
    //
    //

    public:

    IQTensor 
    id(int i) const { return op("Id",i); }

    IQTensor 
    proj(int i,int n) const;

    IQTensor
    tReverse(int i) const { return makeTReverse(i); }

    IQTensor 
    sx(int i) const { return makeSx(i); }

    IQTensor 
    isy(int i) const { return makeISy(i); }

    IQTensor 
    sz(int i) const { return makeSz(i); }

    IQTensor 
    sp(int i) const { return makeSp(i); }

    IQTensor 
    sm(int i) const { return makeSm(i); }

    IQTensor
    sz2(int i) const { return makeSz2(i); }

    IQTensor
    sx2(int i) const { return makeSx2(i); }

    IQTensor
    sy2(int i) const { return makeSy2(i); }

    IQTensor
    n(int i) const { return makeN(i); }

    IQTensor
    C(int i) const { return makeC(i); }

    IQTensor
    Cdag(int i) const { return makeCdag(i); }

    IQTensor
    A(int i) const { return makeA(i); }

    IQTensor
    Adag(int i) const { return makeAdag(i); }

    IQTensor
    fermiPhase(int i) const { return makeFermiPhase(i); }

    IQTensor
    Nup(int i) const { return makeNup(i); }

    IQTensor
    Ndn(int i) const { return makeNdn(i); }

    IQTensor
    Nupdn(int i) const { return makeNupdn(i); }

    IQTensor
    Ntot(int i) const { return makeNtot(i); }

    IQTensor
    Cup(int i) const { return makeCup(i); }

    IQTensor
    Cdagup(int i) const { return makeCdagup(i); }

    IQTensor
    Cdn(int i) const { return makeCdn(i); }

    IQTensor
    Cdagdn(int i) const { return makeCdagdn(i); }

    IQTensor
    Aup(int i) const { return makeAup(i); }

    IQTensor
    Adagup(int i) const { return makeAdagup(i); }

    IQTensor
    Adn(int i) const { return makeAdn(i); }

    IQTensor
    Adagdn(int i) const { return makeAdagdn(i); }

    //
    // Implementations for older interface (deprecated)
    //
    private:

    virtual IQTensor 
    makeProj(int i, int n) const { return proj(i,n); }

    virtual IQTensor 
    makeTReverse(int i) const { return getOp(i,"tReverse"); }

    virtual IQTensor 
    makeSx(int i) const { return getOp(i,"Sx"); }

    virtual IQTensor 
    makeISy(int i) const { return getOp(i,"ISy"); }

    virtual IQTensor 
    makeSz(int i) const { return getOp(i,"Sz"); }

    virtual IQTensor 
    makeSp(int i) const { return getOp(i,"Sp"); }

    virtual IQTensor 
    makeSm(int i) const { return getOp(i,"Sm"); }

    virtual IQTensor
    makeSz2(int i) const { return getOp(i,"Sz2"); }

    virtual IQTensor
    makeSx2(int i) const { return getOp(i,"Sx2"); }

    virtual IQTensor
    makeSy2(int i) const { return getOp(i,"Sy2"); }

    virtual IQTensor 
    makeN(int i) const { return getOp(i,"N"); }

    virtual IQTensor 
    makeC(int i) const { return getOp(i,"C"); }

    virtual IQTensor 
    makeCdag(int i) const { return getOp(i,"Cdag"); }

    virtual IQTensor 
    makeA(int i) const { return getOp(i,"A"); }

    virtual IQTensor 
    makeAdag(int i) const { return getOp(i,"Adag"); }

    virtual IQTensor 
    makeFermiPhase(int i) const { return getOp(i,"FermiPhase"); }

    virtual IQTensor 
    makeNup(int i) const { return getOp(i,"Nup"); }

    virtual IQTensor 
    makeNdn(int i) const { return getOp(i,"Ndn"); }

    virtual IQTensor 
    makeNupdn(int i) const { return getOp(i,"Nupdn"); }

    virtual IQTensor 
    makeNtot(int i) const { return getOp(i,"Ntot"); }

    virtual IQTensor 
    makeCup(int i) const { return getOp(i,"Cup"); }

    virtual IQTensor 
    makeCdagup(int i) const { return getOp(i,"Cdagup"); }

    virtual IQTensor 
    makeCdn(int i) const { return getOp(i,"Cdn"); }

    virtual IQTensor 
    makeCdagdn(int i) const { return getOp(i,"Cdagdn"); }

    virtual IQTensor 
    makeAup(int i) const { return getOp(i,"Aup"); }

    virtual IQTensor 
    makeAdagup(int i) const { return getOp(i,"Adagup"); }

    virtual IQTensor 
    makeAdn(int i) const { return getOp(i,"Adn"); }

    virtual IQTensor 
    makeAdagdn(int i) const { return getOp(i,"Adagdn"); }


    };

inline IQTensor Model::
op(const String& opname, int i, 
   const OptSet& opts) const
    { 
    if(opname == "Id")
        {
        IQIndex s = conj(si(i));
        IQIndex sP = siP(i);
        IQTensor id_(s,sP);
        for(int j = 1; j <= s.m(); ++j)
            {
            id_(s(j),sP(j)) = 1;
            }
        return id_;
        }
    else
        {
        return getOp(i,opname,opts);
        }
    }

IQIndexVal inline Model::
getState(int i, const String& state) const
    {
    Error("This Model class only supports older interface - getState not defined.");
    return IQIndexVal();
    }

IQTensor inline Model::
getOp(int i, const String& opname, 
      const OptSet& opts) const
    {
    Error("This Model class only supports older interface - getOp not defined.");
    return IQTensor();
    }

IQTensor inline Model::
proj(int i, int n) const
    {
    IQIndex s = conj(si(i));
    IQIndex sP = siP(i);
    IQTensor proj_(s,sP);
    proj_(s(n),sP(n)) = 1;
    return proj_;
    }

inline std::ostream& 
operator<<(std::ostream& s, const Model& M)
    {
    s << "Model:\n";
    for(int j = 1; j <= M.N(); ++j) 
        s << boost::format("si(%d) = ")%j << M.si(j) << "\n";
    return s;
    }

#endif
