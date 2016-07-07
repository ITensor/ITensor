//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SITESET_H
#define __ITENSOR_SITESET_H
#include "itensor/iqtensor.h"

namespace itensor {

//
// Classes derived from SiteSet 
// represent the Hilbert space of a 
// system as a set of Site indices.
//
// Classes derived from SiteSet are
// responsible for implementing 
// site operators such as Sz for 
// spin models, Cdag for particle
// models, etc. whereas the SiteSet
// base class is reponsible for
// enforcing a consistent interface.
//
// The convention for operators is
// that they are 2-index IQTensors
// with the Site IQIndex pointing
// In and the Site' IQIndex pointing
// Out. This is so we can compute expectation
// values by doing dag(prime(A,Site)) * Op * A.
// (assuming the tensor A is an ortho center 
// of our MPS)
//

class SiteSet;
using Model = SiteSet; //for backwards compatibility

class SiteSet
    {
    public:

    using String = std::string;
    using DefaultOpsT = std::vector<String>;

    SiteSet() { }

    SiteSet(std::ifstream& s) { }

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
    siP(int i) const { return prime(getSi(i)); }

    //Index at site i set to a certain state
    //indicated by the string "state"
    //e.g. sites("Up",5) returns the IQIndexVal
    //representing the spin up state on site 5
    //(assuming a spin SiteSet such as SpinHalf)
    IQIndexVal
    operator()(int i, const String& state) const
        { return getState(i,state); }

    IQIndexVal
    st(int i, const String& state) const
        { return getState(i,state); }

    IQIndexVal
    stP(int i, const String& state) const
        { return prime(getState(i,state)); }


    //Get the operator indicated by
    //"opname" located at site i
    IQTensor
    op(const String& opname, int i,
       const Args& args = Global::args()) const;

    DefaultOpsT
    defaultOps(const Args& args = Global::args()) const 
        { return getDefaultOps(args); }

    void 
    read(std::istream& s) { doRead(s); }

    void 
    write(std::ostream& s) const { doWrite(s); }

    virtual 
    ~SiteSet() { }

    //Implementations (To Be Overridden by Derived Classes) 

    private:

    virtual int
    getN() const = 0;

    virtual const IQIndex&
    getSi(int i) const = 0;

    virtual IQIndexVal
    getState(int i, const String& state) const = 0;

    virtual IQTensor
    getOp(int i, const String& opname, const Args& args) const = 0;

    virtual DefaultOpsT
    getDefaultOps(const Args& args) const { return DefaultOpsT(); }

    protected:

    virtual void
    doRead(std::istream& s) = 0;

    virtual void
    doWrite(std::ostream& s) const = 0;

    private:

    std::string
    op1(const std::string& opname, size_t n) const;

    std::string
    op2(const std::string& opname, size_t n) const;

    };

inline IQTensor SiteSet::
op(const String& opname, int i, 
   const Args& args) const
    { 
    if(opname == "Id")
        {
        IQIndex s = dag(si(i));
        IQIndex sP = siP(i);
        IQTensor id_(s,sP);
        for(int j = 1; j <= s.m(); ++j)
            {
            id_(s(j),sP(j)) = 1;
            }
        return id_;
        }
    else
    if(opname == "Proj")
        {
        const int n = args.getInt("State");
        IQIndexVal v = si(i)(n);
        return IQTensor(dag(v),prime(v));
        }
    else
        {
        //If opname of the form "Name1*Name2",
        //return product of Name1 operator times Name2 operator
        auto found = opname.find_first_of('*');
        if(found != std::string::npos)
            {
            try {
            return multSiteOps(getOp(i,op1(opname,found),args),
                               op(op2(opname,found),i,args));
                }
            catch(const ITError& e)
                {
                Print(opname);
                Print(found);
                Error("More than one * in operator name string?");
                }
            }

        return getOp(i,opname,args);
        }
    }

std::string inline SiteSet::
op1(const std::string& opname, size_t n) const
    {
    return opname.substr(0,n);
    }

std::string inline SiteSet::
op2(const std::string& opname, size_t n) const
    {
    return opname.substr(n+1);
    }

inline std::ostream& 
operator<<(std::ostream& s, const SiteSet& M)
    {
    s << "SiteSet:\n";
    for(int j = 1; j <= M.N(); ++j) 
        s << format("si(%d) = ",j,M.si(j),"\n");
    return s;
    }

} //namespace itensor

#endif
