//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MODEL_H
#define __ITENSOR_MODEL_H
#include "iqtensor.h"

namespace itensor {

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
// values by doing conj(primed(A,Site)) * Op * A.
// (assuming the tensor A is an ortho center 
// of our MPS)
//

class Model
    {
    public:

    typedef std::string
    String;

    typedef std::vector<String>
    DefaultOpsT;

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

    DefaultOpsT
    defaultOps(const OptSet& opts = Global::opts()) const 
        { return getDefaultOps(opts); }

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
    getState(int i, const String& state) const = 0;

    virtual IQTensor
    getOp(int i, const String& opname, const OptSet& opts) const = 0;

    virtual DefaultOpsT
    getDefaultOps(const OptSet& opts) const { return DefaultOpsT(); }

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
    if(opname == "Proj")
        {
        const int n = opts.getInt("State");
        IQIndexVal v = si(i)(n);
        return IQTensor(conj(v),primed(v));
        }
    else
        {
        //If opname of the form "Name1*Name2",
        //return product of Name1 operator times Name2 operator
        const
        std::size_t found = opname.find_first_of('*');
        if(found != std::string::npos)
            {
            try {
            return multSiteOps(getOp(i,op1(opname,found),opts),
                               getOp(i,op2(opname,found),opts));
                }
            catch(const ITError& e)
                {
                Print(opname);
                Print(found);
                Error("More than one * in operator name string?");
                }
            }

        return getOp(i,opname,opts);
        }
    }

std::string inline Model::
op1(const std::string& opname, size_t n) const
    {
    return opname.substr(0,n);
    }

std::string inline Model::
op2(const std::string& opname, size_t n) const
    {
    return opname.substr(n+1);
    }

inline std::ostream& 
operator<<(std::ostream& s, const Model& M)
    {
    s << "Model:\n";
    for(int j = 1; j <= M.N(); ++j) 
        s << boost::format("si(%d) = ")%j << M.si(j) << "\n";
    return s;
    }

}; //namespace itensor

#endif
