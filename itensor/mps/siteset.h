//
// Distributed under the ITensor Library License, Version 1.2
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
    using storage = std::vector<IQIndex>;
    std::shared_ptr<storage> sites_;
    public:

    using String = std::string;
    using DefaultOpsT = std::vector<String>;

    SiteSet() { }

    SiteSet(int N) 
        : sites_{std::make_shared<storage>(N+1)}
        { }

    explicit operator bool() const { return bool(sites_); }

    //Number of Sites
    int 
    N() const { return getN(); }

    //Index at Site i
    IQIndex const&
    operator()(int i) const { return getSi(i); }

    //Index at Site i, alternate version
    IQIndex const& 
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
    operator()(int i, String const& state) const
        { return getState(i,state); }

    IQIndexVal
    st(int i, String const& state) const
        { return getState(i,state); }

    IQIndexVal
    stP(int i, String const& state) const
        { return prime(getState(i,state)); }


    //Get the operator indicated by
    //"opname" located at site i
    IQTensor
    op(String const& opname, int i,
       Args const& args = Global::args()) const;

    DefaultOpsT
    defaultOps(Args const& args = Args::global()) const 
        { return getDefaultOps(args); }

    void 
    read(std::istream & s);

    void 
    write(std::ostream & s) const;

    virtual 
    ~SiteSet() { }

    //Index at Site i, alternate version
    void
    set(int i, IQIndex I) 
        { 
        if(!*this) Error("Cannot call set on default-initialized SiteSet");
        if(not sites_.unique()) 
            {
            sites_ = std::make_shared<storage>(*sites_);
            }
        sites_->at(i) = I;
        }

    private:

    //Implementations (To Be Overridden by Derived Classes) 

    virtual int
    getN() const { return sites_ ? static_cast<int>(sites_->size())-1 : 0; }

    virtual IQIndex const&
    getSi(int i) const 
        { 
        if(not *this) Error("Default initialized SiteSet");
        return (*sites_).at(i); 
        }

    virtual IQIndexVal
    getState(int i, String const& state) const
        {
        Error("getState not defined SiteSet or derived class");
        return IQIndexVal();
        }

    virtual IQTensor
    getOp(int i, 
          String const& opname, 
          Args const& args) const
        {
        Error("getOp not defined in SiteSet or derived class");
        return IQTensor{};
        }

    virtual DefaultOpsT
    getDefaultOps(Args const& args) const 
        { 
        return DefaultOpsT(); 
        }

    protected:

    virtual void
    doRead(std::istream& s) { }

    virtual void
    doWrite(std::ostream& s) const { }

    private:

    std::string
    op1(std::string const& opname, size_t n) const;

    std::string
    op2(std::string const& opname, size_t n) const;

    public:

    //for backwards compatibility with ITensor v1
    SiteSet(std::ifstream& s) { }

    };

inline IQTensor SiteSet::
op(String const& opname, int i, 
   Args const& args) const
    { 
    if(not *this) Error("Cannot call .op(..) on default-initialized SiteSet");
    if(opname == "Id")
        {
        IQIndex s = dag(si(i));
        IQIndex sP = siP(i);
        auto id_ = IQTensor(s,sP);
        for(int j = 1; j <= s.m(); ++j)
            {
            id_.set(s(j),sP(j),1);
            }
        return id_;
        }
    else
    if(opname == "Proj")
        {
        auto n = args.getInt("State");
        auto v = si(i)(n);
        return setElt(dag(v),prime(v));
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
                               getOp(i,op2(opname,found),args));
                }
            catch(ITError const& e)
                {
                println("opname = ",opname);
                printfln("found = %s",found);
                Error("More than one * in operator name string?");
                }
            }

        return getOp(i,opname,args);
        }
    }

std::string inline SiteSet::
op1(std::string const& opname, size_t n) const
    {
    return opname.substr(0,n);
    }

std::string inline SiteSet::
op2(std::string const& opname, size_t n) const
    {
    return opname.substr(n+1);
    }

void inline SiteSet::
read(std::istream & s)
    {
    int N = 0;
    s.read((char*) &N,sizeof(N));
    if(N > 0)
        {
        sites_ = std::make_shared<storage>(N+1);
        for(int j = 1; j <= N; ++j) 
            {
            sites_->at(j).read(s);
            }
        }
    doRead(s);
    }

void inline SiteSet::
write(std::ostream & s) const
    {
    auto N = getN();
    s.write((char*) &N,sizeof(N));
    if(sites_)
        {
        for(int j = 1; j <= N; ++j) 
            {
            sites_->at(j).write(s);
            }
        }
    doWrite(s);
    }

inline std::ostream& 
operator<<(std::ostream& s, SiteSet const& M)
    {
    s << "SiteSet:\n";
    for(int j = 1; j <= M.N(); ++j) 
        {
        s << format("si(%d) = %s\n",j,M.si(j));
        }
    return s;
    }

} //namespace itensor


#endif
