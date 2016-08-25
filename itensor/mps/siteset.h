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
// The convention for operators is
// that they are 2-index IQTensors
// with the Site IQIndex pointing
// In and the Site' IQIndex pointing
// Out. This is so we can compute expectation
// values by doing dag(prime(A,Site)) * Op * A.
// (assuming the tensor A is an ortho center 
// of our MPS)
//

class GenericSite;
struct BaseSiteStore;

class SiteSet
    {
    std::shared_ptr<BaseSiteStore> sites_;
    public:

    using String = std::string;

    SiteSet() { }

    explicit operator bool() const { return bool(sites_); }

    int 
    N() const;

    //Index at Site i
    IQIndex
    operator()(int i) const;

    //Index at site i set to a certain state
    //indicated by the string "state"
    //e.g. sites(5,"Up") returns the IQIndexVal
    //representing the spin up state on site 5
    //(assuming a spin SiteSet such as SpinHalf)
    IQIndexVal
    operator()(int i, String const& state) const;

    //Get the operator indicated by
    //"opname" located at site i
    IQTensor
    op(String const& opname, int i,
       Args const& args = Args::global()) const;

    void 
    read(std::istream & s) { initStream<SiteStore<>>(s); }

    void 
    write(std::ostream & s) const;

    //
    // Older interface for backwards compatibility
    //

    //Index at Site i, alternate name
    IQIndex
    si(int i) const { return operator()(i); }

    //Primed Index at Site i
    IQIndex 
    siP(int i) const { return prime(operator()(i)); }

    //Alternate name for operator()(i,state)
    IQIndexVal
    st(int i, String const& state) const
        { return operator()(i,state); }

    IQIndexVal
    stP(int i, String const& state) const
        { return prime(operator()(i,state)); }


    protected:

    template<class SiteType>
    void
    init(std::vector<SiteType> && sites);

    template<typename SiteType>
    void
    initStream(std::istream & s);

    };

class GenericSite 
    { 
    IQIndex i;
    public:

    GenericSite() { }

    GenericSite(IQIndex i_) : i(i_) { }

    IQIndex
    index() const { return i; }

    IQTensor
    op(std::string const& opname,
       Args const& args) const
        {
        Error("\'op\' method not defined for generic site");
        return IQTensor{};
        }

    IQIndexVal
    state(std::string const& state)
        {
        Error("\'state\' method not defined for generic site");
        return IQIndexVal{};
        }
    };

//Placeholder for SiteStorage below
struct NoSite
    { 
    NoSite() { }

    IQIndex
    index() const { return IQIndex{}; }

    IQTensor
    op(std::string const& opname,
       Args const& args) const
        {
        return IQTensor{};
        }

    IQIndexVal
    state(std::string const& state)
        {
        return IQIndexVal{};
        }
    };

struct BaseSiteStore
    {
    int virtual
    N() const = 0;

    IQIndex virtual
    si(int j) const = 0;

    IQIndexVal virtual
    state(int j,
          std::string const& state) = 0;

    IQTensor virtual
    op(int j,
       std::string const& opname,
       Args const& args) const = 0;
    };

template<class SiteType1 = GenericSite, class SiteType2 = GenericSite>
struct SiteStore : public BaseSiteStore
    {
    struct Element
        {
        int tag = 0;
        union s
            {
            SiteType1 s1;
            SiteType2 s2;
            };
        };
    using storage = std::vector<Element>;
    private:
    storage sites_;
    public:

    SiteStore() { }

    SiteStore(int N) : sites_(1+N) { }

    void
    set(int i, SiteType1 const& s1) 
        {
        sites_.at(i).tag = 1;
        sites_.at(i).s = s1;
        }

    void
    set(int i, SiteType2 const& s2) 
        {
        sites_.at(i).tag = 2;
        sites_.at(i).s = s2;
        }

    void
    set(int i, IQIndex const& i, int tag = 1) 
        {
        if(not (tag == 1 || tag == 2)) Error("Invalid tag passed to set");
        sites_.at(i).tag = tag;
        if(tag == 1)
            {
            sites_[i].s = SiteType1(i);
            }
        else
            {
            sites_[i].s = SiteType2(i);
            }
        }

    int
    tag(int i) const { return sites_.at(i).tag; }

    int
    N() const { return sites_.empty() ? 0 : sites_.size()-1ul; }

    IQIndex
    si(int j) const 
        { 
        auto& el = sites_.at(j);
        return el.tag==1 ? el.s1.index() : el.s2.index();
        }

    IQIndexVal
    state(int j,
          std::string const& state)
        {
        auto& el = sites_.at(j);
        return el.tag==1 ? el.s1.state(state) : el.s2.state(state);
        }

    IQTensor
    op(int j,
       std::string const& opname,
       Args const& args) const
        {
        auto& el = sites_.at(j);
        return el.tag==1 ? el.s1.op(opname,args) : el.s2.op(opname,args);
        }
    };



int inline SiteSet::
N() const { return sites_ ? sites_->N() : 0; }

inline IQIndex SiteSet::
operator()(int i) const
    {
    if(not *this) Error("Cannot retrieve site from default-initialized SiteSet");
    return sites_->si(i);
    }

IQIndexVal inline SiteSet::
operator()(int i, String const& state) const
    {
    if(not *this) Error("Cannot retrieve state from default-initialized SiteSet");
    return sites_->state(i,state);
    }

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
        auto op1 = [](std::string const& opname, size_t n)
            {
            return opname.substr(0,n);
            };
        auto op2 = [](std::string const& opname, size_t n)
            {
            return opname.substr(n+1);
            };

        //If opname of the form "Name1*Name2",
        //return product of Name1 operator times Name2 operator
        auto found = opname.find_first_of('*');
        if(found != std::string::npos)
            {
            try {
            return multSiteOps(sites_->op(i,op1(opname,found),args),
                               sites_->op(i,op2(opname,found),args));
                }
            catch(ITError const& e)
                {
                println("opname = ",opname);
                printfln("found = %s",found);
                Error("More than one * in operator name string?");
                }
            }
        return sites_->op(i,opname,args);
        }
    }

template<class StoreType>
void SiteSet::
init(StoreType && store)
    { 
    static_assert(std::is_base_of<BaseSiteStore,StoreType>::value,
                  "StoreType must be derived from BaseSiteStore");
    sites_ = std::make_shared<StoreType>(std::move(store));
    }

template<typename StoreType>
void SiteSet::
initStream(std::istream & s)
    {
    static_assert(std::is_base_of<BaseSiteStore,StoreType>::value,
                  "StoreType must be derived from BaseSiteStore");
    int N = 0;
    s.read((char*) &N,sizeof(N));
    if(N > 0)
        {
        auto store = StoreType(N);
        for(int j = 1; j <= N; ++j) 
            {
            auto I = IQIndex{};
            I.read(s);
            int tag = I.primeLevel();
            I.primeLevel(0);
            store.set(j,I,tag);
            }
        init(std::move(store));
        }
    }

void inline SiteSet::
write(std::ostream & s) const
    {
    auto N_ = N();
    s.write((char*) &N_,sizeof(N_));
    if(sites_)
        {
        for(int j = 1; j <= N_; ++j) 
            {
            auto tag = sites_->tag(j);
            auto I = sites_->si(j);
            I.primeLevel(tag);
            I.write(s);
            }
        }
    }

inline std::ostream& 
operator<<(std::ostream& s, SiteSet const& sites)
    {
    s << "SiteSet:\n";
    for(int j = 1; j <= sites.N(); ++j) 
        {
        s << format("site %d = %s\n",j,sites(j));
        }
    return s;
    }

} //namespace itensor


#endif
