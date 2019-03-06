////
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SITESET_H
#define __ITENSOR_SITESET_H
#include "itensor/itensor.h"

namespace itensor {

//
// Classes derived from SiteSet 
// represent the Hilbert space of a 
// system as a set of Site indices.
//
// The convention for operators is
// that they are 2-index ITensors
// with the Site Index pointing
// In and the Site' Index pointing
// Out. This is so we can compute expectation
// values by doing dag(prime(A,Site)) * Op * A.
// (assuming the tensor A is an ortho center 
// of our MPS)
//

class GenericSite;
struct SiteStore;

class SiteSet
    {
    std::shared_ptr<SiteStore> sites_;
    public:

    using String = std::string;

    SiteSet() { }

    //Create generic SiteSet of N sites with local dimension d
    SiteSet(int N, int d);

    //Create SiteSet from vector of provided Indices (0-indexed)
    SiteSet(std::vector<Index> const& inds);

    explicit operator bool() const { return bool(sites_); }

    int 
    length() const;

		// Deprecated in favor of .length()
    int 
    N() const;

    //Index at Site i
    Index
    operator()(int i) const;

    //Index at site i set to a certain state
    //indicated by the string "state"
    //e.g. sites(5,"Up") returns the IndexVal
    //representing the spin up state on site 5
    //(assuming a spin SiteSet such as SpinHalf)
    IndexVal
    operator()(int i, String const& state) const;

    //Get the operator indicated by
    //"opname" located at site i
    ITensor
    op(String const& opname, int i,
       Args const& args = Args::global()) const;

    void 
    read(std::istream & s) { readType<GenericSite>(s); }

    void 
    write(std::ostream & s) const;

    //
    // Older interface for backwards compatibility
    //

    //Index at Site i, alternate name
    Index
    si(int i) const { return operator()(i); }

    //Primed Index at Site i
    Index 
    siP(int i) const { return prime(operator()(i)); }

    //Alternate name for operator()(i,state)
    IndexVal
    st(int i, String const& state) const
        { return operator()(i,state); }

    IndexVal
    stP(int i, String const& state) const
        { return prime(operator()(i,state)); }


    protected:

    void
    init(SiteStore && sites);

    template<typename SiteType>
    void
    readType(std::istream & s);

    };

//
// "Base" type for virtual mechanism
//
class SiteBase 
    { 
    public:

    SiteBase() { }

    virtual ~SiteBase() { }

    Index virtual
    index() const = 0;

    ITensor virtual
    op(std::string const& opname,
       Args const& args) const = 0;

    IndexVal virtual
    state(std::string const& state) = 0;
    };

//
// Derived "box" type with virtual methods
// Wraps any object implementing the "SiteType" interface
//
template<typename SiteType>
class SiteHolder : public SiteBase
    { 
    SiteType s;
    public:

    SiteHolder() { }

    SiteHolder(SiteType && s_) : s(std::move(s_)) { }

    virtual ~SiteHolder() { }

    Index virtual
    index() const { return s.index(); }

    ITensor virtual
    op(std::string const& opname,
       Args const& args) const
        {
        return s.op(opname,args);
        }

    IndexVal virtual
    state(std::string const& state)
        {
        return s.state(state);
        }
    };

class GenericSite
    { 
    Index i;
    public:

    GenericSite() { }

    GenericSite(Index i_) : i(i_) { }

    Index
    index() const { return i; }

    ITensor
    op(std::string const& opname,
       Args const& args) const
        {
        Error("\'op\' method not defined for generic site");
        return ITensor{};
        }

    IndexVal
    state(std::string const& state)
        {
        Error("\'state\' method not defined for generic site");
        return IndexVal{};
        }
    };


struct SiteStore
    {
    using sptr = std::unique_ptr<SiteBase>;
    using storage = std::vector<sptr>;
    private:
    storage sites_;
    public:

    SiteStore() { }

    SiteStore(int N) : sites_(1+N) { }

    template<typename SiteType>
    void
    set(int i, SiteType && s) 
        {
        sites_.at(i) = sptr(new SiteHolder<SiteType>(std::move(s)));
        }

    int
    length() const { return sites_.empty() ? 0 : sites_.size()-1ul; }

		// Deprecated in favor of .length()
    int
    N() const 
        {
				Global::warnDeprecated(".N() is deprecated in favor of length(SiteStore)");
        return this->length();
        }

    Index
    si(int j) const 
        { 
        if(not sites_.at(j)) Error("Unassigned site in SiteStore");
        return sites_[j]->index();
        }

    IndexVal
    state(int j,
          std::string const& state)
        {
        if(not sites_.at(j)) Error("Unassigned site in SiteStore");
        return sites_[j]->state(state);
        }

    ITensor
    op(int j,
       std::string const& opname,
       Args const& args) const
        {
        if(not sites_.at(j)) Error("Unassigned site in SiteStore");
        return sites_[j]->op(opname,args);
        }
    };


inline SiteSet::
SiteSet(int N, int d)
    {
    auto sites = SiteStore(N);
    for(int j = 1; j <= N; ++j)
        {
        auto I = Index(d,format("Site,n=%d",j));
        sites.set(j,GenericSite(I));
        }
    SiteSet::init(std::move(sites));
    }

inline SiteSet::
SiteSet(std::vector<Index> const& inds)
    {
    auto sites = SiteStore(inds.size());
    for(auto j : range(inds))
        {
        sites.set(1+j,GenericSite(inds.at(j)));
        }
    SiteSet::init(std::move(sites));
    }


int inline SiteSet::
length() const { return sites_ ? sites_->length() : 0; }

int inline
length(SiteSet const& sites) { return sites.length(); }

// Deprecated in favor of .length()
int inline SiteSet::
N() const 
    {
    Global::warnDeprecated(".N() is deprecated in favor of length(SiteSet)");
    return this->length();
    }

inline Index SiteSet::
operator()(int i) const
    {
    if(not *this) Error("Cannot retrieve site from default-initialized SiteSet");
    return sites_->si(i);
    }

IndexVal inline SiteSet::
operator()(int i, String const& state) const
    {
    if(not *this) Error("Cannot retrieve state from default-initialized SiteSet");
    return sites_->state(i,state);
    }

bool inline
hasQNs(SiteSet const& sites)
    {
    for(auto i : range1(length(sites)))
        if(not hasQNs(sites(i))) return false;
    return true;
    }

inline ITensor SiteSet::
op(String const& opname, 
   int i, 
   Args const& args) const
    { 
    if(not *this) Error("Cannot call .op(..) on default-initialized SiteSet");
    if(opname == "Id")
        {
        Index s = dag(si(i));
        Index sP = siP(i);
        auto id_ = ITensor(s,sP);
        for(int j = 1; j <= dim(s); ++j)
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
            return multSiteOps(op(op1(opname,found),i,args),
                               op(op2(opname,found),i,args));
            }
        return sites_->op(i,opname,args);
        }
    }

void inline SiteSet::
init(SiteStore && store)
    { 
    sites_ = std::make_shared<SiteStore>(std::move(store));
    }

template<typename SiteType>
void SiteSet::
readType(std::istream & s)
    {
    int N = itensor::read<int>(s);
    if(N > 0)
        {
        auto store = SiteStore(N);
        for(int j = 1; j <= N; ++j) 
            {
            auto I = Index{};
            I.read(s);
            store.set(j,SiteType(I));
            }
        init(std::move(store));
        }
    }

void inline SiteSet::
write(std::ostream & s) const
    {
    itensor::write(s,length());
    if(sites_)
        {
        for(int j = 1; j <= length(); ++j) 
            {
            sites_->si(j).write(s);
            }
        }
    }

inline std::ostream& 
operator<<(std::ostream& s, SiteSet const& sites)
    {
    s << "SiteSet:\n";
    for(int j = 1; j <= length(sites); ++j) 
        {
        s << format("site %d = %s\n",j,sites(j));
        }
    return s;
    }

template<typename SiteType>
class BasicSiteSet : public SiteSet
    {
    public:

    BasicSiteSet() { }

    BasicSiteSet(int N, 
                 Args const& args = Args::global())
        {
        auto sites = SiteStore(N);
        for(int j = 1; j <= N; ++j)
            {
            sites.set(j,SiteType(j,args));
            }
        SiteSet::init(std::move(sites));
        }

    BasicSiteSet(std::vector<Index> const& inds)
        {
        int N = inds.size();
        auto sites = SiteStore(N);
        for(int j = 1, i = 0; j <= N; ++j, ++i)
            {
            auto& Ii = inds.at(i);
            sites.set(j,SiteType(Ii));
            }
        SiteSet::init(std::move(sites));
        }

    void
    read(std::istream& s)
        {
        SiteSet::readType<SiteType>(s);
        }

    };

template<typename ASiteType, typename BSiteType>
class MixedSiteSet : public SiteSet
    {
    public:

    MixedSiteSet() { }

    MixedSiteSet(int N, 
                 Args const& args = Args::global())
        {
        auto sites = SiteStore(N);
        for(int j = 1; j <= N; ++j)
            {
            if(j%2 == 1) sites.set(j,ASiteType(j,args));
            else         sites.set(j,BSiteType(j,args));
            }
        SiteSet::init(std::move(sites));
        }

    MixedSiteSet(std::vector<Index> const& inds)
        {
        int N = inds.size();
        auto sites = SiteStore(N);
        for(int j = 1, i = 0; j <= N; ++j, ++i)
            {
            auto& Ii = inds.at(i);
            if(j%2 == 1) sites.set(j,ASiteType(Ii));
            else         sites.set(j,BSiteType(Ii));
            }
        SiteSet::init(std::move(sites));
        }

    void
    read(std::istream& s)
        {
        int N = itensor::read<int>(s);
        if(N > 0)
            {
            auto store = SiteStore(N);
            for(int j = 1; j <= N; ++j) 
                {
                auto I = Index{};
                I.read(s);
                if(j%2==1) store.set(j,ASiteType(I));
                else       store.set(j,BSiteType(I));
                }
            init(std::move(store));
            }
        }
    };

template<typename SiteSetT,
         class = stdx::enable_if_t<std::is_base_of<SiteSet,SiteSetT>::value>>
SiteSetT
merge(SiteSetT const& sites1,
      SiteSetT const& sites2,
      Args const& args = Args::global())
    {
    auto N1 = length(sites1);
    auto N2 = length(sites2);
    auto start1 = args.getInt("Start1",1);
    auto end1 = args.getInt("End1",N1);
    auto start2 = args.getInt("Start2",1);
    auto end2 = args.getInt("End2",N2);
    if(end1 <= start1)
        {
        println("Start1=",start1);
        println("End1=",end1);
        Error("end1 must be greater than start1");
        }
    if(end2 <= start2)
        {
        println("Start2=",start2);
        println("End2=",end2);
        Error("end2 must be greater than start2");
        }
    auto nsize = (end1-start1+1)+(end2-start2+1);
    auto inds = std::vector<Index>(nsize);
    auto i = 0;
    for(auto j1 : range1(start1,end1))
        {
        inds.at(i) = sites1(j1);
        ++i;
        }
    for(auto j2 : range1(start2,end2))
        {
        inds.at(i) = sites2(j2);
        ++i;
        }
    return SiteSetT{inds};
    }

} //namespace itensor


#endif
