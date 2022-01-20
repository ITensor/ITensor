//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef __ITENSOR_SITESET_H
#define __ITENSOR_SITESET_H
#include "itensor/itensor.h"
#include "itensor/util/str.h"

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

class SiteBase;
class GenericSite;

class SiteSet
    {
 
    private:
    typedef std::shared_ptr<SiteBase> SitePtr; //use shared pointer to get copy on write semantics.
    typedef std::vector    <SitePtr>  SitePtrs; //Should be safe to copy.
    SitePtrs pm_sites_; //Polymorphic array of sites;

    public:

    using String = std::string;

    SiteSet() { }

    //Create generic SiteSet of N sites with local dimension d
    SiteSet(int N, int d);

    //Create SiteSet from vector of provided Indices (0-indexed)
    SiteSet(IndexSet const& is);

    //Create SiteSet of length N (can be set using .set method)
    SiteSet(int N);

    explicit operator bool() const { return pm_sites_.size()>1; } //compensate for dummy elements at index 0.

    int 
    length() const;

		// Deprecated in favor of .length()
    int 
    N() const;

    //Index at Site i
    Index
    operator()(int i) const;

    //Return an IndexSet of the indices of the SiteSet
    IndexSet
    inds() const;

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


    template<typename SiteType>
    void
    set(int i, SiteType && s);

    protected:

    void 
    allocate(int N);
    
    void
    insert(int j, SiteBase*);

    template<typename SiteType>
    void
    readType(std::istream & s);

    }; //class SiteSet

//----------------------------------------------------------------------
//
//  Abstract SiteBase is a virtual base class for concrete site types.
//  This allows us to store a polymorphic array similar to std::vector<SiteBase*>
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

class GenericSite : public virtual SiteBase
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

//---------------------------------------------------------------------------
//
//  SiteSet inline implementations.
//
inline SiteSet::
SiteSet(int N, int d)
    : pm_sites_(N+1) // 1 based vector with unused [0] element
    {
    for(int j = 1; j <= N; ++j)
        {
        auto I = Index(d,"Site,n="+str(j));
        insert(j,new GenericSite(I));
        }
    }

inline SiteSet::
SiteSet(int N)
    : pm_sites_(N+1) // 1 based vector with unused [0] element
    {
    // We don;t know d here so we can't create any sites.
    }

inline SiteSet::
SiteSet(IndexSet const& is)
    : pm_sites_(is.length()+1) // 1 based vector with unused [0] element
    {
    auto N = is.length();
    for(auto j : range1(N))
        insert(j,new GenericSite(is(j)));
    }


int inline SiteSet::
length() const { return operator bool() ? pm_sites_.size()-1 : 0; } //1 based vector with dummy element at index 0

int inline
length(SiteSet const& sites) { return sites.length(); }

IndexSet inline
inds(SiteSet const& sites) { return sites.inds(); }

// Deprecated in favor of .length()
int inline SiteSet::
N() const 
    {
    Global::warnDeprecated(".N() is deprecated in favor of length(SiteSet)");
    return this->length();
    }

Index inline SiteSet::
operator()(int i) const
    {
    if(not *this) Error("Cannot retrieve site from default-initialized SiteSet");
    return pm_sites_[i]->index();
    }

IndexSet inline SiteSet::
inds() const
    {
    if(not *this) Error("Cannot retrieve sites from default-initialized SiteSet");
    auto N = length();
    auto is = IndexSetBuilder(N);
    for( auto n : range1(N) )
        {
        auto I = (*this)(n); //Reuse op() function
        is.nextIndex(std::move(I));
        }
    return is.build();
    }

IndexVal inline SiteSet::
operator()(int i, String const& state) const
    {
    if(not *this) Error("Cannot retrieve state from default-initialized SiteSet");
    return pm_sites_[i]->state(state);
    }

bool inline
hasQNs(SiteSet const& sites)
    {
    for(auto i : range1(length(sites)))
        if(not hasQNs(sites(i))) return false;
    return true;
    }

ITensor inline SiteSet::
op(String const& opname, 
   int i, 
   Args const& args) const
    { 
    if(not *this) Error("Cannot call .op(..) on default-initialized SiteSet");
    if(opname == "Id")
        {
        auto s = si(i);
        auto id_ = ITensor(dag(s),prime(s));
        for(auto j : range1(dim(s))) id_.set(j,j,1.0);
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
    if(opname == "F")
        {
        try {
            return pm_sites_[i]->op(opname,args);
        } catch(...) {
            //
            // If no "F" operator defined by site set
            // then return a 'trivial' F operator
            // equal to the identity:
            auto s = si(i);
            auto trivF = ITensor(dag(s),prime(s));
            for(auto j : range1(dim(s))) trivF.set(j,j,1.0);
            return trivF;
        }
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
        return pm_sites_[i]->op(opname,args);
        }
    }

// Toss out any old sites and replace with an array of N empty SitePtr .
void inline SiteSet::allocate(int N)
{
    assert(N>1);
    pm_sites_=SitePtrs(N+1); //Fake one based indexing.
}
        
void inline SiteSet::
insert(int j, SiteBase* s)
{
    assert(s);
    assert(j>0);
    assert(static_cast<SitePtrs::size_type>(j)<=pm_sites_.size());
    pm_sites_[j]=SitePtr(s);
}
    
    
template<typename SiteType>
void SiteSet::
readType(std::istream & s)
    {
    int N = itensor::read<int>(s);
    allocate(N); //Allocate array of site pointers
    if(N > 0)
        {
        for(int j = 1; j <= N; ++j) 
            {
            auto I = Index{};
            I.read(s);
            insert(j,new SiteType(I)); //JRTODO support factory for unpickling polymorphic sites.
            }
        }
    }


void inline SiteSet::
write(std::ostream & s) const
    {
    itensor::write(s,length());
    if(pm_sites_.size())
        {
        for(int j = 1; j <= length(); ++j) 
            (*this)(j).write(s);
        }
    }

//---------------------------------------------------------------------------
//
//  SiteSet inline friend functions.
//
ITensor inline
op(SiteSet const& sites,
   std::string const& opname,
   int i,
   Args const& args = Args::global())
    {
    return sites.op(opname,i,args);
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

//---------------------------------------------------------------------------
//
//  Templated base class for concrete SiteSet classes.
//
template<typename SiteType>
class BasicSiteSet : public SiteSet
    {
    public:

    BasicSiteSet() { }

    BasicSiteSet(int N, 
                 Args const& args = Args::global())
         : SiteSet(N)
        {
        for(int j = 1; j <= N; ++j)
            insert (j,new SiteType({args,"SiteNumber=",j})); 
        }

    BasicSiteSet(IndexSet const& is)
         : SiteSet(is.length())
        {
        int N = is.length();
        for(auto j : range1(N))
            insert (j,new SiteType(is(j)));
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
        : SiteSet(N)
        {
        for(int j = 1; j <= N; ++j)
            {
            if(j%2 == 1) insert(j,new ASiteType({args,"SiteNumber=",j}));
            else         insert(j,new BSiteType({args,"SiteNumber=",j}));
            }
        }

    MixedSiteSet(IndexSet const& is)
        : SiteSet(is.length())
        {
        int N = is.length();
        for(auto j : range1(N))
            {
            if(j%2 == 1) insert(j,new ASiteType(is(j)));
            else         insert(j,new BSiteType(is(j)));
            }
        }

    void
    read(std::istream& s)
        {
        int N = itensor::read<int>(s);
        allocate(N); //Allocate array of site pointers
        if(N > 0)
            {
            for(int j = 1; j <= N; ++j) 
                {
                auto I = Index{};
                I.read(s);
                if(j%2==1) insert(j,new ASiteType(I));
                else       insert(j,new BSiteType(I));
                }
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

//
// For site types that define member
// functions .index(), .state(), and .op(),
// these are the free function versions
//

template <typename SiteType>
Index
index(SiteType const& s)
    {
    return s.index();
    }

template <typename SiteType>
IndexVal
state(SiteType const& s,
      std::string const& st)
    {
    return s.state(st);
    }

template <typename SiteType>
ITensor
op(SiteType const& s,
   std::string const& opname,
   Args const& args = Args::global())
    {
    return s.op(opname,args);
    }


} //namespace itensor


#endif
