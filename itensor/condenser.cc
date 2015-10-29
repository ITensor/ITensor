//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#include "itensor/condenser.h"
#include <map>
#include <algorithm>

namespace itensor {

using std::ostream;
using std::vector;
using std::map;


Condenser::
Condenser(const IQIndex& bigindex, IQIndex& smallindex)
    : bigind_(bigindex) //Use connections in bigind to create groupings
    { 
    init(smallindex.name()); 
    smallindex = smallind_; 
    }

Condenser::
Condenser(const IQIndex& bigindex, const std::string& smallind_name)
    : bigind_(bigindex)
    { 
    init(smallind_name); 
    }

void Condenser::
dag()
    {
    bigind_.dag();
    smallind_.dag();
    }


// Use connections in t to create groupings; big = uncondensed, small = cond

void Condenser::
init(const std::string& smallind_name)
    {
    static std::vector<QN> qns(10);
    qns.resize(0);
    for(const IndexQN& x : bigind_.indices()) 
        qns.push_back(x.qn);

    sort(qns.begin(),qns.end());

    std::vector<QN>::iterator ue = unique(qns.begin(),qns.end());

    IQIndex::Storage iq;
    for(std::vector<QN>::iterator qi = qns.begin(); qi != ue; ++qi)
        {
        const QN& q = *qi;

        int totm = 0;
        for(const IndexQN& x : bigind_.indices())
            if(x.qn == q) totm += x.m();

        Index small_qind(smallind_name,totm,bigind_.type(),bigind_.primeLevel());
        int start = 0;
        for(const IndexQN& x : bigind_.indices())
            if(x.qn == q)
                {
                maps_.push_back(IndexMap(small_qind,start,x));
                start += x.m();
                }
        iq.push_back(IndexQN(small_qind,q));
        }

    smallind_ = IQIndex(smallind_name,iq,bigind_.dir(),bigind_.primeLevel());

    bigind_.dag();
    }

void Condenser::
prime(IndexType type, int inc)
    {
    bigind_.prime(type,inc);
    smallind_.prime(type,inc);

    for(IndexMap& m : maps_)
        {
        m.big.prime(type,inc);
        m.small.prime(type,inc);
        }
    }

static const Index&
findBig(const Index& small, int j,
        const vector<IndexMap>& maps)
    {
    for(const IndexMap& m : maps)
        {
        if(m.i == j && m.small == small)
            return m.big;
        }
    Error("small,j pair not found");
    return maps.front().big;
    }

static const IndexMap&
findSmall(const Index& big,
        const vector<IndexMap>& maps)
    {
    for(const IndexMap& m : maps)
        {
        if(m.big == big)
            return m;
        }
    Error("'big' Index not found");
    return maps.front();
    }

void Condenser::
product(const IQTensor& t, IQTensor& res) const
    {
    if(&t == &res)
        Error("Cannot condense into same IQTensor");

    std::vector<IQIndex> iqinds; 
    iqinds.reserve(t.r());

    int smallind_pos = -2;
    int bigind_pos   = -2;
    int j = 0;
    for(const IQIndex& J : t.indices())
        {
        iqinds.push_back(J);

        if(iqinds.back() == smallind_) 
            {
            if(iqinds.back().dir() == smallind_.dir())
                {
                Print(smallind_);
                Error("Incompatible Arrow for smallind");
                }
            smallind_pos = j;
            }
        else if(iqinds.back() == bigind_) 
            {
            if(iqinds.back().dir() == bigind_.dir())
                {
                Print(bigind_);
                Error("Incompatible Arrow for bigind");
                }
            bigind_pos = j;
            }
        ++j;
        }

    if(smallind_pos != -2) //expand condensed form into uncondensed
        {
        iqinds.at(smallind_pos) = bigind_;

        res = IQTensor(iqinds);

        for(const ITensor& b : t.blocks())
            {
            Index sind;
            for(const Index& I : b.indices())
                if(hasindex(smallind_,I))
                    {
                    sind = I;
                    break;
                    }

            for(int start = 0; start < sind.m(); )
                {
                const Index& bind = findBig(sind,start,maps_);
                Matrix C(sind.m(),bind.m()); 
                C = 0;
                for(int kk = 1; kk <= bind.m(); ++kk) 
                    { 
                    C(start+kk,kk) = 1; 
                    }
                ITensor converter(sind,bind,C);
                converter *= b;
                res += converter;
                start += bind.m();
                }
            }
        }
    else //contract regular form into condensed
        {
        if(bigind_pos == -2)
            {
            Print(t); 
            Print(*this);
            Error("Condenser::product: couldn't find bigind");
            }
        iqinds.at(bigind_pos) = smallind_;

        res = IQTensor(iqinds);

        for(ITensor tt : t.blocks())
            {
            bool gotit = false;

            for(const Index& K : tt.indices())
                if(hasindex(bigind_,K))
                    {
                    const IndexMap& m = findSmall(K,maps_);
                    tt.expandIndex(K,m.small,m.i);
                    res += tt;
                    gotit = true;
                    break;
                    }

            if(!gotit)
                {
                Print(*this);
                Print(tt);
                Error("Combiner::product: Can't find common Index");
                }
            }
        }
    }

std::ostream& 
operator<<(std::ostream & s, const Condenser & c)
    {
    s << "bigind_ is " << c.bigind() << "\n";
    s << "smallind_ is " << c.smallind() << "\n";
    s << "index maps \n";
    for(const IndexMap& m : c.maps())
        {
        s << "(" << m.small << "," << m.i << ") " << m.big << "\n";
        }
    return s << std::endl;
    }

} //namespace itensor

