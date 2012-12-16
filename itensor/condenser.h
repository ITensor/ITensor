//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_CONDENSER_H
#define __ITENSOR_CONDENSER_H
#include "iqtensor.h"

//
// Condenser
//
// Within one IQIndex, combine all Index's having the same QNs
//

class Condenser
    {
    public:

    Condenser() { }

    Condenser(const IQIndex& bigindex, IQIndex& smallindex);

    Condenser(const IQIndex& bigindex, const std::string& smallind_name);

    const IQIndex& 
    bigind() const   { return bigind_; }

    const IQIndex& 
    smallind() const { return smallind_; }

    IQTensor 
    operator*(const IQTensor& t) { IQTensor res; product(t,res); return res; }

    friend inline IQTensor 
    operator*(const IQTensor& t, const Condenser& c) 
        { IQTensor res; c.product(t,res); return res; }

    void 
    conj()
        {
        bigind_.conj();
        smallind_.conj();
        }

    void 
    prime(int inc = 1) { prime(All,inc); }

    void 
    prime(IndexType type = All, int inc = 1);

    void product(const IQTensor& t, IQTensor& res) const;

    friend std::ostream& 
    operator<<(std::ostream & s, const Condenser & c);

private:

    ///////////////
    //
    // Data Members
    //

    IQIndex bigind_,   //uncondensed
            smallind_; //condensed

    mutable std::map<Index, std::pair<Index,int> > 
    big_to_small;

    mutable std::map<std::pair<Index,int>,Index > 
    small_to_big;

    //
    //////////////

    typedef std::map<Index, std::pair<Index,int> >::const_iterator
    bts_const_it;

    typedef std::map<std::pair<Index,int>,Index >::const_iterator
    stb_const_it;

    void 
    init(const std::string& smallind_name);

    }; //class Condenser

void inline Condenser::
prime(IndexType type, int inc)
    {
    bigind_.prime(type,inc);
    smallind_.prime(type,inc);

    small_to_big.clear();

    std::map<Index,std::pair<Index,int> >
    new_bts;
    for(bts_const_it it = big_to_small.begin();
        it != big_to_small.end(); ++it)
        {
        Index pindex = it->first;
        pindex.prime(type,inc);

        std::pair<Index,int> newpair = it->second;
        newpair.first.prime(type,inc);

        new_bts[pindex] = newpair;
        small_to_big[newpair] = pindex;
        }
    big_to_small.swap(new_bts);
    }


inline Condenser::
Condenser(const IQIndex& bigindex, IQIndex& smallindex)
    : bigind_(bigindex) //Use connections in bigind to create groupings
    { 
    init(smallindex.name()); 
    smallindex = smallind_; 
    }

inline Condenser::
Condenser(const IQIndex& bigindex, const std::string& smallind_name)
    : bigind_(bigindex)
    { 
    init(smallind_name); 
    }


// Use connections in t to create groupings; big = uncondensed, small = cond
inline
void Condenser::
init(const std::string& smallind_name)
    {
    /* May not be appropriate when orthogonalizing MPOs
       Could be a hint that there is a better way...
    if(bigind_.dir() != _smallind.dir())
    {
        std::cerr << "bigind_ = " << bigind_ << std::endl;
        std::cerr << "_smallind = " << _smallind << std::endl;
        Error("Arrow dirs not the same in Condenser.");
    }
    */
    static std::vector<QN> qns(1000);
    qns.resize(0);
    Foreach(const inqn& x, bigind_.iq()) 
        qns.push_back(x.qn);

    sort(qns.begin(),qns.end());

    std::vector<QN>::iterator ue = unique(qns.begin(),qns.end());

    std::vector<inqn> iq;
    for(std::vector<QN>::iterator qi = qns.begin(); qi != ue; ++qi)
        {
        const QN& q = *qi;

        int totm = 0;
        Foreach(const inqn& x, bigind_.iq())
            if(x.qn == q) totm += x.index.m();

        Index small_qind("condensed",totm);
        int start = 0;
        Foreach(const inqn& x, bigind_.iq())
            if(x.qn == q)
                {
                const Index &xi = x.index;
                small_to_big[std::make_pair(small_qind,start)] = xi;
                big_to_small[xi] = std::make_pair(small_qind,start);
                start += xi.m();
                }
        iq.push_back(inqn(small_qind,q));
        }

    smallind_ = IQIndex(smallind_name,iq,bigind_.dir(),bigind_.primeLevel());

    bigind_.conj();
    }

inline
void Condenser::
product(const IQTensor& t, IQTensor& res) const
    {
    if(&t == &res)
        Error("Cannot condense into same IQTensor");

    std::vector<IQIndex> iqinds; iqinds.reserve(t.r());

    int smallind_pos = -2;
    int bigind_pos   = -2;
    for(int j = 1; j <= t.r(); ++j)
        {
        iqinds.push_back(t.index(j));

        if(iqinds.back() == smallind_) 
            {
            if(iqinds.back().dir() == smallind_.dir())
                {
                Print(smallind_);
                Error("Incompatible Arrow for smallind");
                }
            smallind_pos = (j-1);
            }
        else if(iqinds.back() == bigind_) 
            {
            if(iqinds.back().dir() == bigind_.dir())
                {
                Print(bigind_);
                Error("Incompatible Arrow for bigind");
                }
            bigind_pos = (j-1);
            }
        }

    if(smallind_pos != -2) //expand condensed form into uncondensed
        {
        iqinds.at(smallind_pos) = bigind_;

        res = IQTensor(iqinds);

        for(IQTensor::const_iten_it i = t.const_iten_begin(); i != t.const_iten_end(); ++i)
            {
            int k;
            for(k = 1; k <= i->r(); ++k)
                if(smallind_.hasindex(i->index(k))) break;

            Index sind = i->index(k);
            for(int start = 0; start < sind.m(); )
                {
                Index bind = small_to_big[std::make_pair(sind,start)];
                Matrix C(sind.m(),bind.m()); C = 0;
                for(int kk = 1; kk <= bind.m(); ++kk) 
                    { 
                    C(start+kk,kk) = 1; 
                    }
                ITensor converter(sind,bind,C);
                converter *= (*i);
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

        Foreach(ITensor tt, t.blocks())
            {
            bool gotit = false;

            for(int k = 1; k <= tt.r(); ++k)
                if(bigind_.hasindex(tt.index(k)))
                    {
                    std::pair<Index,int> Ii = big_to_small[tt.index(k)];
                    tt.expandIndex(tt.index(k),Ii.first,Ii.second);
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

inline
std::ostream& 
operator<<(std::ostream & s, const Condenser & c)
    {
    s << "bigind_ is " << c.bigind_ << "\n";
    s << "smallind_ is " << c.smallind_ << "\n";
    s << "big_to_small is " << "\n";
    for(std::map<Index, std::pair<Index,int> >::const_iterator kk = c.big_to_small.begin();
        kk != c.big_to_small.end(); ++kk)
        { s << kk->first SP kk->second.first SP kk->second.second << "\n"; }
    s << "small_to_big is " << std::endl;
    for(std::map<std::pair<Index,int>,Index>::const_iterator kk = c.small_to_big.begin();
        kk != c.small_to_big.end(); ++kk)
        { s << kk->first.first SP kk->first.second SP kk->second << "\n"; }
    return s << std::endl;
    }

#endif
