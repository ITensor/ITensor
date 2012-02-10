#ifndef __ITENSOR_CONDENSER_H
#define __ITENSOR_CONDENSER_H
#include "iqtensor.h"

class Condenser	// Within one IQIndex, combine indices, presumably with same QNs
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
    doprime(PrimeType pt = primeBoth, int inc = 1)
        {
        bigind_.doprime(pt,inc);
        smallind_.doprime(pt,inc);
        }

    void product(const IQTensor& t, IQTensor& res) const;

    friend std::ostream& 
    operator<<(std::ostream & s, const Condenser & c);

private:

    IQIndex bigind_, smallind_;		// uncondensed, condensed
    mutable std::map<Index, std::pair<Index,int> > big_to_small;
    mutable std::map<std::pair<Index,int>,Index > small_to_big;

    void 
    init(const std::string& smallind_name);

}; //class Condenser


inline Condenser::
Condenser(const IQIndex& bigindex, IQIndex& smallindex)
    : bigind_(bigindex) //Use connections in bigind to create groupings
    { init(smallindex.name()); smallindex = smallind_; }

inline Condenser::
Condenser(const IQIndex& bigindex, const std::string& smallind_name)
    : bigind_(bigindex)
    { init(smallind_name); }


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
        Foreach(const inqn& x, bigind_.iq()) qns.push_back(x.qn);
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
        assert(&t != &res);
        assert(smallind_.isNotNull());
        assert(bigind_.isNotNull());
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
            GET(iqinds,smallind_pos) = bigind_;
            res = IQTensor(iqinds);
            for(IQTensor::const_iten_it i = t.const_iten_begin(); i != t.const_iten_end(); ++i)
                {
                int k;
                for(k = 1; k <= i->r(); k++)
                if(smallind_.hasindex(i->index(k))) break;

                Index sind = i->index(k);
                for(int start = 0; start < sind.m(); )
                    {
                    Index bind = small_to_big[std::make_pair(sind,start)];
                    Matrix C(sind.m(),bind.m()); C = 0;
                    for(int kk = 1; kk <= bind.m(); ++kk) { C(start+kk,kk) = 1; }
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
                Print(t); Print(*this);
                Error("Condenser::product: couldn't find bigind");
                }
            GET(iqinds,bigind_pos) = smallind_;
            res = IQTensor(iqinds);
            ITensor tt;
            Foreach(const ITensor& it, t.itensors())
                {
                bool gotit = false;
                for(int k = 1; k <= it.r(); ++k)
                if(bigind_.hasindex(it.index(k)))
                    {
                    std::pair<Index,int> Ii = big_to_small[it.index(k)];
                    tt = it;
                    tt.expandIndex(it.index(k),Ii.first,Ii.second);
                    res += tt;
                    gotit = true;
                    break;
                    }
                if(!gotit)
                    {
                    Print(*this);
                    Print(it);
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
