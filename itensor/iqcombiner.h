#ifndef __IQCOMBINER_H
#define __IQCOMBINER_H
#include "combiner.h"
#include "iqtensor.h"

class Condenser	// Within one IQIndex, combine indices, presumably with same QNs
{
    IQIndex bigind_, smallind_;		// uncondensed, condensed
    mutable std::map<Index, std::pair<Index,int> > big_to_small;
    mutable std::map<std::pair<Index,int>,Index > small_to_big;

    // Use connections in t to create groupings; big = uncondensed, small = cond
    void init(const IQIndex& _bigind, const std::string& smallind_name)
    {
        /* May not be appropriate when orthogonalizing MPOs
           Could be a hint that there is a better way...
        if(_bigind.dir() != _smallind.dir())
        {
            std::cerr << "_bigind = " << _bigind << std::endl;
            std::cerr << "_smallind = " << _smallind << std::endl;
            Error("Arrow dirs not the same in Condenser.");
        }
        */
        static std::vector<QN> qns(1000);
        qns.resize(0);
        foreach(const inqn& x, bigind_.iq()) qns.push_back(x.qn);
        sort(qns.begin(),qns.end());
        std::vector<QN>::iterator ue = unique(qns.begin(),qns.end());

        std::vector<inqn> iq;
        for(std::vector<QN>::iterator qi = qns.begin(); qi != ue; ++qi)
        {
            const QN& q = *qi;
            int totm = 0;
            foreach(const inqn& x, bigind_.iq())
            if(x.qn == q) totm += x.index.m();

            Index small_qind("condensed",totm);
            int start = 0;
            foreach(const inqn& x, bigind_.iq())
            if(x.qn == q)
            {
                const Index &xi = x.index;
                small_to_big[std::make_pair(small_qind,start)] = xi;
                big_to_small[xi] = std::make_pair(small_qind,start);
                start += xi.m();
            }
            iq.push_back(inqn(small_qind,q));
        }
        smallind_ = IQIndex(smallind_name,iq,_bigind.dir(),_bigind.primeLevel());
    }
public:
    const IQIndex& bigind() const   { return bigind_; }
    const IQIndex& smallind() const { return smallind_; }
    
    Condenser() { }

    Condenser(const IQIndex& _bigind, IQIndex& _smallind)
        : bigind_(_bigind) //Use connections in bigind to create groupings
    { init(_bigind,_smallind.name()); _smallind = smallind_; }

    Condenser(const IQIndex& _bigind, const std::string& smallind_name)
    : bigind_(_bigind)
    { init(_bigind,smallind_name); }

    IQTensor operator*(const IQTensor& t) { IQTensor res; product(t,res); return res; }
    friend inline IQTensor operator*(const IQTensor& t, const Condenser& c) { IQTensor res; c.product(t,res); return res; }

    void product(const IQTensor& t, IQTensor& res) const
    {
        assert(&t != &res);
        assert(smallind_.is_not_null());
        assert(bigind_.is_not_null());
        std::vector<IQIndex> iqinds; iqinds.reserve(t.r());
        int smallind_pos = -2;
        int bigind_pos   = -2;
        for(int j = 1; j <= t.r(); ++j)
        {
            iqinds.push_back(t.index(j));
            if(iqinds.back() == smallind_) smallind_pos = (j-1);
            else if(iqinds.back() == bigind_) bigind_pos = (j-1);
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
            foreach(const ITensor& it, t.itensors())
            {
                bool gotit = false;
                for(int k = 1; k <= it.r(); ++k)
                if(bigind_.hasindex(it.index(k)))
                {
                    std::pair<Index,int> Ii = big_to_small[it.index(k)];
                    //doconvert(it,pp.first,it.index(k),pp.second,tt);
                    it.expandIndex(it.index(k),Ii.first,Ii.second,tt);
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

    inline friend std::ostream& operator<<(std::ostream & s, const Condenser & c)
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
}; //class Condenser

namespace { //Anonymous namespace means don't expose outside of this header

class QCounter
{
public:
    std::vector<int> n;
    std::vector<int> ind;
    bool don;
    QCounter(const std::vector<IQIndex>& v)
	{
        foreach(const IQIndex& I,v)
        {
            n.push_back(I.nindex());
            ind.push_back(0);
        }
        don = false;
	}
    bool notdone() const { return !don; }
    QCounter& operator++()
	{
        int nn = n.size();
        ind[0]++;
        if(ind[0] >= n[0])
        {
            for(int j = 1; j < nn; j++)
            {
                ind[j-1] = 0;
                ++ind[j];
                if(ind[j] < n[j]) break;
            }
        }
        if(ind[nn-1] >= n[nn-1])
        {
            ind = std::vector<int>(nn,0);
            don = true;
        }

        return *this;
	}

    void getVecInd(const std::vector<IQIndex>& v, std::vector<Index>& vind, QN& q) const
	{
        q = QN(); vind.clear();
        for(unsigned int i = 0; i < ind.size(); ++i)
        {
            const int j = ind[i]+1;
            if(GET(v,i).nindex() < j)
            {
                for(unsigned int k = 0; k < n.size(); ++k) std::cerr << boost::format("n[%d] = %d\n")%k%n[k];
                std::cout << boost::format("i=%d, j=%d, v[i].nindex()=%d\n")%i%j%v[i].nindex();
                Error("bad v[i].iq in getVecInd");
            }
            vind.push_back(v[i].index(j));
            q += v[i].qn(j)*v[i].dir();
        }
	} //void QCounter::getVecInd
};

}


/*
   Combine several indices into one index without loss of states.
   If the IQCombiner C is created from indices of IQTensor T, then
   an identity is   T * C * conj(C).  This looks like (for example)

        ___ __          __
       /      \        /
   ---T---- ---C-- --cC---
   where cC is conj(C).  Use of IQCombiners is efficient, whereas
   use of IQTensors for this purpose would not be.
*/
class IQCombiner
{
    std::vector<IQIndex> left;
    mutable IQIndex _right;
    mutable std::map<ApproxReal, Combiner> setcomb;
    mutable std::map<Index, Combiner> rightcomb;
    mutable bool initted;

    mutable Condenser cond;
    mutable IQIndex cindex;
    bool do_condense;
public:

    void doCondense(bool val) 
    {
        if(initted) 
            Error("IQCombiner: can't set doCondense after already initted.");
        do_condense = val;
    }

    IQCombiner() : initted(false), do_condense(false) { }
    IQCombiner(
	    const IQIndex& l1, const IQIndex& l2 = IQIndex::Null(), 
        const IQIndex& l3 = IQIndex::Null(), const IQIndex& l4 = IQIndex::Null(), 
	    const IQIndex& l5 = IQIndex::Null(), const IQIndex& l6 = IQIndex::Null() )
        : initted(false), do_condense(false)
	{
        if(l1 != IQIndex::Null()) left.push_back(l1); 
        if(l2 != IQIndex::Null()) left.push_back(l2);
        if(l3 != IQIndex::Null()) left.push_back(l3); 
        if(l4 != IQIndex::Null()) left.push_back(l4);
        if(l5 != IQIndex::Null()) left.push_back(l5); 
        if(l6 != IQIndex::Null()) left.push_back(l6);
        foreach(IQIndex& L, left) L.conj();
	}
    void addleft(const IQIndex& l) 	// Include another left index
	{ 
        left.push_back(l);
        //Flip arrows to make combiner compatible with
        //the IQTensor from which it got its left indices
        left.back().conj();
        //assert(l.dir() != left.front().dir());
        initted = false;
	}

    inline bool isInit() const { return initted; }

    // Initialize after all lefts are there and before being used
    void init(std::string rname = "combined", IndexType = Link, 
              int primelevel = -1) const 
	{
        if(initted) return;
        if(left.size() == 0)
            Error("No left indices in IQCombiner.");

        Arrow rdir = Switch*left.back().dir();
        int plev = 0;
        if(primelevel == -1) { plev = left.back().primeLevel(); }
        else                 { plev = primelevel; }

        //Prefer to derive right Arrow, primelevel from Link indices
        for(size_t j = 0; j < left.size(); ++j)
        if(left[j].type() == Link) 
        { 
            rdir = Switch*left[j].dir(); 
            if(primelevel == -1) { plev = left[j].primeLevel(); }
            break;
        }

        setcomb.clear();
        rightcomb.clear();

        //Construct individual Combiners
        QCounter c(left);
        std::vector<inqn> iq;
        for( ; c.notdone(); ++c)
        {
            std::vector<Index> vind;
            QN q;
            c.getVecInd(left, vind, q);		// updates vind and q
            q *= -rdir;

            Combiner co; Real rss = 0.0;
            foreach(const Index& i, vind)
            { 
                co.addleft(i); 
                rss += i.unique_Real(); 
            }
            co.init(rname+q.toString());

            iq.push_back(inqn(co.right(),q));
            setcomb[ApproxReal(rss)] = co;
            rightcomb[co.right()] = co;
        }
        if(do_condense) 
        {
            cindex = IQIndex(rname,iq,rdir,plev);
            std::string cname = "cond::" + rname;
            cond = Condenser(cindex,cname);
            _right = cond.smallind();
        }
        else _right = IQIndex(rname,iq,rdir,plev);

        initted = true;
	}
    
    operator IQTensor() const
    {
        if(!initted) Error("IQCombiner::operator IQTensor(): IQCombiner not initialized.");

        //if(_right.m() > 16) 
        //{ std::cerr << std::endl << std::endl << "WARNING: too large of an m in IQCombiner::operator IQTensor(). May be inefficient!" << std::endl << std::endl; }

        std::vector<IQIndex> iqinds(left);
        iqinds.push_back((do_condense ? cindex : _right));
        IQTensor res(iqinds);
        for(std::map<ApproxReal,Combiner>::const_iterator it = setcomb.begin();
            it != setcomb.end(); 
            ++it)
            { res.insert(it->second); }

        //Combiners should always have the 
        //structure of zero divergence IQTensors
        assert(checkQNs(res));

        if(do_condense) { IQTensor rcopy(res); cond.product(rcopy,res); }

        return res;
    }

    const IQIndex& right() const 
    { 
        if(!initted) Error("IQCombiner::right(): IQCombiner not initialized.");
        return _right; 
    }

    int findindex(const IQIndex& i) const
	{
        for(size_t j = 0; j < left.size(); ++j)
        if(left[j] == i) return j;
        return -1;
	}
    bool hasindex(const IQIndex& I) const
	{
        for(size_t j = 0; j < left.size(); ++j)
            if(left[j] == I) return true;
        return false;
	}
    bool hasindex(const Index& i) const
	{
        for(size_t j = 0; j < left.size(); ++j)
            if(left[j].hasindex(i)) return true;
        return false;
	}
    int num_left() const { return int(left.size()); }

    void conj() 
    { 
        if(!initted) Error("IQCombiner::conj(): IQCombiner not initialized.");
        foreach(IQIndex& I, left) I.conj(); 
        (do_condense ? cindex : _right).conj(); 
    }

    inline friend std::ostream& operator<<(std::ostream & s, const IQCombiner & c)
    {
        if(c.isInit())
            { s << std::endl << "right is " << c.right() << "\n"; }
        else
            { s << std::endl << "right is not initialized\n"; }
        s << "lefts are \n";
        foreach(const IQIndex& I, c.left) s << I << std::endl;
        return s << "\n\n";
    }
    IQTensor operator*(const IQTensor& t) const { IQTensor res; product(t,res); return res; }
    friend inline IQTensor operator*(const IQTensor& t, const IQCombiner& c) { return c.operator*(t); }

    void product(const IQTensor& t, IQTensor& res) const
    {
        init();
        std::vector<IQIndex> iqinds;

        int j;
        //t has right IQIndex, expand it
        if((j = t.findindex(_right)) != 0)
        {
            IQTensor t_uncondensed;
            if(do_condense) 
            { 
                cond.product(t,t_uncondensed); 
                j = t_uncondensed.findindex(cindex);
            }
            const IQTensor& t_ = (do_condense ? t_uncondensed : t);
            const IQIndex& r = (do_condense ? cindex : _right);

            if(t_.index(j).dir() == r.dir())
            {
                std::cerr << "IQTensor = " << t_ << std::endl;
                std::cerr << "IQCombiner = " << *this << std::endl;
                std::cerr << "IQIndex from IQTensor = " << t_.index(j) << std::endl;
                std::cerr << "(Right) IQIndex from IQCombiner = " << r << std::endl;
                Error("Incompatible arrow directions in operator*(IQTensor,IQCombiner).");
            }
            copy(t_.const_iqind_begin(),t_.const_iqind_begin()+j-1,std::back_inserter(iqinds));
            copy(left.begin(),left.end(),std::back_inserter(iqinds));
            copy(t_.const_iqind_begin()+j,t_.const_iqind_end(),std::back_inserter(iqinds));

            res = IQTensor(iqinds);

            foreach(const ITensor& it, t_.itensors())
            for(int k = 1; k <= it.r(); ++k)
            if(r.hasindex(it.index(k)))
            { res += (it * rightcomb[it.index(k)]); }

        }
        else
        {
            //t has left IQIndex's, combine them

            //res will have all IQIndex's of t not in the left of c
            foreach(const IQIndex& I, t.iqinds()) 
            { if(!hasindex(I)) iqinds.push_back(I); }
            //and res will have c's right IQIndex
            if(do_condense) iqinds.push_back(cindex);
            else            iqinds.push_back(_right);

            res = IQTensor(iqinds);

            for(std::vector<IQIndex>::const_iterator I = left.begin(); I != left.end(); ++I)
            {
                if((j = t.findindex(*I)) == 0)
                {
                    t.printIQInds("t");
                    std::cerr << "Left indices\n";
                    for(size_t j = 0; j < left.size(); ++j)
                    { std::cerr << j SP left[j] << "\n"; }
                    
                    Error("bad IQCombiner IQTensor product");
                }
                else //IQIndex is in left
                if(t.index(j).dir() == I->dir())
                {
                    std::cerr << "IQTensor = " << t << std::endl;
                    std::cerr << "IQCombiner = " << *this << std::endl;
                    std::cerr << "IQIndex from IQTensor = " << t.index(j) << std::endl;
                    std::cerr << "(Left) IQIndex from IQCombiner = " << *I << std::endl;
                    Error("Incompatible arrow directions in operator*(IQTensor,IQCombiner).");
                }
            }

            for(IQTensor::const_iten_it i = t.const_iten_begin(); i != t.const_iten_end(); ++i)
            {
                Real rse = 0;
                for(int k = 1; k <= i->r(); ++k)
                {
                if(hasindex(i->index(k))) 
                { rse += i->index(k).unique_Real(); }
                }

                if(setcomb.count(rse) == 0)
                {
                    Print(*i);
                    std::cerr << "\nleft indices \n";
                    for(size_t j = 0; j < left.size(); ++j)
                        { std::cerr << j << " " << left[j] << "\n"; }
                    std::cerr << "\n\n";
                    for(std::map<ApproxReal, Combiner>::const_iterator uu = setcomb.begin();
                        uu != setcomb.end(); ++uu)
                    {
                        std::cout << "Combiner: " << std::endl;
                        std::cout << uu->second << std::endl;
                    }
                    Error("no setcomb for rse in IQCombiner prod");
                }
                res += (*i * setcomb[rse]);
            }
            if(do_condense) { IQTensor rcopy(res); cond.product(rcopy,res); }
        }
    } //void product(const IQTensor& t, IQTensor& res) const


};


#endif
