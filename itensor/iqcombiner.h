#ifndef __IQCOMBINER_H
#define __IQCOMBINER_H
#include "combiner.h"
#include "iqtensor.h"

class Condenser	// Within one IQIndex, combine indices, presumably with same QNs
{
    IQIndex bigind_, smallind_;		// uncondensed, condensed
    mutable map<Index, pair<Index,int> > big_to_small;
    mutable map<pair<Index,int>,Index > small_to_big;

    // Use connections in t to create groupings; big = uncondensed, small = cond
    void init(const IQIndex& _bigind, const string& smallind_name)
    {
        /* May not be appropriate when orthogonalizing MPOs
           Could be a hint that there is a better way...
        if(_bigind.dir() != _smallind.dir())
        {
            cerr << "_bigind = " << _bigind << endl;
            cerr << "_smallind = " << _smallind << endl;
            Error("Arrow dirs not the same in Condenser.");
        }
        */
        static vector<QN> qns(1000);
        qns.resize(0);
        foreach(const inqn& x, bigind_.iq()) qns.push_back(x.qn);
        sort(qns.begin(),qns.end());
        vector<QN>::iterator ue = unique(qns.begin(),qns.end());

        vector<inqn> iq;
        for(vector<QN>::iterator qi = qns.begin(); qi != ue; ++qi)
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
                small_to_big[make_pair(small_qind,start)] = xi;
                big_to_small[xi] = make_pair(small_qind,start);
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

    Condenser(const IQIndex& _bigind, const string& smallind_name)
    : bigind_(_bigind)
    { init(_bigind,smallind_name); }

    IQTensor operator*(const IQTensor& t) { IQTensor res; product(t,res); return res; }
    friend inline IQTensor operator*(const IQTensor& t, const Condenser& c) { IQTensor res; c.product(t,res); return res; }

    void product(const IQTensor& t, IQTensor& res) const
    {
        assert(&t != &res);
        assert(smallind_.is_not_null());
        assert(bigind_.is_not_null());
        vector<IQIndex> iqinds; iqinds.reserve(t.r());
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
                    Index bind = small_to_big[make_pair(sind,start)];
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
                    pair<Index,int> Ii = big_to_small[it.index(k)];
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
        res.addindex1(t.virtual_ind());
    }

    inline friend ostream& operator<<(ostream & s, const Condenser & c)
    {
        s << "bigind_ is " << c.bigind_ << "\n";
        s << "smallind_ is " << c.smallind_ << "\n";
        s << "big_to_small is " << "\n";
        for(map<Index, pair<Index,int> >::const_iterator kk = c.big_to_small.begin();
            kk != c.big_to_small.end(); ++kk)
            { s << kk->first SP kk->second.first SP kk->second.second << "\n"; }
        s << "small_to_big is " << endl;
        for(map<pair<Index,int>,Index>::const_iterator kk = c.small_to_big.begin();
            kk != c.small_to_big.end(); ++kk)
            { s << kk->first.first SP kk->first.second SP kk->second << "\n"; }
        return s << endl;
    }
}; //class Condenser

namespace { //Anonymous namespace means don't expose outside of this header

class QCounter
{
public:
    vector<int> n;
    vector<int> ind;
    bool don;
    QCounter(const vector<IQIndex>& v)
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
            ind = vector<int>(nn,0);
            don = true;
        }

        return *this;
	}

    void getVecInd(const vector<IQIndex>& v, vector<Index>& vind, QN& q) const
	{
        q = QN(); vind.clear();
        for(unsigned int i = 0; i < ind.size(); ++i)
        {
            const int j = ind[i]+1;
            if(GET(v,i).nindex() < j)
            {
                for(unsigned int k = 0; k < n.size(); ++k) cerr << boost::format("n[%d] = %d\n")%k%n[k];
                cout << boost::format("i=%d, j=%d, v[i].nindex()=%d\n")%i%j%v[i].nindex();
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
    vector<IQIndex> left;
    mutable IQIndex _right;
    mutable map<ApproxReal, Combiner> setcomb;
    mutable map<Index, Combiner> rightcomb;
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
	    const IQIndex& l1, const IQIndex& l2 = IQIndNull, 
        const IQIndex& l3 = IQIndNull, const IQIndex& l4 = IQIndNull, 
	    const IQIndex& l5 = IQIndNull, const IQIndex& l6 = IQIndNull )
        : initted(false), do_condense(false)
	{
        if(l1 != IQIndNull) left.push_back(l1); 
        if(l2 != IQIndNull) left.push_back(l2);
        if(l3 != IQIndNull) left.push_back(l3); 
        if(l4 != IQIndNull) left.push_back(l4);
        if(l5 != IQIndNull) left.push_back(l5); 
        if(l6 != IQIndNull) left.push_back(l6);
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

    inline bool check_init() const { return initted; }

    // Initialize after all lefts are there and before being used
    void init(string rname = "combined", IndexType = Link, 
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
        vector<inqn> iq;
        for( ; c.notdone(); ++c)
        {
            vector<Index> vind;
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
            string cname = "cond::" + rname;
            cond = Condenser(cindex,cname);
            _right = cond.smallind();
        }
        else _right = IQIndex(rname,iq,rdir,plev);

        initted = true;
	}
    
    operator IQTensor() const
    {
        if(!initted) Error("IQCombiner::operator IQTensor(): IQCombiner not initialized.");

        if(_right.m() > 16) 
        { cerr << endl << endl << "WARNING: too large of an m in IQCombiner::operator IQTensor(). May be inefficient!" << endl << endl; }

        vector<IQIndex> iqinds(left);
        iqinds.push_back((do_condense ? cindex : _right));
        IQTensor res(iqinds);
        for(map<ApproxReal,Combiner>::const_iterator it = setcomb.begin();
            it != setcomb.end(); ++it)
        { res.insert(it->second); }

        //Combiners should always have the 
        //structure of zero divergence IQTensors
        assert(check_QNs(res));

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

    inline friend ostream& operator<<(ostream & s, const IQCombiner & c)
    {
        s << endl << "right is " << c.right() << endl;
        s << "lefts are " << endl;
        foreach(const IQIndex& I, c.left) s << I << endl;
        return s << "\n\n";
    }
    IQTensor operator*(const IQTensor& t) const { IQTensor res; product(t,res); return res; }
    friend inline IQTensor operator*(const IQTensor& t, const IQCombiner& c) { return c.operator*(t); }

    void product(const IQTensor& t, IQTensor& res) const
    {
        init();
        vector<IQIndex> iqinds;

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
                cerr << "IQTensor = " << t_ << endl;
                cerr << "IQCombiner = " << *this << endl;
                cerr << "IQIndex from IQTensor = " << t_.index(j) << endl;
                cerr << "(Right) IQIndex from IQCombiner = " << r << endl;
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

            res.addindex1(t_.virtual_ind());
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

            for(vector<IQIndex>::const_iterator I = left.begin(); I != left.end(); ++I)
            {
                if((j = t.findindex(*I)) == 0)
                {
                    t.printIQInds("t");
                    cerr << "Left indices\n";
                    for(size_t j = 0; j < left.size(); ++j)
                    { cerr << j SP left[j] << "\n"; }
                    
                    Error("bad IQCombiner IQTensor product");
                }
                else //IQIndex is in left
                if(t.index(j).dir() == I->dir())
                {
                    cerr << "IQTensor = " << t << endl;
                    cerr << "IQCombiner = " << *this << endl;
                    cerr << "IQIndex from IQTensor = " << t.index(j) << endl;
                    cerr << "(Left) IQIndex from IQCombiner = " << *I << endl;
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
                    cerr << "\nleft indices \n";
                    for(size_t j = 0; j < left.size(); ++j)
                        { cerr << j << " " << left[j] << "\n"; }
                    cerr << "\n\n";
                    for(map<ApproxReal, Combiner>::const_iterator uu = setcomb.begin();
                        uu != setcomb.end(); ++uu)
                    {
                        cout << "Combiner: " << endl;
                        cout << uu->second << endl;
                    }
                    Error("no setcomb for rse in IQCombiner prod");
                }
                res += (*i * setcomb[rse]);
            }
            res.addindex1(t.virtual_ind());
            if(do_condense) { IQTensor rcopy(res); cond.product(rcopy,res); }
        }
    } //void product(const IQTensor& t, IQTensor& res) const


};


#endif
