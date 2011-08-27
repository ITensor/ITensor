#ifndef __IQCOMBINER_H
#define __IQCOMBINER_H
#include "iq.h"

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
            if(v.at(i).nindex() < j)
            {
                for(unsigned int k = 0; k < n.size(); ++k) cerr << format("n[%d] = %d\n")%k%n[k];
                cout << format("i=%d, j=%d, v[i].nindex()=%d\n")%i%j%v[i].nindex();
                Error("bad v[i].iq in getVecInd");
            }
            vind.push_back(v[i].index(j));
            q += v[i].qn(j)*v[i].dir();
        }
	} //void QCounter::getVecInd
};

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
public:
    vector<IQIndex> left;
    IQIndex _right;
    mutable map<ApproxReal, Combiner> setcomb;
    mutable map<Index, Combiner> rightcomb;
    static IQIndex spec;
    bool initted;

    IQCombiner(
	    IQIndex l1 = spec, IQIndex l2 = spec, IQIndex l3 = spec, IQIndex l4 = spec, 
	    IQIndex l5 = spec, IQIndex l6 = spec )
	    : initted(false)
	{
        if(l1 != spec) left.push_back(l1); if(l2 != spec) left.push_back(l2);
        if(l3 != spec) left.push_back(l3); if(l4 != spec) left.push_back(l4);
        if(l5 != spec) left.push_back(l5); if(l6 != spec) left.push_back(l6);
	}
    void addleft(const IQIndex& l) 	// Include another left index
	{ 
        left.push_back(l);
        initted = false;
	}
    void check_init() const
	{
        if(!initted) Error("not initted");
	}
    void init(IQIndex& r) 		// Initialize after all lefts are there and before being used
	{
        if(initted) return;
        setcomb.clear();
        rightcomb.clear();

        //Flip around arrows of left IQIndices
        //This automatically makes combiner compatible
        //with IQTensor from which it got its left IQIndices
        foreach(IQIndex& l, left) l.conj();

        if(r.is_null()) Error("IQCombiner::init: Uninitialized right IQIndex.");

        //Construct individual Combiners
        QCounter c(left);
        vector<inqn> iq;
        for( ; c.notdone(); ++c)
        {
            vector<Index> vind;
            Index ii("combined");
            ii.primelevel = r.primelevel;
            QN q;
            c.getVecInd(left, vind, q);		// updates vind and q
            Combiner co(ii);
            Real rss = 0.0;
            foreach(const Index& j, vind)
            {
                co.addleft(ii,j);
                //cerr << format("Adding index (ur=%f) ")%j.unique_Real() << j << "\n";
                rss += j.unique_Real();
            }
            q *= -r.dir();
            assert(ii.primelevel == r.primelevel);
            iq.push_back(inqn(ii,q));
            //cerr << format("Inserting the following combiner with an rss=%f\n")%rss; cerr << co << "\n";
            setcomb.insert(make_pair(ApproxReal(rss),co));
            rightcomb.insert(make_pair(ii,co));
        }
        r = IQIndex(r,iq);
        _right = r;

        initted = true;
	}
    
    IQTensor toIQTensor() const
    {
        check_init();

        if(_right.m() > 16) 
        { cerr << endl << endl << "WARNING: too large of an m in IQCombiner::toIQTensor(). May be inefficient!" << endl << endl; }

        //Create a Kronecker delta IQTensor
        //between _right and _right.primed()
        //IQTensor Delta(_right.conj(),_right.primed());
        //foreach(inqn iq, _right.qindex)
        //{
            //ITensor del(iq.index,iq.index.primed());
            //for(int i = 1; i <= iq.index.m(); ++i)
                //del(i,i) = 1.0;
            //Delta.insert(del);
        //}
        IQIndex rc = _right; rc.conj();
        IQIndex rp = _right.primed();
        IQTensor Delta(rc,rp);
        for(int n = 1; n <= rc.nindex(); ++n)
        { Delta.insert(ITensor(rc.index(n),rp.index(n),1)); }

        //Use the delta tensor to convert
        //this IQCombiner into an IQTensor
        IQTensor res = Delta * (*this);

        //Remove the prime
        res.ind_inc_prime(_right,-1);

        //Combiners should always have the 
        //structure of zero divergence IQTensors
        assert(res.checkDivZero());

        return res;
    }

    const IQIndex& right() const { check_init(); return _right; }

    int findindex(const IQIndex& i) const
	{
        check_init();
        for(int j = 0; j < (int)left.size(); j++)
            if(left[j] == i) return j;
        return -1;
	}
    bool hasindex(const IQIndex& i) const
	{
        check_init();
        return findindex(i) != -1;
	}
    bool in_left(Index i) const
	{
        check_init();
        for(int j = 0; j < (int)left.size(); j++)
            if(left[j].hasindex(i)) return true;
        return false;
	}
    int num_left() const { return int(left.size()); }

    void conj() { _right.conj(); foreach(IQIndex& I, left) I.conj(); }

    inline friend ostream& operator<<(ostream & s, const IQCombiner & c)
    {
        s << endl << "right is " << c.right() << endl;
        s << "lefts are " << endl;
        foreach(const IQIndex& I, c.left) s << I << endl;
        return s << "\n\n";
    }
};

inline IQTensor operator*(const IQTensor& t, const IQCombiner& c)
{
    //cerr << "IQTensor*IQCombiner, multiplying\n";
    //t.print("t"); cerr << "c = " << c << "\n";
    int j;
    IQTensor res;

    //Ensure that Virtual IQIndex gets copied
    res.addindex1(t.virtual_ind());

    //t has right IQIndex, expand it
    if((j = t.findindex(c.right())) != -1)
	{
        if(t.iqindex(j).dir() == c.right().dir())
        {
            cerr << "IQTensor = " << t << endl;
            cerr << "IQCombiner = " << c << endl;
            cerr << "IQIndex from IQTensor = " << t.iqindex(j) << endl;
            cerr << "(Right) IQIndex from IQCombiner = " << c.right() << endl;
            Error("Incompatible arrow directions in operator*(IQTensor,IQCombiner).");
        }
        copy(t.iqindex_.begin(),t.iqindex_.begin()+j,std::back_inserter(res.iqindex_));
        copy(c.left.begin(),c.left.end(),std::back_inserter(res.iqindex_));
        copy(t.iqindex_.begin()+j+1,t.iqindex_.end(),std::back_inserter(res.iqindex_));

        foreach(const ITensor& it, t.itensors())
        for(int k = 1; k <= it.r(); ++k)
        if(c.right().hasindex(it.index(k)))
        { res += (it * c.rightcomb[it.index(k)]); }
        return res;
	}

    //t has left IQIndex's, combine them

    //res will have all IQIndex's of t not in the left of c
    foreach(const IQIndex& i, t.iqindex_) 
    { if(!c.hasindex(i)) res.iqindex_.push_back(i); }
    //and res will have c's right IQIndex
    res.iqindex_.push_back(c.right());

    foreach(const IQIndex& i, c.left)
    {
        if((j = t.findindex(i)) == -1)
        {
            cout << t << endl << c << endl;
            Error("bad IQCombiner IQTensor product");
        }
        else //IQIndex is in c.left
        {
            if(t.iqindex_[j].dir() == i.dir())
            {
                cerr << "IQTensor = " << t << endl;
                cerr << "IQCombiner = " << c << endl;
                cerr << "IQIndex from IQTensor = " << t.iqindex_[j] << endl;
                cerr << "(Left) IQIndex from IQCombiner = " << i << endl;
                Error("Incompatible arrow directions in operator*(IQTensor,IQCombiner).");
            }
        }
    }

    for(IQTensor::const_iten_it i = t.const_iten_begin(); i != t.const_iten_end(); ++i)
	{
        Real rse = 0.0;
        int d = i->r();
        for(int k = 1; k <= d; ++k)
        if(c.in_left(i->index(k)))
        { 
        //cerr << format("Adding index (ur=%f) ")%(i->index(k).unique_Real()) << i->index(k) << "\n";
        rse += i->index(k).unique_Real(); 
        }
        //else cerr << "Not adding index " << i->index(k) << "\n";

        //cerr << "rse = " << rse << "\n";

        if(c.setcomb.count(rse) == 0)
        {
            Printit<Index> pr(cout," ");
            cout << "c.left[0] is " << c.left[0];
            //cout << "se is " << endl;
            //for_all(se,pr); cout << endl;
            //for(map<set<Index>, Combiner>::const_iterator uu = c.setcomb.begin();
            for(map<ApproxReal, Combiner>::const_iterator uu = c.setcomb.begin();
                uu != c.setcomb.end(); ++uu)
            {
                //cout << "set members: " << endl;
                //for_all(uu->first,pr); 
                cout << "Combiner: " << endl;
                cout << uu->second << endl;
            }
            Error("no setcomb for se in IQCombiner prod");
        }
        //cout << "c.setcomb[se] is " << c.setcomb[se] << endl;
        res += (*i * c.setcomb[rse]);
	}
    //cout << "here is the result of IQTensor * IQCombiner" << endl << res;
    return res;
}

inline IQTensor operator*(const IQCombiner& c, const IQTensor& t) { return t * c; }


class Condenser	// Within one IQIndex, combine indices, presumably with same QNs
{
public:
    IQIndex bigind, smallind;		// uncondensed, condensed
    mutable map<Index, pair<Index,int> > big_to_small;
    mutable map<pair<Index,int>,Index > small_to_big;

    // Use connections in t to create groupings; big = uncondensed, small = cond
    //Condenser(const IQIndex& _bigind, IQIndex& _smallind,const IQTensor& t);

    Condenser(const IQIndex& _bigind, IQIndex& _smallind,const IQTensor& t)
        : bigind(_bigind) // Use connections in t to create groupings
    {
        if(_bigind.dir() != _smallind.dir())
        {
            cerr << "_bigind = " << _bigind << endl;
            cerr << "_smallind = " << _smallind << endl;
            Error("Arrow dirs not the same in Condenser.");
        }
        static vector<QN> qns(1000);
        qns.resize(0);
        foreach(const inqn& x, bigind.iq()) qns.push_back(x.qn);
        sort(qns.begin(),qns.end());
        vector<QN>::iterator ue = unique(qns.begin(),qns.end());

        vector<inqn> iq;
        for(vector<QN>::iterator qi = qns.begin(); qi != ue; ++qi)
        {
            QN& q = *qi;
            int totm = 0;
            foreach(const inqn& x, bigind.iq())
            if(x.qn == q) totm += x.index.m();

            Index small_qind("condensed",totm);
            int start = 0;
            foreach(const inqn& x, bigind.iq())
            if(x.qn == q)
            {
                const Index &bj(x.index);
                small_to_big[make_pair(small_qind,start)] = bj;
                big_to_small[bj] = make_pair(small_qind,start);
                start += bj.m();
            }
            iq.push_back(inqn(small_qind,q));
        }
        smallind = IQIndex(_smallind.name(),iq,_smallind.dir());
        _smallind = smallind;
    }

    inline friend ostream& operator<<(ostream & s, const Condenser & c)
    {
        s << "bigind is " << c.bigind << endl;
        s << "smallind is " << c.smallind << endl;
        s << "big_to_small is " << endl;
        for(map<Index, pair<Index,int> >::const_iterator kk = c.big_to_small.begin();
            kk != c.big_to_small.end(); ++kk)
        s << kk->first SP kk->second.first SP kk->second.second << endl;
        s << "small_to_big is " << endl;
        for(map<pair<Index,int>,Index>::const_iterator kk = c.small_to_big.begin();
            kk != c.small_to_big.end(); ++kk)
        s << kk->first.first SP kk->first.second SP kk->second << endl;
        return s << endl;
    }
};

inline ITensor converter(Index sind, Index bind, int start)	// sind is condensed from lots of binds
{
    if(sind.m() < bind.m()) 
	Error("sind.m() < bind.m()");
    ITensor temp(sind,bind);
    for(int k = 1; k <= bind.m(); k++)
	temp(start+k,k) = 1.0;
    return temp;
}

inline void doconvert(const ITensor& t, Index cond, Index uncond, int start, ITensor& res)
{
    vector<Index> indices; indices.reserve(t.r());
    for(int j = 1; j <= t.r(); ++j)
    {
    if(t.index(j) == uncond) { indices.push_back(cond); }
    else indices.push_back(t.index(j));
    }
    const bool do_allocate = true;
    res = ITensor(indices,do_allocate);

    int i = res.findindexn(cond);

    Counter c(t);
    int inc[NMAX+1];
    for(int j = 0; j <= NMAX; j++) inc[j] = 0;
    inc[i] = start;
    for( ; c != Counter::done ; ++c)
	{
	res(c.i[1]+inc[1],c.i[2]+inc[2],c.i[3]+inc[3],c.i[4]+inc[4],
	    c.i[5]+inc[5],c.i[6]+inc[6],c.i[7]+inc[7],c.i[8]+inc[8]) 
        = t.dat()(c.ind);
	}
}

inline IQTensor operator*(const IQTensor& t, const Condenser& c)
{
    int j;
    IQTensor res;
    res.iqindex_ = t.iqindex_; res.viqindex = t.viqindex;

    if((j = t.findindex(c.smallind)) != -1)	// expand condensed form into uncondensed
	{
        res.iqindex_[j] = c.bigind;
        for(IQTensor::const_iten_it i = t.const_iten_begin(); i != t.const_iten_end(); ++i)
        {
            int k;
            for(k = 1; k <= i->r(); k++)
            if(c.smallind.hasindex(i->index(k)))
                break;

            Index sind = i->index(k);
            for(int start = 0; start < sind.m(); )
            {
                Index bind = c.small_to_big[make_pair(sind,start)];
                res += (*i * converter(sind,bind,start));
                start += bind.m();
            }
        }
        return res;
	}

    //else contract regular form into condensed

    if((j = t.findindex(c.bigind)) == -1)
	{
        cout << "t = " << t;
        cout << "c = " << c;
        Error("bad call to IQTensor * Condenser");
	}
    //cout << "t = " << t;
    //cout << "c.bigind is " << c.bigind << endl;
    //cout << "j which is location of c.bigind is " << j << endl;
    res.iqindex_[j] = c.smallind;
    foreach(const ITensor& i, t.itensors())
	{
        bool gotit = false;
        int k;
        for(k = 1; k <= i.r(); ++k)
        if(c.bigind.hasindex(i.index(k)))
        {
            pair<Index,int> pp = c.big_to_small[i.index(k)];
            ITensor tt;
            doconvert(i,pp.first,i.index(k),pp.second,tt);
            res += tt;
            gotit = true;
            break;
        }
        if(!gotit)
        {
            cout << c << endl;
            cout << "i = " << i << endl;
            Error("cant find index");
        }
	}
    return res;
}

inline IQTensor operator*(const Condenser& c, const IQTensor& t)
{ return t * c; }


#ifdef THIS_IS_MAIN
IQIndex IQCombiner::spec("spec");
#endif

#endif
