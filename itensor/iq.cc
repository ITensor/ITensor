#include "iq.h"

ostream& operator<<(ostream &o, const QN& q)
{ return o<< format("sz = %d, Nf = %d, fp = %s") % q.sz() % q.Nf() % (q.fp() < 0 ? "-" : "+"); }

ostream& operator<<(ostream & s, const IQIndex& I)
{
    s << "IQIndex: " << (const Index&) I << " <" << I.dir() << ">" << endl;
    for(int j = 1; j <= I.nindex(); ++j) s << " " << I.index(j) SP I.qn(j) << "\n";
    return s;
}

ostream& operator<<(ostream & s, const IQTensor &t)
{
    s << "\n----- IQTensor -----\nIQIndices: " << endl;
    for(unsigned int i = 0; i < t.iqindex.size(); i++)
        s << "  " << t.iqindex[i] << endl;
    if(t.has_virtual()) s << "  " << t.viqindex << endl;
    s << "ITensors: " << endl;
    for(IQTensor::const_iten_it i = t.const_iten_begin(); 
	    i != t.const_iten_end(); ++i)
	s <<"	" << *i << endl;
    s << "-------------------" << endl << endl;
    return s;
}

ostream& operator<<(ostream & s, const IQCombiner & c)
    {
    s << endl << "right is " << c.right() << endl;
    s << "lefts are " << endl;
    foreach(IQIndex I, c.left)
        s << I << endl;
    //Printit<IQIndex> p(s,"   ");
    //for_all(c.left,p);
    s << endl;
    return s << endl;
    }

ostream& operator<<(ostream & s, const Condenser & c)
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

void DoPrimer::operator()(IQIndex &iqi) const { iqi.doprime(pt,inc); }
void MapPrimer::operator()(IQIndex &iqi) const { iqi.mapprime(plevold,plevnew,pt); }

IQIndexVal IQIndex::operator()(int n) const { return IQIndexVal(*this,n); }

void IQTensor::SplitReIm(IQTensor& re, IQTensor& im) const
{
    if(!hasindex(IQIndReIm))
	{
	IQTensor cop(*this);
	re = cop;
	im = cop;
	im *= 0.0;
	return;
	}
    re = IQTensor();
    remove_copy_if(iqindex.begin(),iqindex.end(),std::back_inserter(re.iqindex),
		    bind2nd(std::equal_to<IQIndex>(),IQIndReIm));
    im = re;
    ITensor a,b;
    for(const_iten_it i = itensor.begin(); i != itensor.end(); ++i)
	{
	i->SplitReIm(a,b);
	re.insert(a);
	im.insert(b);
	}
}

IQTensor& IQTensor::operator*=(const IQTensor& other)
{
    if(hasindex(IQIndReIm) && other.hasindex(IQIndReIm) && !other.hasindex(IQIndReImP)
	    && !other.hasindex(IQIndReImPP))
	{
        static ITensor primer(IndReIm,IndReImP,1.0);
        static ITensor primerP(IndReIm,IndReImPP,1.0);
        static ITensor prod(IndReIm,IndReImP,IndReImPP);
        static IQTensor iqprimer(IQIndReIm,IQIndReImP);
        static IQTensor iqprimerP(IQIndReIm,IQIndReImPP);
        static IQTensor iqprod(IQIndReIm,IQIndReImP,IQIndReImPP);
        static int first = 1;
        if(first)
        {
            IndexVal iv0(IndReIm,1), iv1(IndReImP,1), iv2(IndReImPP,1);
            iv0.i = 1; iv1.i = 1; iv2.i = 1; prod(iv0,iv1,iv2) = 1.0;
            iv0.i = 1; iv1.i = 2; iv2.i = 2; prod(iv0,iv1,iv2) = -1.0;
            iv0.i = 2; iv1.i = 2; iv2.i = 1; prod(iv0,iv1,iv2) = 1.0;
            iv0.i = 2; iv1.i = 1; iv2.i = 2; prod(iv0,iv1,iv2) = 1.0;
            first = 0;
            iqprimer += primer;
            iqprimerP += primerP;
            iqprod += prod;
        }
        return *this = (*this * iqprimer) * iqprod * (other * iqprimerP);
	}

        /*
        if(0)
        {
        int numsame = 0;
        order_iqindex_for_product(*this,other,numsame);
        if(0)
            {
            for(vector<IQIndex>::const_iterator i = iqindex.begin(); i != iqindex.end(); ++i)
            cout << i->unique_Real() << " ";
            cout << " | ";
            for(vector<IQIndex>::const_iterator i = other.iqindex.begin(); i != other.iqindex.end(); ++i)
            cout << i->unique_Real() << " ";
            cout << endl;
            }
        match_order();
        other.match_order();
        }
        */

    IQTensor res;

    //Combine virtual indices
    if(viqindex == IQEmptyV)
    {
        if(other.viqindex == IQEmptyV)
            res.viqindex = IQEmptyV;
        else
            res.viqindex = other.viqindex;
    }
    else 
    {
        if(other.viqindex == IQEmptyV)
            res.viqindex = viqindex;
        else
        {
            //Add virtual IQIndex's
            res.viqindex = viqindex; //Inherits its dir from *this
            QN newq = 
                res.viqindex.dir()*(viqindex.dir()*viqindex.qn(1)+other.viqindex.dir()*other.viqindex.qn(1));
            res.viqindex.set_qn(1,newq);
        }
    }
    
    
    //Load res.iqindex with those IQIndex's *not* common to *this and other
    static vector<IQIndex> riqind_holder(1000);
    riqind_holder.resize(0);
    foreach(const IQIndex& I, iqindex)
    {   
        if(!has_element<IQIndex,vector<IQIndex> >(I,other.iqindex))
            riqind_holder.push_back(I);
    }                
    foreach(const IQIndex& I, other.iqindex)
    {   
        if(!has_element<IQIndex,vector<IQIndex> >(I,iqindex))
            riqind_holder.push_back(I);
    }
    if(riqind_holder.size() > 1000) cerr << endl << "WARNING: in IQTensor::operator* riqind_holder had to reallocate." << endl << endl;
    res.iqindex = riqind_holder;


    //Fill common_inds with IQIndex's common to both
    vector<IQIndex> common_inds;
    for(vector<IQIndex>::const_iterator i = iqindex.begin(); i != iqindex.end(); ++i)
    {
        vector<IQIndex>::const_iterator f = find(other.iqindex.begin(),other.iqindex.end(),*i);
        if(f != other.iqindex.end()) //*i is an element of other.iqindex
        {
            common_inds.push_back(*i);
            //Check that arrow directions are compatible
            if(f->dir() == i->dir() && f->type() != ReIm && i->type() != ReIm)
            {
                cerr << "*this = " << *this << endl;
                cerr << "other = " << other << endl;
                cerr << "IQIndex from *this = " << *i << endl;
                cerr << "IQIndex from other = " << *f << endl;
                Error("Incompatible arrow directions in IQTensor::operator*.");
            }
        }
    }

    //is_common returns true for the unique_Real of an Index 
    //if it is part of a common IQIndex
    map<ApproxReal,bool> is_common;
    foreach(const IQIndex& i, iqindex)
	{
        bool iscom = has_element(i,common_inds);
        foreach(const inqn& x, i.iq()) is_common[ApproxReal(x.index.unique_Real())] = iscom;
	}
    foreach(const IQIndex& i, other.iqindex)
	{
        bool iscom = has_element(i,common_inds);
        foreach(const inqn& x, i.iq()) is_common[ApproxReal(x.index.unique_Real())] = iscom;
	}

    multimap<ApproxReal,const_iten_it> com_this;
    static vector<ApproxReal> keys(3000);
    keys.resize(0);

    //com_this maps the unique_Real of a set of Index's to be contracted over together
    //to those ITensors in *this.itensor having all Index's in that set
    for(const_iten_it tt = const_iten_begin(); tt != const_iten_end(); ++tt)
	{
        Real r = 0.0;
        for(int a = 1; a <= tt->r(); ++a)
        {
            if(is_common[ApproxReal(tt->index(a).unique_Real())])
            {
                r += tt->index(a).unique_Real();
            }
        }
        com_this.insert(make_pair(ApproxReal(r),tt));
        keys.push_back(ApproxReal(r));
	}

    //com_other is the same as com_this but for other
    multimap<ApproxReal,const_iten_it> com_other;
    for(const_iten_it ot = other.const_iten_begin(); ot != other.const_iten_end(); ++ot)
	{
        Real r = 0.0;
        for(int b = 1; b <= ot->r(); ++b)
        {
            if(is_common[ApproxReal(ot->index(b).unique_Real())])
            {
                r += ot->index(b).unique_Real();
            }
        }
        com_other.insert(make_pair(ApproxReal(r),ot));
        keys.push_back(ApproxReal(r));
	}
    if(keys.size() > 3000) cerr << endl << format("WARNING: in IQTensor::operator* keys had to reallocate (new size = %d).")%keys.size() << endl << endl;

    //Remove redundant keys and return an iterator pointing to the 
    //new end of the sequence
    sort(keys.begin(),keys.end());
    vector<ApproxReal>::iterator kend = unique(keys.begin(),keys.end());
    typedef multimap<ApproxReal,const_iten_it>::iterator mit;

    pair<mit,mit> lrange,rrange;
    ITensor tt;
    for(vector<ApproxReal>::iterator k = keys.begin(); k != kend; ++k)
	{
        //Equal range returns the begin and end iterators for the sequence
        //corresponding to multimap[key] as a pair
        lrange = com_this.equal_range(*k);
        rrange = com_other.equal_range(*k);

        //Iterate over all ITensors in *this and other sharing
        //the set of contracted Index's corresponding to k
        for(mit ll = lrange.first; ll != lrange.second; ++ll)
        for(mit rr = rrange.first; rr != rrange.second; ++rr)
        {
            //Multiply the ITensors and add into res
            tt = *(ll->second); tt *= *(rr->second);
            res += tt;
        }
	}

    return (*this = res);
} //IQTensor& IQTensor::operator*=(const IQTensor& other) const

//Extracts the real and imaginary parts of the 
//component of a rank 0 tensor (scalar)
void IQTensor::GetSingComplex(Real& re, Real& im) const
{
    IQTensor tre,tim;
    SplitReIm(tre,tim);

    //Only IQIndex should be IQIndReIm
    /*
    if(tre.iqindex.size() != 1)
	{
        cout << *this;
        cout << tre;
        Error("bad tre size");
	}
    if(tim.iqindex.size() != 1) Error("bad tim size");
    */
    foreach(const IQIndex& I, tre.iqindex)
    {
        if(I.type() != Virtual && I != IQTSing.iqindex[0])
        {
            cout << *this;
            cout << tre;
            Error("bad tre size");
        }
    }
    foreach(const IQIndex& I, tim.iqindex)
    if(I.type() != Virtual && I != IQTSing.iqindex[0])
    { Error("bad tim size"); }

    if(tre.iten_size() == 0)
    { re = 0.0; }
    else
	{
        if(tre.iten_begin()->dat().Length() != 1) 
        {
            cout << "tre is\n" << tre << endl;
            Error("bad tre dat size");
        }
        re = tre.itensor.begin()->val0() * exp(tre.itensor.begin()->logfac());
	}
    if(tim.iten_size() == 0)
	{ im = 0.0; }
    else
	{
        if(tim.itensor.begin()->dat().Length() != 1) Error("bad tim dat size");
        im = tim.itensor.begin()->val0() * exp(tim.itensor.begin()->logfac());
	}
}

void IQTensor::Assign(const IQTensor& other)
{
    //other.print("other");
    //Permutation P = getpermBtoA(*this,other);
    //Reshape(other,P,*this);
    if(0) //catch this later for a specific ITensor
    if(iten_size() < other.iten_size())
	{
        cout << "iten sizes: " << iten_size() SP other.iten_size() << endl;
        cout << *this;
        cout << other;
        Error("bad Assign sizes");
	}
    //map<set<Index>,list<ITensor>::iterator> semap;
    map<ApproxReal,list<ITensor>::iterator> semap;
    for(iten_it i = itensor.begin(); i != itensor.end(); ++i)
	{
        semap[ApproxReal(i->unique_Real())] = i;
	}
    for(const_iten_it i = other.itensor.begin(); i != other.itensor.end(); ++i)
	{
        ApproxReal se = ApproxReal(i->unique_Real());
        if(semap.count(se) == 0)
        {
            cout << "warning Assign semap.count is 0" << endl;
            cerr << "offending ITensor is " << *i << "\n";
            Error("bad Assign count se");
        }
        else
            semap[se]->Assign(*i);
	}
}

IQTensor& IQTensor::operator+=(const IQTensor& other)
{
    if(iqindex.size() == 0)		// Automatic initializing a summed IQTensor in a loop
        return *this = other;
    bool complex_this = hasindex(IQIndReIm); 
    bool complex_other = other.hasindex(IQIndReIm); 
    IQTensor& This(*this);
    if(!complex_this && complex_other)
        return (This = (This * IQComplex_1) + other);
    if(complex_this && !complex_other)
        return (This += other * IQComplex_1);
    Real ur1 = This.unique_Real();
    Real ur2 = other.unique_Real();
    if(fabs(ur1-ur2) > 1.0e-11) 
	{
        cout << "This is " << This;
        cout << "other is " << other;
        Error("bad match unique real in IQTensor::operator+=");
	}
    for(const_iten_it i = other.const_iten_begin(); i != other.const_iten_end(); ++i)
        operator+=(*i);
    return *this;
}

IQIndex index_in_common(const IQTensor& A, const IQTensor& B, IndexType t)
{
    foreach(const IQIndex& I, A.iqindex)
    if(I.type() == t) if(B.hasindex(I)) return I;

    return IQIndex();
}

IQTensor operator*(const IQTensor& t, const IQCombiner& c)
{
    //cerr << "IQTensor*IQCombiner, multiplying\n";
    //t.print("t"); cerr << "c = " << c << "\n";
    int j;
    IQTensor res;

    //Ensure that Virtual IQIndex gets copied
    res.viqindex = t.viqindex;

    //t has right IQIndex, expand it
    if((j = t.findindex(c.right())) != -1)
	{
        if(t.iqindex[j].dir() == c.right().dir())
        {
            cerr << "IQTensor = " << t << endl;
            cerr << "IQCombiner = " << c << endl;
            cerr << "IQIndex from IQTensor = " << t.iqindex[j] << endl;
            cerr << "(Right) IQIndex from IQCombiner = " << c.right() << endl;
            Error("Incompatible arrow directions in operator*(IQTensor,IQCombiner).");
        }
        copy(t.iqindex.begin(),t.iqindex.begin()+j,std::back_inserter(res.iqindex));
        copy(c.left.begin(),c.left.end(),std::back_inserter(res.iqindex));
        copy(t.iqindex.begin()+j+1,t.iqindex.end(),std::back_inserter(res.iqindex));
        for(IQTensor::const_iten_it i = t.const_iten_begin(); i != t.const_iten_end(); ++i)
        {
            int d = i->r();
            for(int k = 1; k <= d; ++k)
            if(c.right().hasindex(i->index(k)))
                res += (*i * c.rightcomb[i->index(k)]);
        }
        return res;
	}

    //t has left IQIndex's, combine them

    //res will have all IQIndex's of t not in the left of c
    foreach(const IQIndex& i, t.iqindex) 
    { if(!c.hasindex(i)) res.insert(i); }
    //and res will have c's right IQIndex
    res.insert(c.right());

    foreach(const IQIndex& i, c.left)
    {
        if((j = t.findindex(i)) == -1)
        {
            cout << t << endl << c << endl;
            Error("bad IQCombiner IQTensor product");
        }
        else //IQIndex is in c.left
        {
            if(t.iqindex[j].dir() == i.dir())
            {
                cerr << "IQTensor = " << t << endl;
                cerr << "IQCombiner = " << c << endl;
                cerr << "IQIndex from IQTensor = " << t.iqindex[j] << endl;
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

IQTensor kronecker_delta(const IQIndex& L, const IQIndex& R)
{
    if(L.m() != R.m()) Error("kronecker_delta: IQIndices must have the same m.");
    if(L.nindex() != R.nindex()) Error("kronecker_delta: IQIndices must have the same number of QN sectors.");
    IQTensor Delta(L,R);
    for(int n = 1; n <= L.nindex(); ++n)
        Delta.insert(ITensor(L.index(n),R.index(n),1));
    return Delta;
}

IQTensor IQCombiner::toIQTensor() const
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
    IQTensor Delta = kronecker_delta(rc,_right.primed());

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

template<class T, class C>
int get_loc(const C& c, const T& t)
    { return find(c.begin(),c.end(),t) - c.begin(); }

int num_diff(set<Index>& sa, set<Index>& sb)
    {
    vector<Index> counter;
    set_intersection(sa.begin(),sa.end(),sb.begin(),sb.end(),
	    std::back_inserter(counter));
    return (int) sa.size() - (int)counter.size();
    }

Condenser::Condenser(const IQIndex& _bigind, IQIndex& _smallind,const IQTensor& t)
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
	qns.resize(ue-qns.begin());

    vector<inqn> iq;
    foreach(QN& q, qns)
    {
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


ITensor converter(Index sind, Index bind, int start)	// sind is condensed from lots of binds
{
    if(sind.m() < bind.m()) 
	Error("sind.m() < bind.m()");
    ITensor temp(sind,bind);
    for(int k = 1; k <= bind.m(); k++)
	temp(start+k,k) = 1.0;
    return temp;
}

void doconvert(const ITensor& t, Index cond, Index uncond, int start, ITensor& res)
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

IQTensor operator*(const IQTensor& t, const Condenser& c)
{
    int j;
    IQTensor res;
    res.iqindex = t.iqindex; res.viqindex = t.viqindex;

    if((j = t.findindex(c.smallind)) != -1)	// expand condensed form into uncondensed
	{
        res.iqindex[j] = c.bigind;
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
    res.iqindex[j] = c.smallind;
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
