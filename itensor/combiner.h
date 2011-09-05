#ifndef __COMBINER_H
#define __COMBINER_H

#include "itensor.h"

/*
Combine several indices into one, use * to convert tensors efficiently
   \
    \
  ---C====
    /
   /

*/
class Combiner
{
    array<Index,NMAX+1> _leftn; // max dim is 8
    vector<Index> _left1;
    mutable Index _right;
    int _rln; //Number of m>1 'left' indices (indices to be combined into one)
    mutable bool initted;
public:

    //Accessor Methods ----------------------------------------------

    void doCondense(bool val) { } //For interface compatibility with IQCombiner

    const Index& right() const 
    { 
        if(!initted) Error("Combiner: right ind requested prior to init");
        return _right; 
    }
    int rln() const { return _rln; }
    const Index& leftn(int j) const { return GET(_leftn,j); }

    typedef array<Index,NMAX+1>::const_iterator leftn_it;
    const pair<leftn_it,leftn_it> leftn() const { return make_pair(_leftn.begin()+1,_leftn.begin()+_rln+1); }
    const vector<Index>& left1() const { return _left1; }

    //Constructors --------------------------------------------------

    Combiner() : _rln(0), initted(false) {}
    Combiner(const Index& l1 , const Index& l2 = IndNull, const Index& l3 = IndNull, const Index& l4 = IndNull, 
	    const Index& l5 = IndNull, const Index& l6 = IndNull, const Index& l7 = IndNull, const Index& l8 = IndNull )
    : initted(false)
	{
        _leftn[0] = IndNull;

        //Split given left indices into m==1 and m>1
        _rln = 0;
        if(l1 != IndNull) { if(l1.m() == 1) _left1.push_back(l1); else _leftn[++_rln] = l1; }
        if(l2 != IndNull) { if(l2.m() == 1) _left1.push_back(l2); else _leftn[++_rln] = l2; }
        if(l3 != IndNull) { if(l3.m() == 1) _left1.push_back(l3); else _leftn[++_rln] = l3; }
        if(l4 != IndNull) { if(l4.m() == 1) _left1.push_back(l4); else _leftn[++_rln] = l4; }
        if(l5 != IndNull) { if(l5.m() == 1) _left1.push_back(l5); else _leftn[++_rln] = l5; }
        if(l6 != IndNull) { if(l6.m() == 1) _left1.push_back(l6); else _leftn[++_rln] = l6; }
        if(l7 != IndNull) { if(l7.m() == 1) _left1.push_back(l7); else _leftn[++_rln] = l7; }
        if(l8 != IndNull) { if(l8.m() == 1) _left1.push_back(l8); else _leftn[++_rln] = l8; }
        init();
	}
    
    //Operators -----------------------------------------------------

    ITensor operator*(const ITensor& t) const { ITensor res; product(t,res); return res; }
    friend inline ITensor operator*(const ITensor& t, const Combiner& c) { return c.operator*(t); }

    //Index Methods -------------------------------------------------

    inline bool check_init() const { return initted; }

    void addleft(const Index& l)// Include another left index
	{ 
        initted = false;
        if(_rln == NMAX) Error("Combiner: already reached max number of left indices.");
        if(l.m() == 1) { _left1.push_back(l); return; } else GET(_leftn,++_rln) = l; 
	}

    //Initialize after all lefts are added and before being used
    void init(string rname = "combined", IndexType type = Link, int primelevel = 0) const
	{
        if(initted) return;
        int m = 1; for(int i = 1; i <= _rln; ++i) { m *= GET(_leftn,i).m(); }
        _right = Index(rname,m,type,primelevel); 
        initted = true;
    }

    int findindexn(Index i) const
	{
        for(int j = 1; j <= _rln; ++j)
            if(GET(_leftn,j) == i) return j;
        return 0;
	}
    bool hasindex(Index i) const
	{
        if(i.m() == 1)
        {
            foreach(const Index& L, _left1)
            if(i == L) return true;
            return false;
        }
        for(int j = 1; j <= _rln; ++j) if(GET(_leftn,j) == i) return true;
        return false;
	}

    //Other Methods -------------------------------------------------

    operator ITensor() const
    {
        if(_right.m() > 16) 
        { 
            Print(*this);
            Error(""); }
        //{ cerr << "\n\n" << "WARNING: too large of an m in Combiner::operator ITensor(). May be inefficient!\n\n"; }

        //Use a kronecker delta tensor to convert this Combiner into an Tensor
        ITensor res = operator*(ITensor(_right,_right.primed(),1));
        res.noprimeind(_right.primed());
        return res;
    }

    friend inline ostream & operator<<(ostream & s, const Combiner & c)
    {
        s << "\nRight index: " << c.right() << "\n";
        s << "Left indices:\n";
        foreach(const Index& l, c.leftn()) s << "	" << l << "\n";
        foreach(const Index& l, c.left1()) s << "	" << l << "\n";
        return s;
    }

    void conj() { }

    void product(const ITensor& t, ITensor& res) const
    {
        init();

        int j;
        if((j = t.findindex1(_right)) != 0)
        {
            res = t;
            res.removeindex1(j);
            //All of c's left indices must be m==1, so add them all
            foreach(const Index& I, _left1) res.addindex1(I);
            return;
        }
        else if((j = t.findindexn(_right)) != 0)
        {
            vector<Index> nindices; nindices.reserve(t.rn()+_rln-1);
            for(int i = 1; i < j; ++i)
                nindices.push_back(t.index(i));
            for(int i = 1; i <= _rln; ++i)
                nindices.push_back(GET(_leftn,i));
            for(int i = j+1; i <= t.rn(); ++i)
                nindices.push_back(t.index(i));
            foreach(const Index& I, _left1)     nindices.push_back(I);
            foreach(const Index& I, t.index1()) nindices.push_back(I);
            res = ITensor(nindices,t);
            return;
        }

        vector<Index> nindices; nindices.reserve(t.rn()-_rln+1);
        Permutation P;
        for(int i = 1; i <= _rln; ++i)
        {
            if((j = t.findindexn(GET(_leftn,i))) == 0)
            {
                Print(t); Print(*this);
                cerr << "Couldn't find 'left' Index " << GET(_leftn,i) << " in ITensor t.\n";
                Error("operator*(ITensor,Combiner): bad Combiner ITensor product");
            }
            P.from_to(j,t.rn()-_rln+i);
        }

        int k = 1;
        for(int i = 1; i <= t.rn(); ++i)
        if(findindexn(t.index(i)) == 0) 
        {
            P.from_to(i,k++);
            nindices.push_back(t.index(i));
        }

        nindices.push_back(_right);

        foreach(const Index& L, _left1)
        {
#ifdef NDEBUG
            if(!t.hasindex1(L))
            {
                Print(t); Print(*this);
                cout << "Couldn't find 'left' Index " << L << " in ITensor t.\n";
                Error("operator*(ITensor,Combiner): bad Combiner ITensor product");
            }
#endif
            nindices.push_back(L);
        }

        res = ITensor(nindices,t,P);
    }

}; //class Combiner

#endif
