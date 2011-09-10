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
    boost::array<Index,NMAX+1> left_; // max dim is 8
    mutable Index right_;
    int rl_; //Number of m>1 'left' indices (indices to be combined into one)
    mutable bool initted;
public:

    //Accessor Methods ----------------------------------------------

    const Index& right() const 
    { 
        if(!initted) Error("Combiner: right ind requested prior to init");
        return right_; 
    }
    int rl() const { return rl_; }
    const Index& left(int j) const { return GET(left_,j); }

    typedef boost::array<Index,NMAX+1>::const_iterator left_it;
    const pair<left_it,left_it> left() const 
    { 
        return make_pair(left_.begin()+1,left_.begin()+rl_+1); 
    }

    //Constructors --------------------------------------------------

    Combiner() : rl_(0), initted(false) {}
    Combiner(const Index& l1 , const Index& l2 = IndNull, 
    const Index& l3 = IndNull, const Index& l4 = IndNull, 
    const Index& l5 = IndNull, const Index& l6 = IndNull, 
    const Index& l7 = IndNull, const Index& l8 = IndNull )
    : rl_(0), initted(false)
	{
        boost::array<const Index*,NMAX+1> ll 
        = {{ &IndNull, &l1, &l2, &l3, &l4, &l5, &l6, &l7, &l8 }};

        do { ++rl_; left_[rl_] = *ll[rl_]; } 
        while(rl_ < NMAX && *ll[rl_+1] != IndNull);

        assert(rl_ == NMAX || left_[rl_+1] == IndNull);
        assert(left_[rl_] != IndNull);
	}
    
    //Operators -----------------------------------------------------

    ITensor operator*(const ITensor& t) const 
        { ITensor res; product(t,res); return res; }
    friend inline ITensor operator*(const ITensor& t, const Combiner& c) 
        { return c.operator*(t); }

    //Index Methods -------------------------------------------------

    inline bool check_init() const { return initted; }

    void addleft(const Index& l)// Include another left index
	{ 
        initted = false;
        if(rl_ == NMAX) 
            Error("Combiner: already reached max number of left indices.");
        left_[++rl_] = l; 
	}

    //Initialize after all lefts are added and before being used
    void init(string rname = "combined", IndexType type = Link, 
    int primelevel = 0) const
	{
        if(initted) return;
        int m = 1; for(int i = 1; i <= rl_; ++i) { m *= left_[i].m(); }
        right_ = Index(rname,m,type,primelevel); 
        initted = true;
    }

    int findindex(const Index& I) const
	{
        for(int j = 1; j <= rl_; ++j)
            { if(left_[j] == I) return j; }
        return 0;
	}
    bool hasindex(const Index& I) const
	{
        for(int j = 1; j <= rl_; ++j) if(left_[j] == I) return true;
        return false;
	}

    //Other Methods -------------------------------------------------

    operator ITensor() const
    {
        if(right_.m() > 16) 
        { 
            cerr << "\n\n" 
            << "WARNING: too large of an m in Combiner to ITensor!\n\n"; 
        }

        //Use a kronecker delta tensor to convert this Combiner into an Tensor
        ITensor res = operator*(ITensor(right_,right_.primed(),1));
        res.noprimeind(right_.primed());
        return res;
    }

    friend inline ostream & operator<<(ostream & s, const Combiner & c)
    {
        s << "\nRight index: " << c.right() << "\n";
        s << "Left indices:\n";
        foreach(const Index& l, c.left()) s << " " << l << "\n";
        return s;
    }

    void conj() { }

    void product(const ITensor& t, ITensor& res) const
    {
        init();

        int j;
        if((j = t.findindex(right_)) != 0)
        {
            vector<Index> nindices; 
            nindices.reserve(t.r()+rl_-1);
            for(int i = 1; i < j; ++i)
                nindices.push_back(t.index(i));
            for(int i = 1; i <= rl_; ++i)
                nindices.push_back(left_[i]);
            for(int i = j+1; i <= t.r(); ++i)
                nindices.push_back(t.index(i));
            res = ITensor(nindices,t);
            return;
        }

        t.groupIndices(left_,rl_,right_,res);
    }

    //For interface compatibility with IQCombiner
    void doCondense(bool val) { } 

}; //class Combiner

#endif
