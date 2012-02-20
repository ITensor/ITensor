//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __IQCOMBINER_H
#define __IQCOMBINER_H
#include "combiner.h"
#include "condenser.h"

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

    IQCombiner();

    IQCombiner(
	    const IQIndex& l1, const IQIndex& l2 = IQIndex::Null(), 
        const IQIndex& l3 = IQIndex::Null(), const IQIndex& l4 = IQIndex::Null(), 
	    const IQIndex& l5 = IQIndex::Null(), const IQIndex& l6 = IQIndex::Null());

    bool 
    doCondense() const { return do_condense; }
    void 
    doCondense(bool val);

    void 
    reset();

    void 
    addleft(const IQIndex& l); 	// Include another left index

    inline bool 
    isInit() const { return initted; }

    // Initialize after all lefts are there and before being used
    void 
    init(std::string rname = "combined", IndexType type = Link, 
         Arrow dir = Switch, int primelevel = 0) const;
    
    operator IQTensor() const;

    const IQIndex& 
    right() const;

    int 
    findindex(const IQIndex& i) const;

    bool 
    hasindex(const IQIndex& I) const;

    bool 
    hasindex(const Index& I) const;

    int 
    num_left() const { return int(left.size()); }

    void
    doprime(PrimeType pr, int inc = 1);

    friend IQCombiner
    primed(IQCombiner C, int inc = 1);

    void 
    conj();

    friend std::ostream& 
    operator<<(std::ostream & s, const IQCombiner & c);

    IQTensor 
    operator*(const IQTensor& t) const { IQTensor res; product(t,res); return res; }

    friend inline IQTensor 
    operator*(const IQTensor& t, const IQCombiner& c) { return c.operator*(t); }

    void 
    product(IQTensor t, IQTensor& res) const;

    private:

    /////////////
    //
    // Data Members
    //

    std::vector<IQIndex> left;
    mutable IQIndex right_;
    mutable std::vector<Combiner> combs;
    //mutable std::map<ApproxReal, Combiner> setcomb;
    //mutable std::map<Index, Combiner> rightcomb;
    mutable bool initted;

    mutable Condenser cond;
    mutable IQIndex ucright_;
    bool do_condense;

    //
    /////////////

    typedef std::map<ApproxReal, Combiner>::iterator
    setcomb_it;

    typedef std::map<Index, Combiner>::iterator
    rightcomb_it;

    };

class QCounter
    {
    public:

    QCounter(const std::vector<IQIndex>& v)
        : don(false)
        {
        Foreach(const IQIndex& I,v)
            {
            n.push_back(I.nindex());
            ind.push_back(0);
            }
        }

    bool 
    notdone() const { return !don; }

    QCounter& 
    operator++()
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

    void 
    getVecInd(const std::vector<IQIndex>& v, std::vector<Index>& vind, QN& q) const
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

    private:

    bool don;

    std::vector<int> n,
                     ind;
    };

inline IQCombiner::
IQCombiner() 
    : initted(false), 
      do_condense(false) 
    { }

inline IQCombiner::
IQCombiner(
    const IQIndex& l1, const IQIndex& l2, 
    const IQIndex& l3, const IQIndex& l4, 
    const IQIndex& l5, const IQIndex& l6)
    : initted(false), 
      do_condense(false)
    {
    if(l1 == IQIndex::Null()) Error("Null IQIndex");
    if(l1 != IQIndex::Null()) left.push_back(l1); 
    if(l2 != IQIndex::Null()) left.push_back(l2);
    if(l3 != IQIndex::Null()) left.push_back(l3); 
    if(l4 != IQIndex::Null()) left.push_back(l4);
    if(l5 != IQIndex::Null()) left.push_back(l5); 
    if(l6 != IQIndex::Null()) left.push_back(l6);
    Foreach(IQIndex& L, left) L.conj();
    }

inline
void IQCombiner::
doCondense(bool val)
    {
    if(initted) 
        Error("Can't set doCondense after already initted.");
    do_condense = val;
    }

inline 
void IQCombiner::
reset()
    {
    left.clear();
    //setcomb.clear();
    //rightcomb.clear();
    initted = false;
    }

inline
void IQCombiner::
addleft(const IQIndex& l) 	// Include another left index
	{ 
    if(l == IQIndex::Null()) Error("Null IQIndex");
    left.push_back(l);
    //Flip arrows to make combiner compatible with
    //the IQTensor from which it got its left indices
    left.back().conj();
    initted = false;
	}

inline
void IQCombiner::
init(std::string rname, IndexType type, 
     Arrow dir, int primelevel) const 
    {
    if(initted) return;
    if(left.size() == 0)
        Error("No left indices in IQCombiner.");

    Arrow rdir; 
    if(dir == Switch) //determine automatically
        {
        rdir = Switch*left.back().dir();

        //Prefer to derive right Arrow from Link indices
        for(size_t j = 0; j < left.size(); ++j)
        if(left[j].type() == Link) 
            { 
            rdir = Switch*left[j].dir(); 
            break;
            }
        }
    else
        { rdir = dir; }

    //setcomb.clear();
    //rightcomb.clear();

    //Construct individual Combiners
    QCounter c(left);
    std::vector<inqn> iq;
    for( ; c.notdone(); ++c)
        {
        std::vector<Index> vind;
        QN q;
        c.getVecInd(left, vind, q);		// updates vind and q
        q *= -rdir;

        combs.push_back(Combiner());
        Combiner& co = combs.back();
        Foreach(const Index& i, vind)
            { 
            co.addleft(i); 
            //rss += i.uniqueReal(); 
            }
        co.init(rname+q.toString(),type,rdir,primelevel);

        iq.push_back(inqn(co.right(),q));
        //setcomb[ApproxReal(rss)] = co;
        //rightcomb[co.right()] = co;
        }
    if(do_condense) 
        {
        ucright_ = IQIndex(rname,iq,rdir,primelevel);
        std::string cname = "cond::" + rname;
        cond = Condenser(ucright_,cname);
        right_ = cond.smallind();
        }
    else 
        {
        right_ = IQIndex(rname,iq,rdir,primelevel);
        }
    initted = true;
	}

inline IQCombiner::
operator IQTensor() const
    {
    if(!initted) Error("IQCombiner::operator IQTensor(): IQCombiner not initialized.");

    //if(right_.m() > 16) 
    //{ std::cerr << std::endl << std::endl << "WARNING: too large of an m in IQCombiner::operator IQTensor(). May be inefficient!" << std::endl << std::endl; }

    std::vector<IQIndex> iqinds(left);
    iqinds.push_back((do_condense ? ucright_ : right_));
    IQTensor res(iqinds);
    /*
    for(std::map<ApproxReal,Combiner>::const_iterator it = setcomb.begin();
        it != setcomb.end(); 
        ++it)
        { res.insert(it->second); }
    */
    Foreach(const Combiner& co, combs)
        {
        //Here we are using the fact that Combiners
        //can be converted to ITensors
        res.insert(co);
        }

    //Combiners should always have the 
    //structure of zero divergence IQTensors
    DO_IF_DEBUG(checkQNs(res));

    if(do_condense) 
        { 
        IQTensor rcopy(res); 
        cond.product(rcopy,res); 
        }

    return res;
    }

inline
const IQIndex& IQCombiner::
right() const 
    { 
    init();
    return right_;
    }

inline
int IQCombiner::
findindex(const IQIndex& i) const
	{
    for(size_t j = 0; j < left.size(); ++j)
    if(left[j] == i) return j;
    return -1;
	}

inline
bool IQCombiner::
hasindex(const IQIndex& I) const
	{
    for(size_t j = 0; j < left.size(); ++j)
        if(left[j] == I) return true;
    return false;
	}

inline
bool IQCombiner::
hasindex(const Index& i) const
    {
    for(size_t j = 0; j < left.size(); ++j)
        if(left[j].hasindex(i)) return true;
    return false;
    }

void inline IQCombiner::
doprime(PrimeType pr, int inc)
    {
    Foreach(IQIndex& ll, left)
        ll.doprime(pr,inc);
    Foreach(Combiner& co, combs)
        co.doprime(pr,inc);
    if(initted)
        {
        right_.doprime(pr,inc);
        if(do_condense) 
            {
            cond.doprime(pr,inc);
            ucright_.doprime(pr,inc);
            }
        }
    }

IQCombiner inline
primed(IQCombiner C, int inc)
    {
    C.doprime(primeBoth,inc);
    return C;
    }


inline
void IQCombiner::
conj() 
    { 
    init();
    Foreach(IQIndex& I, left) I.conj(); 
    if(do_condense) 
        {
        cond.conj();
        ucright_.conj();
        }
    right_.conj();
    }

inline 
std::ostream& 
operator<<(std::ostream & s, const IQCombiner & c)
    {
    if(c.isInit())
        { s << std::endl << "Right index is " << c.right() << "\n"; }
    else
        { s << std::endl << "Right index is not initialized\n\n"; }
    s << "Left indices: \n";
    Foreach(const IQIndex& I, c.left) s << I << std::endl;
    return s << "\n\n";
    }

inline
void IQCombiner::
product(IQTensor t, IQTensor& res) const
    {
    init();
    std::vector<IQIndex> iqinds;

    int j;
    //t has right IQIndex, expand it
    if((j = t.findindex(right_)) != 0)
        {
        IQTensor t_uncondensed;
        if(do_condense) 
            { 
            cond.product(t,t_uncondensed); 
            j = t_uncondensed.findindex(ucright_);
            }
        const IQTensor& t_ = (do_condense ? t_uncondensed : t);
        const IQIndex& r = (do_condense ? ucright_ : right_);

        if(Globals::checkArrows())
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

        std::map<Index, const Combiner*> rightcomb;
        Foreach(const Combiner& co, combs)
            {
            rightcomb[co.right()] = &co;
            }

        Foreach(const ITensor& tt, t_.itensors())
            {
            for(int k = 1; k <= tt.r(); ++k)
                {
                if(r.hasindex(tt.index(k)))
                    { 
                    res += (*(rightcomb[tt.index(k)]) * tt); 
                    break;
                    }
                } //end for
            } //end Foreach

        }
    else
        {
        //t has left IQIndex's, combine them

        //res will have all IQIndex's of t not in the left of c
        Foreach(const IQIndex& I, t.iqinds()) 
            { 
            if(!hasindex(I)) iqinds.push_back(I); 
            }
        //and res will have c's right IQIndex
        if(do_condense) iqinds.push_back(ucright_);
        else            iqinds.push_back(right_);

        res = IQTensor(iqinds);

        //Check left indices
        Foreach(const IQIndex& I, left)
            {
            if((j = t.findindex(I)) == 0)
                {
                std::cerr << "Could not find left IQIndex " << I << "\n";
                t.printIndices("t");
                std::cerr << "Left indices\n";
                for(size_t j = 0; j < left.size(); ++j)
                    { 
                    std::cerr << j SP left[j] << "\n"; 
                    }
                Error("bad IQCombiner IQTensor product");
                }
            else //IQIndex is in left
                {
                //Check arrow directions
                if(Globals::checkArrows())
                    if(t.index(j).dir() == I.dir())
                        {
                        std::cerr << "IQTensor = " << t << std::endl;
                        std::cerr << "IQCombiner = " << *this << std::endl;
                        std::cerr << "IQIndex from IQTensor = " << t.index(j) << std::endl;
                        std::cerr << "(Left) IQIndex from IQCombiner = " << I << std::endl;
                        Error("Incompatible arrow directions in operator*(IQTensor,IQCombiner).");
                        }
                }
            }

        std::map<ApproxReal, const Combiner*> setcomb;
        typedef std::map<ApproxReal, const Combiner*>::const_iterator
        setcomb_const_it;

        Foreach(const Combiner& co, combs)
            {
            setcomb[co.uniqueReal()] = &co;
            }

        for(IQTensor::const_iten_it i = t.const_iten_begin(); i != t.const_iten_end(); ++i)
            {
            Real rse = 0;
            for(int k = 1; k <= i->r(); ++k)
                {
                if(this->hasindex(i->index(k))) 
                    rse += i->index(k).uniqueReal();
                }

            if(setcomb.count(rse) == 0)
                {
                Print(*i);
                std::cerr << "\nleft indices \n";
                for(size_t j = 0; j < left.size(); ++j)
                    { std::cerr << j << " " << left[j] << "\n"; }
                std::cerr << "\n\n";
                for(setcomb_const_it uu = setcomb.begin();
                    uu != setcomb.end(); ++uu)
                    {
                    std::cout << "Combiner: " << std::endl;
                    std::cout << *(uu->second) << std::endl;
                    }
                Error("no setcomb for rse in IQCombiner prod");
                }

            res += (*setcomb[rse] * (*i));

            }

        if(do_condense) 
            { 
            IQTensor rcopy(res); 
            cond.product(rcopy,res); 
            }
        }
    } //void product(const IQTensor& t, IQTensor& res) const


#endif
