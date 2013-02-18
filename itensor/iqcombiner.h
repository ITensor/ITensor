//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __IQCOMBINER_H
#define __IQCOMBINER_H
#include "combiner.h"
#include "condenser.h"
#include "qcounter.h"

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

    typedef std::vector<IQIndex>::const_iterator 
    left_it;

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

    const std::pair<left_it,left_it> 
    left() const 
        { 
        return std::make_pair(left_.begin(),left_.end()); 
        }

    inline bool 
    isInit() const { return initted; }

    // Initialize after all lefts are there and before being used
    void 
    init(std::string rname = "combined", IndexType type = Link, 
         Arrow dir = Switch, int primelevel = 0) const;
    
    operator IQTensor() const;

    const IQIndex& 
    right() const;

    bool 
    hasindex(const IQIndex& I) const;

    bool 
    hasindex(const Index& I) const;

    int 
    num_left() const { return int(left_.size()); }

    void
    prime(int inc = 1) { prime(All,inc); }

    void
    prime(IndexType type, int inc = 1);

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

    std::vector<IQIndex> left_;
    mutable IQIndex right_;
    mutable std::vector<Combiner> combs;
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
    if(l1 != IQIndex::Null()) left_.push_back(l1); 
    if(l2 != IQIndex::Null()) left_.push_back(l2);
    if(l3 != IQIndex::Null()) left_.push_back(l3); 
    if(l4 != IQIndex::Null()) left_.push_back(l4);
    if(l5 != IQIndex::Null()) left_.push_back(l5); 
    if(l6 != IQIndex::Null()) left_.push_back(l6);
    Foreach(IQIndex& L, left_) L.conj();
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
    left_.clear();
    initted = false;
    }

inline
void IQCombiner::
addleft(const IQIndex& l) 	// Include another left index
	{ 
    if(l == IQIndex::Null()) Error("Null IQIndex");
    left_.push_back(l);
    //Flip arrows to make combiner compatible with
    //the IQTensor from which it got its left indices
    left_.back().conj();
    initted = false;
	}

inline
void IQCombiner::
init(std::string rname, IndexType type, 
     Arrow dir, int primelevel) const 
    {
    if(initted) return;
    if(left_.size() == 0)
        Error("No left indices in IQCombiner.");

    Arrow rdir; 
    if(dir == Switch) //determine automatically
        {
        rdir = Switch*left_.back().dir();

        //Prefer to derive right Arrow from Link indices
        for(size_t j = 0; j < left_.size(); ++j)
        if(left_[j].type() == Link) 
            { 
            rdir = Switch*left_[j].dir(); 
            break;
            }
        }
    else
        { rdir = dir; }

    //Construct individual Combiners
    QCounter c(left_);
    std::vector<inqn> iq;
    for( ; c.notDone(); ++c)
        {
        std::vector<Index> vind;
        QN q;
        c.getVecInd(left_,vind, q); // updates vind and q
        q *= -rdir;

        combs.push_back(Combiner());
        Combiner& co = combs.back();
        co.addleft(vind);
        co.init(rname+q.toString(),type,rdir,primelevel);

        iq.push_back(inqn(co.right(),q));
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

    std::vector<IQIndex> iqinds(left_);
    iqinds.push_back((do_condense ? ucright_ : right_));
    IQTensor res(iqinds);

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
bool IQCombiner::
hasindex(const IQIndex& I) const
	{
    for(size_t j = 0; j < left_.size(); ++j)
        if(left_[j] == I) return true;
    return false;
	}

inline
bool IQCombiner::
hasindex(const Index& i) const
    {
    for(size_t j = 0; j < left_.size(); ++j)
        if(left_[j].hasindex(i)) return true;
    return false;
    }

void inline IQCombiner::
prime(IndexType type, int inc)
    {
    Foreach(IQIndex& ll, left_)
        ll.prime(type,inc);
    Foreach(Combiner& co, combs)
        co.prime(type,inc);
    if(initted)
        {
        right_.prime(type,inc);
        if(do_condense) 
            {
            cond.prime(type,inc);
            ucright_.prime(type,inc);
            }
        }
    }

IQCombiner inline
primed(IQCombiner C, int inc)
    {
    C.prime(All,inc);
    return C;
    }


inline
void IQCombiner::
conj() 
    { 
    init();
    Foreach(IQIndex& I, left_) I.conj(); 
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
    Foreach(const IQIndex& I, c.left_) s << I << std::endl;
    return s << "\n\n";
    }

inline
void IQCombiner::
product(IQTensor T, IQTensor& res) const
    {
    init();
    std::vector<IQIndex> iqinds;

    if(T.hasindex(right_))
        {
        //
        //T has right IQIndex, expand it
        //
        IQTensor T_uncondensed;
        if(do_condense) 
            { 
            cond.product(T,T_uncondensed); 
            }
        const IQTensor& T_ = (do_condense ? T_uncondensed : T);
        const IQIndex& r = (do_condense ? ucright_ : right_);

        if(Global::checkArrows())
            if(dir(T_.indices(),r) == r.dir())
                {
                std::cerr << "IQTensor = " << T_ << std::endl;
                std::cerr << "IQCombiner = " << *this << std::endl;
                std::cerr << "(Right) IQIndex from IQCombiner = " << r << std::endl;
                Error("Incompatible arrow directions in operator*(IQTensor,IQCombiner).");
                }

        iqinds.reserve(T_.indices().r()-1+left_.size());

        Foreach(const IQIndex& I, T_.indices())
            {
            if(I == r)
                copy(left_.begin(),left_.end(),std::back_inserter(iqinds));
            else
                iqinds.push_back(I);
            }

        res = IQTensor(iqinds);

        std::map<Index, const Combiner*> rightcomb;
        Foreach(const Combiner& co, combs)
            {
            rightcomb[co.right()] = &co;
            }

        Foreach(const ITensor& tt, T_.blocks())
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
        //
        //T has left IQIndex's, combine them
        //

        iqinds.reserve(T.r()-left_.size()+1);

        //res will have all IQIndex's of T not in the left of c
        Foreach(const IQIndex& I, T.indices()) 
            { 
            if(!hasindex(I)) iqinds.push_back(I); 
            }
        //and res will have c's right IQIndex
        if(do_condense) iqinds.push_back(ucright_);
        else            iqinds.push_back(right_);

        res = IQTensor(iqinds);

        //Check left indices
        Foreach(const IQIndex& I, left_)
            {
            if(!T.hasindex(I))
                {
                std::cerr << "Could not find left IQIndex " << I << "\n";
                T.printIndices("T");
                std::cerr << "Left indices\n";
                for(size_t j = 0; j < left_.size(); ++j)
                    { 
                    std::cerr << j SP left_[j] << "\n"; 
                    }
                Error("bad IQCombiner IQTensor product");
                }
            else //IQIndex is in left
                {
                //Check arrow directions
                if(Global::checkArrows())
                    if(dir(T.indices(),I) == I.dir())
                        {
                        PrintIndices(T);
                        Print((*this));
                        std::cerr << "(Left) IQIndex from IQCombiner = " << I << std::endl;
                        Error("Incompatible arrow directions in operator*(IQTensor,IQCombiner).");
                        }
                }
            }

        //Create map of Combiners using uniqueReal as key
        std::map<ApproxReal, const Combiner*> combmap;
        Foreach(const Combiner& co, combs)
            {
            combmap[co.uniqueReal()] = &co;
            }

        //Loop over each block in T and apply appropriate
        //Combiner (determined by the uniqueReal of the 
        //combined Indices)
        Foreach(const ITensor& t, T.blocks())
            {
            Real block_ur = 0;
            for(int k = 1; k <= t.r(); ++k)
                {
                if(this->hasindex(t.index(k))) 
                    block_ur += t.index(k).uniqueReal();
                }

            if(combmap.count(block_ur) == 0)
                {
                Print(t);
                std::cerr << "\nleft indices \n";
                for(size_t j = 0; j < left_.size(); ++j)
                    { std::cerr << j << " " << left_[j] << "\n"; }
                std::cerr << "\n\n";

                typedef std::map<ApproxReal, const Combiner*>::const_iterator
                combmap_const_it;
                for(combmap_const_it uu = combmap.begin();
                    uu != combmap.end(); ++uu)
                    {
                    std::cout << "Combiner: " << std::endl;
                    std::cout << *(uu->second) << std::endl;
                    }
                Error("no combmap entry for block_ur in IQCombiner prod");
                }

            res += (*combmap[block_ur] * t);
            }

        if(do_condense) 
            { 
            IQTensor rcopy(res); 
            cond.product(rcopy,res); 
            }
        }
    } //void product(const IQTensor& T, IQTensor& res) const


#endif
