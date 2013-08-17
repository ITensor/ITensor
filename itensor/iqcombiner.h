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

    IQCombiner(const IQCombiner& other);

    bool 
    doCondense() const { return do_condense; }
    void 
    doCondense(bool val);

    void 
    reset();

    void 
    addleft(const IQIndex& l); 	// Include another left index

    const std::vector<IQIndex>&
    left() const { return left_; }

    inline bool 
    isInit() const { return initted; }

    // Initialize after all lefts are there and before being used
    void 
    init(std::string rname = "combined", 
         IndexType type = Link, 
         Arrow dir = Neither, 
         int primelevel = 0) const;
    
    operator IQTensor() const;

    const IQIndex& 
    right() const;

    int 
    num_left() const { return int(left_.size()); }

    void
    prime(int inc = 1) { prime(All,inc); }

    void
    prime(IndexType type, int inc = 1);

    void 
    conj();

    IQTensor 
    operator*(const IQTensor& t) const { IQTensor res; product(t,res); return res; }

    IQCombiner&
    operator=(const IQCombiner& other);

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



//
// IQCombiner helper methods
//

IQTensor inline
operator*(const IQTensor& t, const IQCombiner& c) { return c.operator*(t); }

bool
hasindex(const IQCombiner& C, const IQIndex& I);

bool
hasindex(const IQCombiner& C, const Index& i);

std::ostream& 
operator<<(std::ostream & s, const IQCombiner & c);


#endif
