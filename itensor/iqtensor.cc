//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "iqtensor.h"

namespace itensor {

using std::istream;
using std::ostream;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::find;
using std::pair;
using std::make_pair;
using std::make_shared;
using std::move;

IQTensor::
IQTensor(Complex val) 
    :
    scale_(1.)
    { 
    if(val.imag()==0)
        store_ = make_shared<ITDiag<Real>>(val.real());
    else
        store_ = make_shared<ITDiag<Complex>>(val);
    }

IQTensor::
IQTensor(const QN& q, vector<IQIndex>&& iqinds) 
	: 
    is_(move(iqinds)),
    scale_(1.),
    store_(make_shared<IQTDense<Real>>(is_,q))
	{ }

IQTensor& IQTensor::
operator+=(const IQTensor& other)
    {
    Error("Not implemented");

    //
    // Idea to implement:
    // - Compute permutation P from is_ to other.is_
    // - If P trivial, do a straight daxpy on the data
    // - If P not trivial, loop over blocks of *this
    //   and get pointer to corresponding block of other 
    //   by applying P and add those blocks.
    //   (Instead of looping over a GCounter
    //   to visit all non-zero blocks, may be faster
    //   to iterate through offset_ and calculate block indices
    //   from position of non-negative offset_ elements.)
    //

    return *this;
    }

class QContract
    {
    IQIndexSet Nis_;
    public:

    QContract() { }

    operator IQIndexSet() const { return move(Nis_); }
    };

IQTensor& IQTensor::
operator*=(const IQTensor& other)
    {
    Error("Not implemented");

    //
    // Ideas to implement:
    //
    // o Multithread the block-block contractions
    //

    return *this;
    }



IQIndex
findIQInd(const IQTensor& T, const Index& i)
    {
    for(const IQIndex& J : T.indices())
        if(hasindex(J,i)) return J;
    Print(T.indices());
    Print(i);
    throw ITError("Index i not found in any of T's IQIndices");
    return IQIndex();
    }


Arrow
dir(const IQTensor& T, const IQIndex& I)
	{
    for(const IQIndex& J : T.indices())
        if(I == J) return J.dir();
    Error("dir: IQIndex not found");
    return Out;
	}

bool
isZero(const IQTensor& T, const Args& args)
    {
    Error("Not implemented");
    //if(T.empty()) return true;
    ////done with all fast checks
    //if(args.getBool("Fast",false)) return false;
    //for(const ITensor& t : T.blocks())
    //    {
    //    if(!isZero(t)) return false;
    //    }
    return true;
    }


}; //namespace itensor
