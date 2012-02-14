#ifndef __ITENSOR_SPARSEITENSOR_H
#define __ITENSOR_SPARSEITENSOR_H
#include "itensor.h"

//
// SparseITensor
//
class SparseITensor : protected ITensor
    {
    public:

    //
    // Constructors
    //

    SparseITensor();

    SparseITensor(const Index& i1);

    SparseITensor(const Index& i1, Real d);

    SparseITensor(const Index& i1, const Vector& diag);

    SparseITensor(const Index& i1, const Index& i2, Real d);

    //
    // Accessor Methods
    //

    using ITensor::index;
    using ITensor::r;
    using ITensor::m;
    using ITensor::isComplex;
    using ITensor::isNotComplex;
    using ITensor::scale;

    //
    // Operators
    //

    ITensor
    operator*(const ITensor& T) const
        {
        ITensor res;
        product(*this,T,res);
        return res;
        }

    friend inline ITensor&
    operator*=(ITensor& T, const SparseITensor& S)
        {
        ITensor res;
        product(S,T,res);
        T.swap(res);
        return T;
        }

    ITensor friend inline
    operator*(const ITensor& T, const SparseITensor& S)
        { 
        ITensor res; 
        product(S,T,res);
        return res; 
        }

    //
    // Index Methods
    //

    using ITensor::findtype;
    using ITensor::findindex;
    using ITensor::hasindex;
    using ITensor::hasAllIndex;
    using ITensor::notin;
    using ITensor::mapindex;

    //
    // Primelevel Methods
    //

    using ITensor::noprime;
    using ITensor::doprime;
    using ITensor::primeall;
    using ITensor::primesite;
    using ITensor::primelink;
    using ITensor::mapprime;
    using ITensor::primeind;
    using ITensor::noprimeind;

    //
    // Other Methods
    //

    using ITensor::print;
    using ITensor::printIndices;

    typedef Index 
    IndexT;

    typedef IndexVal 
    IndexValT;

    typedef Combiner 
    CombinerT;

    private:

    //////////////
    //
    // Data members
    //

    //diagonal elements
    Vector diag_;

    //boost::array<Index,NMAX+1>
    using ITensor::index_;

    //int
    using ITensor::r_;
    using ITensor::rn_;

    //mutable LogNumber
    using ITensor::scale_;

    //Real
    using ITensor::ur;

    //
    //////////////

    using ITensor::setUniqueReal;

    void
    _construct1(const Index& i1);

    void
    _construct2(const Index& i1, const Index& i2);

    void friend
    product(const SparseITensor& S, const ITensor& T, ITensor& res);

    }; // class SparseITensor

#endif
