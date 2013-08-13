//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQTENSOR_H
#define __ITENSOR_IQTENSOR_H
#include "iqindex.h"
#include <list>
#include <map>

class IQTDat;
class IQCombiner;
class IQTSparse;

typedef boost::shared_ptr<IQTDat>
IQTDatPtr;


//
// IQTensor
//

class IQTensor
    {
    public:

    //Constructors --------------------------------------------------

    //Construct Null ITensor, isNull returns true
    IQTensor();

    //Construct rank 0 IQTensor (scalar), value set to val
    explicit 
    IQTensor(Real val);

    //Construct rank 1 IQTensor, set to zero
    explicit 
    IQTensor(const IQIndex& i1);

    //Construct rank 2 IQTensor, set to zero
    IQTensor(const IQIndex& i1,const IQIndex& i2);

    //Construct rank 3 IQTensor, set to zero
    IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3);

    //Construct rank 4 IQTensor, set to zero
    IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
             const IQIndex& i4);

    IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
             const IQIndex& i4,const IQIndex& i5);

    IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
             const IQIndex& i4,const IQIndex& i5,const IQIndex& i6);

    IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
             const IQIndex& i4,const IQIndex& i5,const IQIndex& i6,
             const IQIndex& i7);

    IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
             const IQIndex& i4,const IQIndex& i5,const IQIndex& i6,
             const IQIndex& i7,const IQIndex& i8);

    //Construct IQTensor with IQIndices given by vector iqinds
    explicit 
    IQTensor(std::vector<IQIndex>& iqinds);

    //
    // IQIndexVal IQTensor Constructors
    //
    // Given a set of IQIndexVals
    // iv1 = (I1,n1), iv2 = (I2,n2), iv3 = (I3,n3), ...
    // construct an IQTensor T such that
    // T(I1(n1),I2(n2),I3(n3),...) == 1
    //
    explicit
    IQTensor(const IQIndexVal& iv);

    IQTensor(const IQIndexVal& iv1, const IQIndexVal& iv2);

    IQTensor(const IQIndexVal& iv1, const IQIndexVal& iv2,
             const IQIndexVal& iv3);

    IQTensor(const IQIndexVal& iv1, const IQIndexVal& iv2,
             const IQIndexVal& iv3, const IQIndexVal& iv4);

    //Accessor Methods ------------------------------------------

    //Rank of this IQTensor (number of IQIndices)
    int 
    r() const;

    //true if IQTensor has no blocks
    bool 
    empty() const;

    //true if IQTensor is default constructed
    bool 
    isNull() const;

    bool
    isComplex() const;

    //Returns object containing ITensor blocks
    //The ITensors can be iterated over using a Foreach
    //For example, given an IQTensor T,
    //Foreach(const ITensor& t, T.blocks()) { ... }
    const IQTDat&
    blocks() const { return dat(); }
    
    const IndexSet<IQIndex>& 
    indices() const { return *is_; }


    //----------------------------------------------------
    //IQTensor operators

    //
    // Contracting product
    //
    IQTensor 
    operator*(IQTensor other) const 
        { other *= *this; return other; }

    IQTensor& 
    operator*=(const IQTensor& other);

    //
    // Non-Contracting product
    //
    IQTensor 
    operator/(IQTensor other) const 
        { other /= *this; return other; }

    IQTensor& 
    operator/=(const IQTensor& other);

    //
    // Addition and subtraction
    //
    IQTensor& 
    operator+=(const IQTensor& o);

    IQTensor 
    operator+(const IQTensor& o) const 
        { IQTensor res(*this); res += o; return res; }

    IQTensor&
    operator-=(const IQTensor& o)
        { 
        if(this == &o) { operator*=(0); return *this; }
        IQTensor oth(o);
        oth *= -1;
        return operator+=(oth);
        }

    IQTensor 
    operator-(const IQTensor& o) const 
        { IQTensor res(*this); res -= o; return res; }

    //
    // Multiplication by a scalar
    //
    IQTensor& 
    operator*=(Real fac);

    IQTensor& 
    operator/=(Real fac);

    IQTensor
    operator-() const { IQTensor T(*this); T *= -1; return T; }

    IQTensor& 
    operator*=(const LogNumber& lgnum);

    IQTensor& 
    operator*=(Complex z);

    
    //
    // Multiplication by an IQIndexVal
    //
    IQTensor& 
    operator*=(const IQIndexVal& iv) { return operator*=(IQTensor(iv)); }

    //Convert to ITensor
    ITensor 
    toITensor() const;

    //Automatic conversion to ITensor
    operator 
    ITensor() const { return toITensor(); }

    //Inserts an ITensor block or adds it to
    //existing one if already present and QNs match
    IQTensor& 
    operator+=(const ITensor& block);

    //Like operator+=(ITensor) but
    //demands that the block is zero/absent
    //before inserting
    void insert(const ITensor& block);

    //Non-const element access
    Real& 
    operator()(const IQIndexVal& iv1, 
               const IQIndexVal& iv2 = IQIndexVal::Null(), 
	           const IQIndexVal& iv3 = IQIndexVal::Null(), 
               const IQIndexVal& iv4 = IQIndexVal::Null(), 
	           const IQIndexVal& iv5 = IQIndexVal::Null(), 
               const IQIndexVal& iv6 = IQIndexVal::Null(),
	           const IQIndexVal& iv7 = IQIndexVal::Null(), 
               const IQIndexVal& iv8 = IQIndexVal::Null());

    //const element access
    Real 
    operator()(const IQIndexVal& iv1, 
               const IQIndexVal& iv2 = IQIndexVal::Null(), 
	           const IQIndexVal& iv3 = IQIndexVal::Null(), 
               const IQIndexVal& iv4 = IQIndexVal::Null(), 
	           const IQIndexVal& iv5 = IQIndexVal::Null(), 
               const IQIndexVal& iv6 = IQIndexVal::Null(),
	           const IQIndexVal& iv7 = IQIndexVal::Null(), 
               const IQIndexVal& iv8 = IQIndexVal::Null()) const;

    //Method for specifically requesting const access
    Real 
    at(const IQIndexVal& iv1, 
       const IQIndexVal& iv2 = IQIndexVal::Null(), 
	   const IQIndexVal& iv3 = IQIndexVal::Null(), 
       const IQIndexVal& iv4 = IQIndexVal::Null(), 
	   const IQIndexVal& iv5 = IQIndexVal::Null(), 
       const IQIndexVal& iv6 = IQIndexVal::Null(),
	   const IQIndexVal& iv7 = IQIndexVal::Null(), 
       const IQIndexVal& iv8 = IQIndexVal::Null()) const;

    Real
    toReal() const;

    Complex 
    toComplex() const;


    //----------------------------------------------------
    //IQTensor: prime methods

    IQTensor& 
    noprime(IndexType type = All);

    IQTensor& 
    noprime(const IQIndex& I);

    IQTensor& 
    prime(int inc = 1) { prime(All,inc); return *this; }

    IQTensor& 
    prime(IndexType type, int inc = 1);

    IQTensor& 
    prime(const IQIndex& I, int inc = 1);

    //no need to keep prime level small
    IQTensor& 
    mapprime(int plevold, int plevnew, IndexType type = All);

    //----------------------------------------------------
    //IQTensor index methods


    void
    tieIndices(const boost::array<IQIndex,NMAX>& indices, int nind, const IQIndex& tied);

    void
    tieIndices(const IQIndex& i1, const IQIndex& i2, const IQIndex& tied);

    IQTensor&
    trace(const boost::array<IQIndex,NMAX>& indices, int niqind = -1);

    IQTensor&
    trace(const IQIndex& i1, 
          const IQIndex& i2 = IQIndex::Null(), 
          const IQIndex& i3 = IQIndex::Null(),
          const IQIndex& i4 = IQIndex::Null(),
          const IQIndex& i5 = IQIndex::Null(),
          const IQIndex& i6 = IQIndex::Null(),
          const IQIndex& i7 = IQIndex::Null(),
          const IQIndex& i8 = IQIndex::Null());

    //----------------------------------------------------
    //IQTensor miscellaneous methods

    IQTensor&
    takeRealPart();

    IQTensor&
    takeImagPart();

    Real 
    norm() const;

    LogNumber 
    normLogNum() const;

    Real 
    sumels() const;

    template <typename Callable> 
    IQTensor&
    mapElems(const Callable& f);

    void 
    scaleOutNorm();

    void 
    scaleTo(const LogNumber& newscale);

    void 
    clean(Real min_norm = MIN_CUT);

    void 
    randomize();

    IQTensor& 
    conj();

    void
    swap(IQTensor& other);

    void 
    read(std::istream& s);

    void 
    write(std::ostream& s) const;

    //Typedefs -----------------------------------------------------

    typedef IQIndex 
    IndexT;

    typedef IQIndexVal 
    IndexValT;

    typedef IQCombiner 
    CombinerT;

    typedef IQTSparse
    SparseT;

    typedef std::vector<ITensor>::iterator 
    iten_it;

    typedef std::vector<ITensor>::const_iterator 
    const_iten_it;

    typedef IndexSet<IQIndex>::const_iterator
    const_iqind_it;

    //Deprecated methods

    //Get the jth IQIndex of this ITensor, j = 1,2,..,r()
    //const IQIndex& 
    //index(int j) const;

    //Copy IQTensor, incrementing IQIndices of matching IndexType by 1
    //IQTensor(IndexType type, const IQTensor& other);


    private:

    //Data struct ensures const-correct access
    //to the IQTDat, preventing unnecessary sorting of blocks
    struct Data
        {
        Data();

        Data(const IQTDatPtr& p_);

        //Const access
        const IQTDat&
        operator()() const { return *p; }

        //Non-const access
        IQTDat&
        nc() { return *p; }

        void inline
        solo();

        void
        swap(Data& othr) { p.swap(othr.p); }

        private: IQTDatPtr p;
        };

    /////////////////
    // 
    // Data Members

    boost::shared_ptr<IndexSet<IQIndex> >
    is_;

    Data dat;

    //
    /////////////////

    void 
    soloIndex();

    void 
    solo();

    friend class IQTSparse;

    friend void 
    product(const IQTSparse& S, const IQTensor& T, IQTensor& res);


    }; //class IQTensor

IQTensor inline
operator*(IQTensor T, Real fac) {  T *= fac; return T; }

IQTensor inline
operator*(Real fac, IQTensor T) { T *= fac; return T; }

IQTensor inline
operator/(IQTensor T, Real fac) {  T /= fac; return T; }

IQTensor inline
operator/(Real fac, IQTensor T) { T /= fac; return T; }

IQTensor inline
operator*(IQTensor T, const LogNumber& fac) {  T *= fac; return T; }

IQTensor inline
operator*(const LogNumber& fac, IQTensor T) { T *= fac; return T; }

IQTensor inline
operator*(IQTensor T, Complex fac) {  T *= fac; return T; }

IQTensor inline
operator*(Complex fac, IQTensor T) { T *= fac; return T; }

IQTensor inline
operator*(IQTensor T, const IQIndexVal& iv) { T *= iv; return T; }

IQTensor inline
operator*(const IQIndexVal& iv, const IQTensor& T) { return IQTensor(iv) * T; }


//
// Multiplication by an IndexVal
// Result is an ITensor
//
ITensor inline
operator*(const IQTensor& T, const IndexVal& iv)
    { 
    ITensor res = T.toITensor(); 
    return res *= iv; 
    }

ITensor inline
operator*(const IndexVal& iv, const IQTensor& T) 
    { 
    return ITensor(iv) * T.toITensor(); 
    }




class IQTDat : public boost::noncopyable
    {
    public:

    typedef std::vector<ITensor>
    StorageT;

    typedef StorageT::const_iterator
    const_iterator;

    typedef StorageT::iterator
    iterator;

    //
    // Constructors
    //

    IQTDat();

    explicit 
    IQTDat(const IQTDat& other);

    //
    // Accessors
    //

    const_iterator
    begin() const { return blocks_.begin(); }
    const_iterator
    end() const { return blocks_.end(); }

    iterator
    begin() { return blocks_.begin(); }
    iterator
    end() { return blocks_.end(); }

    bool 
    hasBlock(const IndexSet<Index>& is) const;

    ITensor&
    get(const IndexSet<Index>& is);
    const ITensor&
    get(const IndexSet<Index>& is) const;

    int
    size() const { return blocks_.size(); }

    bool
    empty() const { return blocks_.empty(); }

    void
    clear() { blocks_.clear(); }

    void 
    insert(const ITensor& t);

    void 
    insert_add(const ITensor& t);

    void 
    clean(Real min_norm);

    void
    swap(StorageT& new_itensor);

    //
    // Other Methods
    //

    void
    scaleTo(const LogNumber& newscale);

    void 
    read(std::istream& s);

    void 
    write(std::ostream& s) const;

    static const boost::shared_ptr<IQTDat>& 
    Null();

    //void* operator 
    //new(size_t size) 
    //    throw(std::bad_alloc)
    //    { return allocator().alloc(); }

    //void operator 
    //delete(void* p) 
    //    throw()
    //    { return allocator().dealloc(p); }

    private:

    //////////////
    //
    // Data Members
    //

    StorageT blocks_;

    //
    //////////////

    iterator
    findBlock(const IndexSet<Index>& is)
        {
        return find(blocks_.begin(),blocks_.end(),is);
        }

    const_iterator
    findBlock(const IndexSet<Index>& is) const
        {
        return find(blocks_.begin(),blocks_.end(),is);
        }

    bool
    validBlock(const_iterator it) const 
        { 
        return it != blocks_.end(); 
        }

    //Not copyable with =
    void operator=(const IQTDat&);

    //static DatAllocator<IQTDat>& allocator()
    //    {
    //    static DatAllocator<IQTDat> allocator_;
    //    return allocator_;
    //    };

    }; //class IQTDat

template <typename Callable> 
IQTensor& IQTensor::
mapElems(const Callable& f)
    {
    solo();
    Foreach(ITensor& t, dat.nc()) 
        t.mapElems(f);
    return *this;
    }


//
// Computes the scalar/inner/dot product of two
// real-valued IQTensors.
//
// Equivalent to the IQTensor contraction x * y 
// except the result is a Real number versus a 
// rank 0 IQTensor.
//
// N.B. even if the IQIndex arrows of the two
// arguments don't match, the method will 
// automatically correct them.
//
Real 
Dot(IQTensor x, const IQTensor& y);

//
// Scalar (inner) product of two
// possibly complex ITensors.
//
// Conjugates the first argument, therefore
// equivalent to the contraction conj(x) * y 
// (except it yields two real numbers, re and im,
// instead of a rank 0 IQTensor).
//
Complex 
BraKet(IQTensor x, const IQTensor& y);

//Compute divergence of IQTensor T
//
//If DEBUG defined and all blocks do not have
//the same divergence, throws an exception
//(since IQTensor is not correctly constructed).
QN 
div(const IQTensor& T);

const IQIndex&
findIQInd(const IQTensor& T, const Index& i);

QN
qn(const IQTensor& T, const Index& i);

Arrow
dir(const IQTensor& T, const Index& i);

//Return true if one of the ITensor blocks of
//T uses this Index
bool 
usesIndex(const IQTensor& T, const Index& i);

//Returns true if T is exactly zero.
//
//If passed the option Opt("Fast",true),
//only performs fast operations such as checking
//whether T contains any blocks, but skips computing
//the norm of the blocks.
//This can cause the return value to be true even
//if T is actually zero.
bool
isZero(const IQTensor& T, const OptSet& opts = Global::opts());

std::ostream& 
operator<<(std::ostream & s, const IQTensor &t);


#endif
