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

    //Typedefs -----------------------------------------------------

    typedef std::list<ITensor>::iterator 
    iten_it;

    typedef std::list<ITensor>::const_iterator 
    const_iten_it;

    typedef IndexSet<IQIndex>::const_iterator
    const_iqind_it;

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

    //Copy IQTensor, incrementing IQIndices of matching IndexType by 1
    IQTensor(IndexType type, const IQTensor& other);

    explicit 
    IQTensor(std::istream& s);

    static 
    const IQTensor& 
    Complex_1();

    static 
    const IQTensor& 
    Complex_i();

    static 
    const IQTensor& 
    ComplexProd();

    //Accessor Methods ------------------------------------------

    //Rank of this IQTensor (number of IQIndices)
    int 
    r() const;

    //Get the jth IQIndex of this ITensor, j = 1,2,..,r()
    const IQIndex& 
    index(int j) const;

    //Number of ITensor blocks
    int 
    iten_size() const;

    bool 
    iten_empty() const;

    //true if IQTensor is default constructed
    bool 
    isNull() const;

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

    IQTensor 
    operator*(Real fac) const
        { IQTensor res(*this); res *= fac; return res; }

    friend inline IQTensor 
    operator*(Real fac, IQTensor T) 
        { T *= fac; return T; }

    IQTensor& 
    operator/=(Real fac);

    IQTensor 
    operator/(Real fac) const 
        { IQTensor res(*this); res /= fac; return res; }

    friend inline IQTensor 
    operator/(Real fac, IQTensor t) 
        { return (t /= fac); }

    IQTensor
    operator-() const { IQTensor T(*this); T *= -1; return T; }

    IQTensor& 
    operator*=(const LogNumber& lgnum);

    IQTensor 
    operator*(const LogNumber& lgnum) const
        { IQTensor res(*this); res *= lgnum; return res; }

    friend inline IQTensor 
    operator*(const LogNumber& lgnum, IQTensor T) 
        { T *= lgnum; return T; }

    //
    // Contracting product with an ITensor
    // Result is an ITensor
    //
    ITensor 
    operator*(const ITensor& t) const
        { ITensor res = this->toITensor(); res *= t; return res; }
    
    //
    // Multiplication by an IQIndexVal
    //
    IQTensor& 
    operator*=(const IQIndexVal& iv)
        { (*this) *= IQTensor(iv); return *this; }

    IQTensor 
    operator*(const IQIndexVal& iv) const
        { IQTensor res(*this); res *= iv; return res; }

    friend inline IQTensor 
    operator*(const IQIndexVal& iv, const IQTensor& T) 
        { return IQTensor(iv) * T; }

    //
    // Multiplication by an IndexVal
    //
    ITensor 
    operator*(const IndexVal& iv) const
        { ITensor res = this->toITensor(); res *= iv; return res; }

    friend inline ITensor 
    operator*(const IndexVal& iv, const IQTensor& T) 
        { return ITensor(iv) * T.toITensor(); }

    //Convert to ITensor
    ITensor 
    toITensor() const;

    //Automatic conversion to ITensor
    operator 
    ITensor() const 
        { return toITensor(); }

    //Inserts an ITensor block or adds it to
    //existing one if already present and QNs match
    IQTensor& 
    operator+=(const ITensor& t);

    //Like operator+=(ITensor) but
    //demands that the block is zero/absent
    //before inserting
    void insert(const ITensor& t);

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


    //----------------------------------------------------
    //IQTensor: prime methods

    void 
    noprime(IndexType type = All);

    void 
    noprime(const IQIndex& I);

    void 
    prime(int inc = 1) { prime(All,inc); }

    void 
    prime(IndexType type, int inc = 1);

    void 
    prime(const IQIndex& I, int inc = 1);

    //no need to keep prime level small
    void 
    mapprime(int plevold, int plevnew, IndexType type = All);

    //----------------------------------------------------
    //IQTensor index methods

    //Return true if one of the ITensors uses this Index
    bool 
    uses_ind(const Index& i) const;

    void
    tieIndices(const boost::array<IQIndex,NMAX>& indices, int nind, const IQIndex& tied);

    void
    tieIndices(const IQIndex& i1, const IQIndex& i2, const IQIndex& tied);

    void
    trace(const boost::array<IQIndex,NMAX>& indices, int niqind = -1);

    void
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

    Real 
    uniqueReal() const;

    Real 
    norm() const;

    LogNumber 
    normLogNum() const;

    Real 
    sumels() const;

    template <typename Callable> void
    mapElems(const Callable& f);

    void 
    scaleOutNorm();

    void 
    scaleTo(const LogNumber& newscale);

    void 
    clean(Real min_norm = MIN_CUT);

    int 
    vecSize() const;

    int 
    maxSize() const;

    void 
    assignToVec(VectorRef v) const;

    void 
    assignFromVec(VectorRef v);

    Real
    toReal() const;

    void 
    toComplex(Real& re, Real& im) const;
    
    void 
    randomize();

    void 
    conj();

    void
    swap(IQTensor& other);

    typedef IQIndex 
    IndexT;

    typedef IQIndexVal 
    IndexValT;

    typedef IQCombiner 
    CombinerT;

    typedef IQTSparse
    SparseT;

    friend class IQTSparse;

    friend void 
    product(const IQTSparse& S, const IQTensor& T, IQTensor& res);

    void 
    read(std::istream& s);

    void 
    write(std::ostream& s) const;

    static 
    const IQIndex& 
    ReImIndex() { return IQIndex::IndReIm(); }

    //Deprecated methods --------------------------

    //Use realPart(T) and imagPart(T) instead
    //void 
    //splitReIm(IQTensor& re, IQTensor& im) const;

    //
    //Renamed to randomize in keeping with code conventions
    //
    //void 
    //Randomize();

    //
    //Renamed to randomize in keeping with code conventions
    //
    //void 
    //SplitReIm(IQTensor& re, IQTensor& im) const;

    //Use prime(I) instead
    //void 
    //ind_inc_prime(const IQIndex& I,int inc);

    //Use primed(A,Site) instead
    //friend inline IQTensor 
    //primesite(IQTensor A, int inc = 1) { A.prime(Site,inc); return A; }

    //Use primed(A,Link) instead
    //friend inline IQTensor 
    //primelink(IQTensor A, int inc = 1) { A.prime(Link,inc); return A; }

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

    }; //class IQTensor

inline
const IQTensor&
IQComplex_1() { return IQTensor::Complex_1(); }

inline
const IQTensor&
IQComplex_i() { return IQTensor::Complex_i(); }


class IQTDat : public boost::noncopyable
    {
    public:

    typedef std::list<ITensor>
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

    explicit 
    IQTDat(std::istream& s);

    //
    // Accessors
    //

    const_iterator
    begin() const { return itensor.begin(); }

    iterator
    begin() { uninit_rmap(); return itensor.begin(); }

    const_iterator
    end() const { return itensor.end(); }

    iterator
    end() { uninit_rmap(); return itensor.end(); }

    const ITensor&
    get(const ApproxReal& r) const { return *rmap[r]; }

    ITensor&
    get(const ApproxReal& r) { return *rmap[r]; }

    int
    size() const { return itensor.size(); }

    bool
    empty() const { return itensor.empty(); }

    void
    clear();

    void 
    insert(const ApproxReal& r, const ITensor& t);

    void 
    insert(const ITensor& t);

    void 
    insert_add(const ApproxReal& r, const ITensor& t);

    void 
    insert_add(const ITensor& t);

    //void 
    //insert_assign(const ITensor& t);

    void 
    clean(Real min_norm);

    bool 
    has_itensor(const ApproxReal& r) const;

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

    //void* operator 
    //new(size_t size) 
    //    throw(std::bad_alloc)
    //    { return allocator().alloc(); }

    //void operator 
    //delete(void* p) 
    //    throw()
    //    { return allocator().dealloc(p); }

    static const boost::shared_ptr<IQTDat>& 
    Null();

    private:

    //////////////
    //
    // Data Members
    //

    mutable
    StorageT itensor;

    mutable std::map<ApproxReal,iterator>
    rmap; //mutable so that const IQTensor methods can use rmap

    mutable 
    bool rmap_init;

    //
    //////////////

    void 
    init_rmap() const;

    void 
    uninit_rmap() const;

    //Not copyable with =
    void operator=(const IQTDat&);

    //static DatAllocator<IQTDat>& allocator()
    //    {
    //    static DatAllocator<IQTDat> allocator_;
    //    return allocator_;
    //    };

    }; //class IQTDat

template <typename Callable> 
void IQTensor::
mapElems(const Callable& f)
    {
    solo();
    Foreach(ITensor& t, dat.nc()) 
        t.mapElems(f);
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
void 
BraKet(IQTensor x, const IQTensor& y, Real& re, Real& im);

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

std::ostream& 
operator<<(std::ostream & s, const IQTensor &t);


#endif
