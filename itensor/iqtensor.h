//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __IQ_H
#define __IQ_H
#include "iqindexset.h"
#include <list>
#include <map>

class IQTDat;
class IQCombiner;
class IQTSparse;


//
// IQTensor
//

class IQTensor
    {
    public:

    typedef std::list<ITensor>::iterator 
    iten_it;

    typedef std::list<ITensor>::const_iterator 
    const_iten_it;

    typedef std::vector<IQIndex>::iterator 
    iqind_it;

    typedef std::vector<IQIndex>::const_iterator 
    const_iqind_it;

    int 
    r() const;

    const IQIndex& 
    index(int j) const;

    int 
    iten_size() const;

    bool 
    iten_empty() const;

    inline bool 
    isNull() const { return p == 0; }

    inline bool 
    isNotNull() const { return p != 0; }

    int 
    num_index() const;

    const IQTDat&
    blocks() const { return *p; }
    
    //----------------------------------------------------
    //IQTensor: iterators 
    const_iten_it 
    const_iten_begin() const;

    const_iten_it 
    const_iten_end() const;

    std::pair<const_iten_it,const_iten_it> 
    itensors() const;

    const_iqind_it 
    const_iqind_begin() const;

    const_iqind_it 
    const_iqind_end() const;

    std::pair<const_iqind_it,const_iqind_it> 
    iqinds() const;

    //----------------------------------------------------
    //IQTensor: Constructors

    IQTensor();

    explicit 
    IQTensor(Real val);

    explicit 
    IQTensor(const IQIndex& i1);

    IQTensor(const IQIndex& i1,const IQIndex& i2);

    IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3);

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

    explicit 
    IQTensor(std::vector<IQIndex>& iqinds_);

    explicit
    IQTensor(const IQIndexVal& iv1);

    IQTensor(const IQIndexVal& iv1, const IQIndexVal& iv2);

    IQTensor(const IQIndexVal& iv1, const IQIndexVal& iv2,
             const IQIndexVal& iv3);

    IQTensor(ITensor::ITmaker itm);

    IQTensor(IQmaker i);

    IQTensor(PrimeType pt, const IQTensor& other);

    explicit 
    IQTensor(std::istream& s);

    static const IQTensor& Sing()
        {
        static const IQTensor Sing_(makeSing);
        return Sing_;
        }

    static const IQTensor& Complex_1()
        {
        static const IQTensor Complex_1_(ITensor::makeComplex_1);
        return Complex_1_;
        }

    static const IQTensor& Complex_i()
        {
        static const IQTensor Complex_i_(ITensor::makeComplex_i);
        return Complex_i_;
        }

    void read(std::istream& s);

    void write(std::ostream& s) const;


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

    IQTensor& 
    operator*=(const LogNumber& lgnum);

    IQTensor 
    operator*(const LogNumber& lgnum) const
        { IQTensor res(*this); res *= lgnum; return res; }

    friend inline IQTensor 
    operator*(const LogNumber& lgnum, IQTensor T) 
        { T *= lgnum; return T; }

    //
    // Multiplication by an ITensor
    //
    ITensor 
    operator*(const ITensor& t) const
        { ITensor res = this->toITensor(); res *= t; return res; }
    
    //
    // Multiplication by an IQIndexVal
    //
    IQTensor 
    operator*(const IQIndexVal& iv) const
        { IQTensor res(*this); res *= IQTensor(iv); return res; }

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
    //operator ITensor() { return toITensor(); }

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
    //IQTensor quantum number methods

    QN 
    div() const;

    friend void 
    checkDiv(const IQTensor& T, QN expected = QN());

    QN 
    qn(const Index& in) const;

    Arrow 
    dir(const Index& in) const;


    //----------------------------------------------------
    //IQTensor: prime methods

    void ind_inc_prime(const IQIndex& i,int inc);

    void noprime(PrimeType pt = primeBoth);

    friend inline IQTensor 
    deprimed(IQTensor A) { A.noprime(); return A; }

    void 
    noprimelink();

    void 
    doprime(PrimeType pt, int inc = 1);

    //no need to keep prime level small
    void 
    mapprime(int plevold, int plevnew, PrimeType pt = primeBoth);

    void 
    primeind(const IQIndex& I, int inc = 1);

    friend inline IQTensor 
    primeind(IQTensor A, const IQIndex& I)
        { A.primeind(I); return A; }

    friend inline IQTensor 
    primeind(IQTensor A, const IQIndex& I, const IQIndex& J)
        { A.primeind(I); A.primeind(J); return A; }

    void 
    noprimeind(const IQIndex& I);

    friend inline IQTensor 
    primed(IQTensor A) { A.doprime(primeBoth); return A; }

    void 
    primesite() { doprime(primeSite); }

    friend inline IQTensor 
    primesite(IQTensor A) { A.doprime(primeSite); return A; }

    void 
    primelink() { doprime(primeLink); }

    friend inline IQTensor 
    primelink(const IQTensor& A)
        { IQTensor res(A); res.doprime(primeLink); return res; }


    //----------------------------------------------------
    //IQTensor index methods

    int 
    find_iqind(const Index& I) const;

    //Return true if one of the ITensors uses this Index
    bool 
    uses_ind(const Index& i) const;

    int 
    findindex(const IQIndex& I) const;

    bool 
    hastype(IndexType t) const;

    const IQIndex& 
    findtype(IndexType t) const;

    const IQIndex& 
    finddir(Arrow dir) const;

    bool 
    hasindex(const IQIndex& I) const;

    bool 
    isComplex() const 
        { return findindex(IQIndex::IndReIm()) != 0; }

    void 
    addindex1(const IQIndex& I);

    void
    tieIndices(const boost::array<IQIndex,NMAX+1>& indices, int nind, const IQIndex& tied);

    void
    tieIndices(const IQIndex& i1, const IQIndex& i2, const IQIndex& tied);

    friend inline IQTensor
    tieIndices(const IQIndex& i1, const IQIndex& i2, 
               const IQIndex& tied, IQTensor T)
        { T.tieIndices(i1,i2,tied); return T; }

    void
    trace(const boost::array<IQIndex,NMAX+1>& indices, int nind);

    void
    trace(const IQIndex& i1, const IQIndex& i2);

    void
    trace(const IQIndex& i1);

    IQTensor friend inline
    trace(const IQIndex& i1, const IQIndex& i2, IQTensor T)
        { T.trace(i1,i2); return T; }

    //
    // Tracing over all indices results in a Real
    //
    Real friend
    trace(IQTensor T);

    //----------------------------------------------------
    //IQTensor miscellaneous methods

    void
    symmetricDiag11(const IQIndex& i1, IQTensor& D, IQTensor& U, IQIndex& mid) const;

    void
    symmetricDiag11(const IQIndex& i1, IQTensor& D, IQTensor& U, IQIndex& mid, int& mink, int& maxk) const;

    Real 
    uniqueReal() const;

    Real 
    norm() const;

    Real 
    sumels() const;

    template <typename Callable> void
    mapElems(const Callable& f);

    void 
    scaleOutNorm() const;

    void 
    scaleTo(LogNumber newscale) const;

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

    void 
    GetSingComplex(Real& re, Real& im) const;
    
    void 
    Randomize();

    void 
    print(std::string name = "",Printdat pdat = HideData) const;

    void 
    printIndices(const std::string& name = "") const;
    inline void 
    printIndices(const boost::format& fname) const
        { printIndices(fname.str()); }

    void 
    assignFrom(const IQTensor& other);

    void 
    SplitReIm(IQTensor& re, IQTensor& im) const;

    void 
    conj();

    void 
    conj(const IQIndex& I);

    void
    swap(IQTensor& other);

    friend std::ostream& 
    operator<<(std::ostream & s, const IQTensor &t);

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

    static const IQIndex& ReImIndex()
        { return IQIndex::IndReIm(); }

    private:

    /////////////////
    // 
    // Data Members

    boost::intrusive_ptr<IQIndexSet>
    is_;

    boost::intrusive_ptr<IQTDat> 
    p;

    //
    /////////////////

    void 
    soloIndex();

    void 
    soloDat();

    void 
    solo();

    //Workaround to ensure
    //that p is treated as const
    //by the compiler
    const IQTDat&
    dat() const { return *p; }

    IQTDat&
    ncdat() const { return *p; }

    }; //class IQTensor

class IQTDat
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

    void 
    insert_assign(const ITensor& t);

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
    write(std::ostream& s) const;

    //void* operator 
    //new(size_t size) 
    //    throw(std::bad_alloc)
    //    { return allocator().alloc(); }

    //void operator 
    //delete(void* p) 
    //    throw()
    //    { return allocator().dealloc(p); }


    friend void 
    intrusive_ptr_add_ref(IQTDat* p);

    friend void 
    intrusive_ptr_release(IQTDat* p);

    int 
    count() const { return numref; }

    static IQTDat* Null()
        {
        //Set initial numref to 1000
        static IQTDat Null_(1000);
        return &Null_;
        }

    private:

    //////////////
    //
    // Data Members
    //

    mutable std::list<ITensor> 
    itensor;

    mutable std::map<ApproxReal,iterator>
    rmap; //mutable so that const IQTensor methods can use rmap

    mutable unsigned int 
    numref;

    mutable bool 
    rmap_init;

    //
    //////////////

    void 
    init_rmap() const;

    void 
    uninit_rmap() const;

    explicit
    IQTDat(int init_numref);

    //Must be dynamically allocated
    ~IQTDat() { }
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
    Foreach(ITensor& t, *p)
        t.mapElems(f);
    }

Real 
ReSingVal(const IQTensor& x);


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
Dot(const IQTensor& x, const IQTensor& y);

//
// Computes the scalar (inner) product of two
// possibly complex ITensors.
//
// The first argument gets conjugated so this method
// is equivalent to the IQTensor contraction conj(x) * y 
// except it yields two real numbers (re and im) instead 
// of a rank 0 IQTensor.
//
void 
BraKet(const IQTensor& x, const IQTensor& y, Real& re, Real& im);

//Checks if all IQTensor blocks have the same divergence
void 
checkQNs(const IQTensor& T);

//ITensor version for compatibility
void inline
checkQNs(const ITensor& t) { }

//ITensor version for compatibility
void inline
checkDiv(const ITensor& t, QN q = QN()) { } 

#endif
