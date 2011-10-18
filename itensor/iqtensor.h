#ifndef __IQ_H
#define __IQ_H
#include "itensor.h"
#include "iqindex.h"
#include <list>
#include <map>

class IQTDat;

class IQCombiner;

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
    is_null() const { return p == 0; }

    inline bool 
    is_not_null() const { return p != 0; }

    int 
    num_index() const;
    
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

    IQTensor(ITmaker itm);

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
        static const IQTensor Complex_1_(makeComplex_1);
        return Complex_1_;
        }

    static const IQTensor& Complex_i()
        {
        static const IQTensor Complex_i_(makeComplex_i);
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

    IQTensor 
    operator-(const IQTensor& o) const 
        { IQTensor res(o); res *= -1; res += *this; return res; }

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

    //
    // Multiplication by an ITensor
    //
    ITensor 
    operator*(const ITensor& t) const
        { ITensor res(*this); res *= t; return res; }
    
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
        { ITensor res(*this); res *= iv; return res; }

    friend inline ITensor 
    operator*(const IndexVal& iv, const IQTensor& T) 
        { return ITensor(iv) * T; }

    //Convert to ITensor
    operator ITensor() const;

    IQTensor& 
    operator+=(const ITensor& t);

    //Like operator+=(ITensor) but
    //demands that the block is zero
    //before inserting
    void insert(const ITensor& t);

    Real& 
    operator()(const IQIndexVal& iv1, 
               const IQIndexVal& iv2 = IQIndexVal::Null(), 
	           const IQIndexVal& iv3 = IQIndexVal::Null(), 
               const IQIndexVal& iv4 = IQIndexVal::Null(), 
	           const IQIndexVal& iv5 = IQIndexVal::Null(), 
               const IQIndexVal& iv6 = IQIndexVal::Null(),
	           const IQIndexVal& iv7 = IQIndexVal::Null(), 
               const IQIndexVal& iv8 = IQIndexVal::Null());

    //----------------------------------------------------
    //IQTensor quantum number methods

    QN 
    div() const;

    bool 
    checkDiv(QN expected) const;

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
    primeind(const IQIndex& I);

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
    findindex(const IQIndex& i) const;

    bool 
    hastype(IndexType t) const;

    const IQIndex& 
    findtype(IndexType t) const;

    const IQIndex& 
    finddir(Arrow dir) const;

    bool 
    hasindex(const IQIndex& i) const;

    bool 
    is_complex() const 
        { return findindex(IQIndex::IndReIm()) != 0; }

    void 
    addindex1(const IQIndex& I);

    //----------------------------------------------------
    //IQTensor miscellaneous methods

    Real 
    unique_Real() const;

    int 
    num_index(IndexType t) const;

    Real 
    norm() const;

    Real 
    sumels() const;

    void 
    scaleOutNorm() const;

    void 
    scaleTo(LogNumber newscale) const;

    inline void 
    clean(Real min_norm = MIN_CUT);

    int 
    vec_size() const;

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
    printIQInds(const std::string& name = "") const;
    inline void 
    printIQInds(const boost::format& fname) const
        { printIQInds(fname.str()); }

    void 
    assignFrom(const IQTensor& other) const;

    void 
    SplitReIm(IQTensor& re, IQTensor& im) const;

    void 
    conj();

    friend std::ostream& 
    operator<<(std::ostream & s, const IQTensor &t);

    typedef IQIndex IndexT;

    typedef IQIndexVal IndexValT;

    typedef IQCombiner CombinerT;


    static const IQIndex& ReImIndex()
        { return IQIndex::IndReIm(); }

private:

    boost::intrusive_ptr<IQTDat> p;

    void solo();

}; //class IQTensor

class IQTDat
{
public:

    IQTDat();

    explicit 
    IQTDat(const IQIndex& i1);

    IQTDat(const IQIndex& i1, const IQIndex& i2);

    IQTDat(const IQIndex& i1, const IQIndex& i2, const IQIndex& i3);

    IQTDat(const IQIndex& i1, const IQIndex& i2, 
	       const IQIndex& i3, const IQIndex& i4,
	       const IQIndex& i5 = IQIndex::Null(), 
	       const IQIndex& i6 = IQIndex::Null(), 
	       const IQIndex& i7 = IQIndex::Null(), 
	       const IQIndex& i8 = IQIndex::Null());

    explicit 
    IQTDat(std::vector<IQIndex>& iqinds_);

    explicit 
    IQTDat(const IQTDat& other);

    explicit 
    IQTDat(std::istream& s);

    void 
    read(std::istream& s);

    void 
    write(std::ostream& s) const;

    void 
    init_rmap() const;

    void 
    uninit_rmap() const;

    bool 
    has_itensor(const ApproxReal& r) const;

    void 
    insert_itensor(const ApproxReal& r, const ITensor& t);

    void 
    clean(Real min_norm);

    inline void* operator 
    new(size_t size) 
        throw(std::bad_alloc)
        { return allocator.alloc(); }

    inline void operator 
    delete(void* p) 
        throw()
        { return allocator.dealloc(p); }

public:

    mutable std::list<ITensor> 
    itensor; // This is mutable to allow reordering

    std::vector<IQIndex> 
    iqindex_;

    mutable std::map<ApproxReal,std::list<ITensor>::iterator> 
    rmap; //mutable so that const IQTensor methods can use rmap

    ENABLE_INTRUSIVE_PTR(IQTDat)

private:

    ~IQTDat() { } //must be dynamically allocated

    void operator=(const IQTDat&);

    static DatAllocator<IQTDat> 
    allocator;

    mutable unsigned int 
    numref;

    mutable bool 
    rmap_init;

}; //class IQTDat


Real 
ReSingVal(const IQTensor& x);

Real 
Dot(const IQTensor& x, const IQTensor& y, bool doconj = true);

void 
Dot(const IQTensor& x, const IQTensor& y, Real& re, Real& im, bool doconj = true);

inline bool 
checkQNs(const ITensor& t) 
    { return true; }

//Checks if all IQTensor blocks have the same divergence
bool 
checkQNs(const IQTensor& T);

template<class T> 
class Printit
    {
public:
    std::ostream& s;
    std::string spacer;
    Printit(std::ostream& _s, std::string _spacer) : s(_s), spacer(_spacer) {}
    void operator()(const T& t) { s << t << spacer; }
    };

#endif
