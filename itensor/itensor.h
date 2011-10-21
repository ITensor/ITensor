#ifndef __ITENSOR_ITENSOR_H
#define __ITENSOR_ITENSOR_H
#include "types.h"
#include "real.h"
#include "allocator.h"
#include "index.h"
#include "permutation.h"
#include "prodstats.h"

//Forward declarations
struct ProductProps;
class Counter;
class Combiner;
class ITDat;

//
// ITensor
//
class ITensor
    {
public:

    //Accessor Methods ----------------------------------------------

    //unique_Real depends on indices only, unordered:
    inline Real 
    unique_Real() const { return ur; } 

    const Index& 
    index(int j) const;

    inline int 
    r() const { return r_; }

    int 
    m(int j) const;

    inline bool 
    is_null() const { return (p == 0); }

    inline bool 
    is_not_null() const { return (p != 0); }

    inline bool 
    is_complex() const { return hasindexn(Index::IndReIm()); }

    inline bool 
    is_not_complex() const { return !hasindexn(Index::IndReIm()); }

    inline LogNumber 
    scale() const { return scale_; }

    //Can be used for iteration over Indices in a foreach loop
    //e.g. foreach(const Index& I, t.index() ) { ... }
    typedef boost::array<Index,NMAX+1>::const_iterator index_it;
    inline const std::pair<index_it,index_it> 
    index() const  
        { return std::make_pair(index_.begin()+1,index_.begin()+r_+1); }


    //Constructors --------------------------------------------------

    ITensor();

    ITensor(Real val);

    explicit 
    ITensor(const Index& i1);

    ITensor(const Index& i1, Real val);

    ITensor(const Index& i1, const Vector& V);

    ITensor(Index i1,Index i2);

    //Create an ITensor as a matrix with 'a' on the diagonal
    ITensor(Index i1,Index i2,Real a);

    ITensor(Index i1,Index i2,const MatrixRef& M);

    ITensor(Index i1, Index i2, Index i3,
            Index i4 = Index::Null(), 
            Index i5 = Index::Null(), 
            Index i6 = Index::Null(),
            Index i7 = Index::Null(), 
            Index i8 = Index::Null());

    explicit 
    ITensor(const IndexVal& iv, Real fac = 1);

    ITensor(const IndexVal& iv1, const IndexVal& iv2);

    ITensor(const IndexVal& iv1, const IndexVal& iv2, 
            const IndexVal& iv3, const IndexVal& iv4 = IndexVal::Null(), 
            const IndexVal& iv5 = IndexVal::Null(), const IndexVal& iv6 = IndexVal::Null(), 
            const IndexVal& iv7 = IndexVal::Null(), const IndexVal& iv8 = IndexVal::Null());

    explicit 
    ITensor(const std::vector<Index>& I);

    ITensor(const std::vector<Index>& I, const Vector& V);

    ITensor(const std::vector<Index>& I, const ITensor& other);

    ITensor(const std::vector<Index>& I, const ITensor& other, Permutation P);

    ITensor(std::istream& s) { read(s); }

    static const ITensor& 
    Complex_1()
        {
        static const ITensor Complex_1_(makeComplex_1);
        return Complex_1_;
        }

    static const ITensor& 
    Complex_i()
        {
        static const ITensor Complex_i_(makeComplex_i);
        return Complex_i_;
        }

    static const ITensor& 
    ConjTensor()
        {
        static const ITensor ConjTensor_(makeConjTensor);
        return ConjTensor_;
        }

    void 
    read(std::istream& s);

    void
    write(std::ostream& s) const;


    //Operators -------------------------------------------------------

    ITensor& 
    operator*=(const ITensor& other);

    ITensor 
    operator*(ITensor other) const { other *= *this; return other; }

    ITensor& 
    operator*=(const IndexVal& iv) 
        { ITensor oth(iv); return operator*=(oth); } 

    ITensor 
    operator*(const IndexVal& iv) const 
        { ITensor res(*this); res *= iv; return res; }

    friend inline ITensor 
    operator*(const IndexVal& iv, ITensor t) 
        { return (t *= iv); }

    ITensor& 
    operator*=(Real fac) { scale_ *= fac; return *this; }

    ITensor 
    operator*(Real fac) const 
        { ITensor res(*this); res *= fac; return res; }

    friend inline ITensor 
    operator*(Real fac, ITensor t) 
        { return (t *= fac); }

    ITensor& 
    operator/=(Real fac) { scale_ /= fac; return *this; }

    ITensor 
    operator/(Real fac) const 
        { ITensor res(*this); res /= fac; return res; }

    friend inline ITensor 
    operator/(Real fac, ITensor t) 
        { return (t /= fac); }

    //operator/=(ITensor) is actually non-contracting product
    ITensor& 
    operator/=(const ITensor& other);

    ITensor 
    operator/(const ITensor& other) const 
        { ITensor res(*this); res /= other; return res; }

    ITensor& 
    operator+=(const ITensor& o);

    ITensor 
    operator+(const ITensor& o) const 
        { ITensor res(*this); res += o; return res; }

    ITensor& 
    operator-=(const ITensor& o)
        {
        if(this == &o) { scale_ = 0; return *this; }
        scale_ *= -1; operator+=(o); scale_ *= -1; return *this; 
        }

    ITensor 
    operator-(const ITensor& o) const 
        { ITensor res(*this); res -= o; return res; }


    //Index Methods ---------------------------------------------------

    Index 
    findtype(IndexType t) const;

    bool 
    findtype(IndexType t, Index& I) const;

    int 
    findindex(const Index& I) const;

    int 
    findindexn(const Index& I) const;

    int 
    findindex1(const Index& I) const;

    bool 
    has_common_index(const ITensor& other) const;
    
    bool 
    hasindex(const Index& I) const;

    bool 
    hasindexn(const Index& I) const;

    bool 
    hasindex1(const Index& I) const;

    inline bool 
    notin(const Index& I) const { return !hasindex(I); }

    void 
    addindex1(const std::vector<Index>& indices);

    void 
    addindex1(const Index& I);

    //Removes the jth index as found by findindex
    void 
    removeindex1(int j);

    inline void 
    removeindex1(const Index& I) 
        { removeindex1(findindex1(I)); }


    //Primelevel Methods ------------------------------------

    void 
    noprime(PrimeType p = primeBoth);

    void 
    doprime(PrimeType pt, int inc = 1);

    inline void 
    primeall() { doprime(primeBoth,1); }

    inline void 
    primesite() { doprime(primeSite,1); }

    inline void 
    primelink() { doprime(primeLink,1); }

    void 
    mapprime(int plevold, int plevnew, PrimeType pt = primeBoth);

    void 
    mapprimeind(const Index& I, int plevold, int plevnew, 
                PrimeType pt = primeBoth);

    inline void 
    primeind(const Index& I, int inc = 1)
        { mapindex(I,primed(I)); }

    void 
    primeind(const Index& I, const Index& J);

    inline void 
    noprimeind(const Index& I) { mapindex(I,I.deprimed()); }

    friend inline ITensor 
    primed(ITensor A, int inc = 1)
        { A.doprime(primeBoth,inc); return A; }

    friend inline ITensor 
    primesite(ITensor A, int inc = 1)
        { A.doprime(primeSite,inc); return A; }

    friend inline ITensor 
    primelink(ITensor A, int inc = 1)
        { A.doprime(primeLink,inc); return A; }

    friend inline ITensor 
    primeind(ITensor A, const Index& I)
        { A.mapindex(I,primed(I)); return A; }

    friend ITensor 
    primeind(ITensor A, const Index& I1, const Index& I2);

    friend inline ITensor 
    deprimed(ITensor A) { A.noprime(); return A; }


    //Element Access Methods ----------------------------------------

    Real 
    val0() const;

    Real 
    val1(int i1) const;

    Real& 
    operator()();

    Real 
    operator()() const;

    Real& 
    operator()(const IndexVal& iv1);

    Real operator()(const IndexVal& iv1) const;

    Real& 
    operator()(const IndexVal& iv1, const IndexVal& iv2);

    Real 
    operator()(const IndexVal& iv1, const IndexVal& iv2) const;

    Real& 
    operator()(const IndexVal& iv1, const IndexVal& iv2, 
               const IndexVal& iv3, const IndexVal& iv4 = IndexVal::Null(), 
               const IndexVal& iv5 = IndexVal::Null(),const IndexVal& iv6 = IndexVal::Null(),
               const IndexVal& iv7 = IndexVal::Null(),const IndexVal& iv8 = IndexVal::Null());

    Real 
    operator()(const IndexVal& iv1, const IndexVal& iv2, 
               const IndexVal& iv3, const IndexVal& iv4 = IndexVal::Null(), 
               const IndexVal& iv5 = IndexVal::Null(),const IndexVal& iv6 = IndexVal::Null(),
               const IndexVal& iv7 = IndexVal::Null(),const IndexVal& iv8 = IndexVal::Null()) 
    const;


    //Methods for Mapping to Other Objects ----------------------------------

    // Assume *this and other have same indices but different order.
    // Copy other into *this, without changing the order of indices in either
    // operator= would put the order of other into *this
    void 
    assignFrom(const ITensor& other);

    void 
    groupIndices(const boost::array<Index,NMAX+1>& indices, int nind, 
                      const Index& grouped, ITensor& res) const;

    void 
    expandIndex(const Index& small, const Index& big, 
                     int start, ITensor& res) const;

    void 
    fromMatrix11(const Index& i1, const Index& i2, const Matrix& res);

    void 
    toMatrix11NoScale(const Index& i1, const Index& i2, 
                           Matrix& res) const;
    void 
    toMatrix11(const Index& i1, const Index& i2, Matrix& res) const;

    /*
    // group i1,i2; i3,i4
    void toMatrix22(const Index& i1, const Index& i2, 
                    const Index& i3, const Index& i4,Matrix& res) const;
    void fromMatrix22(const Index& i1, const Index& i2, 
                      const Index& i3, const Index& i4,const Matrix& res);

    // group i1,i2; i3
    void toMatrix21(const Index& i1, const Index& i2, 
                    const Index& i3, Matrix& res) const;
    void fromMatrix21(const Index& i1, const Index& i2, 
                      const Index& i3, const Matrix& res);

    // group i1; i2,i3
    void toMatrix12(const Index& i1, const Index& i2, 
                    const Index& i3, Matrix& res) const;
    void fromMatrix12(const Index& i1, const Index& i2, 
                      const Index& i3, const Matrix& res);
    */

    int 
    vec_size() const;

    void 
    assignToVec(VectorRef v) const;

    void 
    assignFromVec(const VectorRef& v);

    void 
    reshapeDat(const Permutation& p, Vector& rdat) const;

    void 
    reshapeTo(const Permutation& P, ITensor& res) const;

    void 
    reshape(const Permutation& P);


    //Other Methods -------------------------------------------------

    void 
    Randomize();

    void 
    SplitReIm(ITensor& re, ITensor& im) const;

    inline void 
    conj() 
        { 
        if(!is_complex()) return; 
        operator/=(ITensor::ConjTensor()); 
        }

    inline bool 
    is_zero() const { return (norm() < 1E-20); } 

    Real 
    sumels() const;

    Real 
    norm() const;

    void 
    scaleOutNorm() const;

    void 
    scaleTo(LogNumber newscale) const;

    void 
    print(std::string name = "",Printdat pdat = HideData) const;

    friend std::ostream& 
    operator<<(std::ostream & s, const ITensor & t);

    friend class commaInit;

    typedef Index IndexT;
    typedef IndexVal IndexValT;
    typedef Combiner CombinerT;
    static const Index& ReImIndex() { return Index::IndReIm(); }

private:

    //mutable: const methods may want to reshape data
    mutable boost::intrusive_ptr<ITDat> p; 

    //Indices, maximum of 8 (index_[0] not used), mutable to allow reordering
    mutable boost::array<Index,NMAX+1> index_; 

    int r_,rn_;
    //mutable since e.g. scaleTo is logically const
    mutable LogNumber scale_; 

    Real ur;

    void 
    initCounter(Counter& C) const;

    void 
    allocate(int dim);

    void 
    allocate();

#ifdef DO_ALT
    void 
    newAltDat(const Permutation& P) const;

    PDat& 
    lastAlt() const;
#endif

    //Disattach self from current ITDat and create own copy instead.
    //Necessary because ITensors logically represent distinct
    //objects even though they may share data in reality.
    void 
    solo() const;
    
    void 
    set_unique_Real();

    void 
    _construct1(const Index& i1);

    void 
    _construct2(const Index& i1, const Index& i2);

    template<class Iterable>
    int 
    fillFromIndices(const Iterable& I, int size);

    //Prefer to map via a Combiner
    //Though 'mapindex' is useful for tested, internal calls
    void 
    mapindex(const Index& i1, const Index& i2);

    void 
    getperm(const boost::array<Index,NMAX+1>& oth_index_, Permutation& P) const;

    friend struct ProductProps;

    friend void toMatrixProd(const ITensor& L, const ITensor& R, 
                             const ProductProps& pp,
                             MatrixRefNoLink& lref, MatrixRefNoLink& rref);


    int _ind(int i1, int i2, int i3, int i4, 
             int i5, int i6, int i7, int i8) const;

    int _ind2(const IndexVal& iv1, const IndexVal& iv2) const;

    int _ind8(const IndexVal& iv1, const IndexVal& iv2, 
              const IndexVal& iv3, const IndexVal& iv4 = IndexVal::Null(), 
              const IndexVal& iv5 = IndexVal::Null(),const IndexVal& iv6 = IndexVal::Null(),
              const IndexVal& iv7 = IndexVal::Null(),const IndexVal& iv8 = IndexVal::Null())
        const;

public:

    enum ITmaker { makeComplex_1, makeComplex_i, makeConjTensor };

    ITensor(ITmaker itm);

    }; // class ITensor

//
// Counter
//
class Counter
    {
public:
    boost::array<int,NMAX+1> n;
    boost::array<int,NMAX+1> i;
    int ind;

    Counter();
    Counter(const boost::array<Index,NMAX+1>& ii,int rn,int r);

    void init(const boost::array<Index,NMAX+1>& ii, int rn, int r);

    Counter& operator++();

    bool operator!=(const Counter& other) const;
    bool operator==(const Counter& other) const;

    inline bool notDone() const 
	{ return i[1] != 0; }

    friend inline std::ostream& operator<<(std::ostream& s, const Counter& c);

private:
    void reset(int a)
	{
        i.assign(a);
        ind = 1;
	}
    int rn_,r_;
    };

//#define DO_ALT
#ifdef DO_ALT
struct PDat
    {
    Permutation I; 
    Vector v;
    PDat(const Permutation& P_, const Vector& v_) 
		: I(P_.inverse()), v(v_) { }
    PDat(const Permutation& P_) : I(P_.inverse()) { }
    };
#endif

//
// ITDat
//
class ITDat
    {
private:
    mutable unsigned int numref;
    static DatAllocator<ITDat> allocator;
    void operator=(const ITDat&);
    ~ITDat() { } //must be dynamically allocated
public:
    Vector v;
#ifdef DO_ALT
    std::vector<PDat> alt;
#endif

    ITDat() : numref(0), v(0) { }

    explicit 
    ITDat(int size) : numref(0), v(size)
	{ assert(size > 0); v = 0; }

    explicit 
    ITDat(const Vector& v_) : numref(0), v(v_) { }

    explicit 
    ITDat(Real r) : numref(0), v(1)
	{ v = r; }

    explicit 
    ITDat(std::istream& s) : numref(0) 
	{ read(s); }

    explicit 
    ITDat(const ITDat& other) : numref(0), v(other.v) { }

    void read(std::istream& s)
	{ 
	int size = 0;
	s.read((char*) &size,sizeof(size));
	v.ReDimension(size);
	s.read((char*) v.Store(), sizeof(Real)*size);
	}

    void write(std::ostream& s) const 
	{ 
	const int size = v.Length();
	s.write((char*) &size, sizeof(size));
	s.write((char*) v.Store(), sizeof(Real)*size); 
	}
    
    void print() const 
	{ std::cout << "ITDat: v = " << v; }

    inline void* operator 
    new(size_t) throw(std::bad_alloc)
        { return allocator.alloc(); }

    inline void operator 
    delete(void* p) throw()
        { return allocator.dealloc(p); }

    friend class ITensor;
    ENABLE_INTRUSIVE_PTR(ITDat)
    };


Real Dot(const ITensor& x, const ITensor& y, bool doconj = true);

void Dot(const ITensor& x, const ITensor& y, Real& re, Real& im, 
                bool doconj = true);

inline ITensor 
operator*(const IndexVal& iv1, const IndexVal& iv2) 
    { ITensor t(iv1); return (t *= iv2); }

inline ITensor 
operator*(const IndexVal& iv1, Real fac) 
    { return ITensor(iv1,fac); }

inline ITensor 
operator*(Real fac, const IndexVal& iv) 
    { return ITensor(iv,fac); }

// Given Tensors which represent operators 
//(e.g. A(site-1',site-1), B(site-1',site-1), 
// Multiply them, fixing primes C(site-1',site-1)
// a * b  (a above b in diagram, unprimed = right index of matrix)
template<class Tensor>
inline Tensor 
multSiteOps(Tensor a, Tensor b) 
    {
    a.mapprime(1,2,primeSite);
    a.mapprime(0,1,primeSite);
    Tensor res = a * b;
    res.mapprime(2,1,primeSite);
    return res;
    }

class commaInit
    {
public:
    commaInit(ITensor& T_)
        : T(T_)
        { 
        if(T.is_null()) Error("Can't assign to null ITensor");
        T.solo();
        T.scaleTo(1);
        T.initCounter(c);
        }

    commaInit& operator<<(Real r)
        {
        return operator,(r);
        }

    commaInit& operator,(Real r)
        {
        if(c.notDone()) 
            { T.p->v(c.ind) = r; ++c; }
        else 
            { Error("Comma assignment list too long.\n"); }
        return *this;
        }

private:
    ITensor& T;
    Counter c; 
    };


template<class Iterable>
int ITensor::
fillFromIndices(const Iterable& I, int size)
    {
    r_ = size;
    assert(r_ <= NMAX);
    assert(rn_ == 0);
    int r1_ = 0;
    boost::array<const Index*,NMAX+1> index1_;
    int alloc_size = 1;
    for(int n = 0; n < r_; ++n)
        {
        const Index& i = I[n];
        DO_IF_DEBUG(if(i == Index::Null()) Error("ITensor: null Index in constructor.");)
        if(i.m()==1) 
            { GET(index1_,++r1_) = &i; }
        else         
            { 
            GET(index_, ++rn_) = i; 
            alloc_size *= i.m(); 
            }
        }
    for(int l = 1; l <= r1_; ++l) 
        index_[rn_+l] = *(index1_[l]);
    set_unique_Real();
    return alloc_size;
    }


#endif
