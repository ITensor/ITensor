//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITENSOR_H
#define __ITENSOR_ITENSOR_H
#include "types.h"
#include "real.h"
#include "allocator.h"
#include "index.h"
#include "prodstats.h"
#include "indexset.h"

//Forward declarations
struct ProductProps;
class Counter;
class Combiner;
class ITDat;
class ITSparse;

//
// ITensor
//
class ITensor
    {
    public:

    //Accessor Methods ----------------------------------------------

    //uniqueReal depends on indices only, unordered:
    Real 
    uniqueReal() const { return is_.ur_; } 

    const Index& 
    index(int j) const { return is_.index(j); }

    int 
    r() const { return is_.r_; }

    int 
    rn() const { return is_.rn_; }

    int 
    m(int j) const { return is_.m(j); }

    bool 
    isNull() const { return (p == 0); }

    bool 
    isNotNull() const { return (p != 0); }

    bool 
    isComplex() const { return hasindexn(Index::IndReIm()); }

    bool 
    isNotComplex() const { return !hasindexn(Index::IndReIm()); }

    const LogNumber&
    scale() const { return scale_; }

    //Can be used for iteration over Indices in a Foreach loop
    //e.g. Foreach(const Index& I, t.index() ) { ... }
    const std::pair<IndexSet::index_it,IndexSet::index_it> 
    index() const  
        { return is_.index(); }


    //Constructors --------------------------------------------------

    ITensor();

    ITensor(Real val);

    explicit 
    ITensor(const Index& i1);

    ITensor(const Index& i1, Real val);

    ITensor(const Index& i1, const VectorRef& V);

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


    //
    // Operators
    //

    // Contracting product

    ITensor& 
    operator*=(const ITensor& other);

    ITensor 
    operator*(const ITensor& other) const 
        { ITensor res(*this); res *= other; return res; }

    // Contracting product with IndexVals
    // (sets an Index to a particular value)

    ITensor& 
    operator*=(const IndexVal& iv) 
        { return operator*=(ITensor(iv)); } 

    ITensor 
    operator*(const IndexVal& iv) const 
        { ITensor res(*this); res *= iv; return res; }

    friend inline ITensor 
    operator*(const IndexVal& iv, const ITensor& t) 
        { return (ITensor(iv) *= t); }

    // Multiplication and division by scalars

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

    // Non-contracting product

    ITensor& 
    operator/=(const ITensor& other);

    ITensor 
    operator/(const ITensor& other) const 
        { ITensor res(*this); res /= other; return res; }

    // Addition and subtraction

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
    findtype(IndexType t) const { return is_.findtype(t); }

    bool 
    findtype(IndexType t, Index& I) const { return is_.findtype(t,I); }

    int 
    findindex(const Index& I) const { return is_.findindex(I); }

    int 
    findindexn(const Index& I) const { return is_.findindexn(I); }

    int 
    findindex1(const Index& I) const { return is_.findindex1(I); }

    bool 
    has_common_index(const ITensor& other) const
        { return is_.has_common_index(other.is_); }
    
    bool 
    hasindex(const Index& I) const { return is_.hasindex(I); }

    bool 
    hasindexn(const Index& I) const { return is_.hasindexn(I); }

    bool 
    hasindex1(const Index& I) const { return is_.hasindex1(I); }

    bool
    hasAllIndex(const boost::array<Index,NMAX+1>& I, int nind) const
        { return is_.hasAllIndex(I,nind); }

    bool 
    notin(const Index& I) const { return !hasindex(I); }

    void 
    addindex1(const std::vector<Index>& indices) { is_.addindex1(indices); }

    void 
    addindex1(const Index& I) { is_.addindex1(I); }

    //Removes the jth index as found by findindex
    void 
    removeindex1(int j) { is_.removeindex1(j); }

    void 
    removeindex1(const Index& I) { is_.removeindex1(is_.findindex1(I)); }

    void 
    mapindex(const Index& i1, const Index& i2) { is_.mapindex(i1,i2); }

    //Primelevel Methods ------------------------------------

    void 
    noprime(PrimeType p = primeBoth) { is_.noprime(p); }

    void 
    doprime(PrimeType pt, int inc = 1) { is_.doprime(pt,inc); }

    void 
    primeall() { doprime(primeBoth,1); }

    void 
    primesite() { doprime(primeSite,1); }

    void 
    primelink() { doprime(primeLink,1); }

    void 
    mapprime(int plevold, int plevnew, PrimeType pt = primeBoth)
        { is_.mapprime(plevold,plevnew,pt); }

    void 
    mapprimeind(const Index& I, int plevold, int plevnew, 
                PrimeType pt = primeBoth)
        { is_.mapprimeind(I,plevold,plevnew,pt); }

    void 
    primeind(const Index& I, int inc = 1)
        { mapindex(I,primed(I,inc)); }

    void 
    primeind(const Index& I, const Index& J) { is_.primeind(I,J); }

    void 
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
    primeind(ITensor A, const Index& I, int inc = 1)
        { A.mapindex(I,primed(I,inc)); return A; }

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

    //
    // Assume *this and other have same indices but different order.
    // Copy other into *this, without changing the order of indices in either
    // operator= would put the order of other into *this
    //
    void 
    assignFrom(const ITensor& other);

    //
    // groupIndices combines a set of indices (of possibly different sizes) 
    // together, leaving only single grouped Index.
    //
    // RiJ = Ai(jk) <-- Here J represents the grouped pair of indices (jk)
    //                  If j.m() == 5 and k.m() == 7, J.m() == 5*7.
    //
    void 
    groupIndices(const boost::array<Index,NMAX+1>& indices, int nind, 
                      const Index& grouped, ITensor& res) const;

    //
    // tieIndices locks a set of indices (of the same size) together,
    // leaving only a single tied Index.
    //
    // Rijl = Aijil <-- here we have tied the 1st and 3rd index of A
    //
    void
    tieIndices(const boost::array<Index,NMAX+1>& indices, int nind,
               const Index& tied);

    void
    tieIndices(const Index& i1, const Index& i2,
               const Index& tied);

    friend inline ITensor
    tieIndices(const Index& i1, const Index& i2, 
               const Index& tied, ITensor T)
        { T.tieIndices(i1,i2,tied); return T; }

    void
    tieIndices(const Index& i1, const Index& i2,
               const Index& i3,
               const Index& tied);

    friend inline ITensor
    tieIndices(const Index& i1, const Index& i2, 
               const Index& i3, const Index& tied, ITensor T)
        { T.tieIndices(i1,i2,i3,tied); return T; }

    void
    tieIndices(const Index& i1, const Index& i2,
               const Index& i3, const Index& i4,
               const Index& tied);

    friend inline ITensor
    tieIndices(const Index& i1, const Index& i2, 
               const Index& i3, const Index& i4, 
               const Index& tied, ITensor T)
        { T.tieIndices(i1,i2,i3,i4,tied); return T; }

    //
    // The trace method sums over the given set of indices
    // (which must all have the same dimension).
    //
    // Rik = trace(j,l,m,Aijkml) = \sum_j Aijkjj
    //

    void
    trace(const boost::array<Index,NMAX+1>& indices, int nind);

    void
    trace(const Index& i1, const Index& i2);

    void
    trace(const Index& i1);

    ITensor friend inline
    trace(const Index& i1, const Index& i2, ITensor T)
        { T.trace(i1,i2); return T; }

    void
    trace(const Index& i1, const Index& i2, const Index& i3);

    ITensor friend inline
    trace(const Index& i1, const Index& i2, const Index& i3,
          ITensor T)
        { T.trace(i1,i2,i3); return T; }

    void
    trace(const Index& i1, const Index& i2, const Index& i3, const Index& i4);

    ITensor friend inline
    trace(const Index& i1, const Index& i2, const Index& i3, const Index& i4,
          ITensor T)
        { T.trace(i1,i2,i3,i4); return T; }

    //
    // Tracing over all indices results in a Real
    //
    Real friend inline
    trace(ITensor T)
        {
        if(T.rn() != 0) T.trace(T.is_.index_,T.rn());
        return T.val0();
        }

    //
    // expandIndex replaces a smaller index with a bigger one, padding out
    // the elements of the resulting ITensor with zeros as necessary.
    // Say we have a tensor Aij and j has range m. Now expand j with 
    // a larger Index J. The result is RiJ, where
    //        _
    // RiJ = |  Ai(j=J-start+1) for J = start...start+m
    //       |_ 0               otherwise
    //        

    void 
    expandIndex(const Index& small, const Index& big, int start);

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

    void
    symmetricDiag11(const Index& i1, ITensor& D, ITensor& U, Index& mid) const;

    void
    symmetricDiag11(const Index& i1, ITensor& D, ITensor& U, Index& mid, int& mink, int& maxk) const;

    int 
    vecSize() const;

    int 
    maxSize() const;

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

    //
    // Swap can be used for similar purposes
    // as operator=(const ITensor& other)
    // but is more efficient and has same
    // end result if other is just a temporary
    //
    void
    swap(ITensor& other);


    //Other Methods -------------------------------------------------

    void 
    Randomize();

    void 
    SplitReIm(ITensor& re, ITensor& im) const;

    void 
    conj() 
        { 
        if(!isComplex()) return; 
        operator/=(ITensor::ConjTensor()); 
        }

    void 
    conj(const Index& I) { }

    bool 
    is_zero() const { return (norm() < 1E-20); } 

    Real 
    sumels() const;

    Real 
    norm() const;

    template <typename Callable> void
    mapElems(const Callable& f);

    void
    pseudoInvert(Real cutoff = 0);

    void 
    scaleOutNorm() const;

    void 
    scaleTo(LogNumber newscale) const;

    void 
    print(std::string name = "",Printdat pdat = HideData) const;

    void 
    printIndices(const std::string& name = "") const
        { print(name,HideData); }

    void 
    printIndices(const boost::format& fname) const
        { printIndices(fname.str()); }

    friend std::ostream& 
    operator<<(std::ostream & s, const ITensor & t);

    friend class commaInit;

    typedef Index 
    IndexT;

    typedef IndexVal 
    IndexValT;

    typedef Combiner 
    CombinerT;

    typedef ITSparse
    SparseT;

    static const Index& 
    ReImIndex() { return Index::IndReIm(); }

    protected:

    //////////////
    //
    // Data Members
    //

    //mutable: const methods may want to reshape data
    mutable boost::intrusive_ptr<ITDat> p; 

    //Indices, maximum of 8 (is_.index_[0] not used)
    IndexSet is_;

    //mutable since e.g. scaleTo is logically const
    mutable LogNumber scale_; 

    //
    //
    //////////////

    void 
    initCounter(Counter& C) const;

    void 
    allocate(int dim);

    void 
    allocate();

    //Disattach self from current ITDat and create own copy instead.
    //Necessary because ITensors logically represent distinct
    //objects even though they may share data in reality.
    void 
    solo() const;
    
    friend struct ProductProps;

    friend void toMatrixProd(const ITensor& L, const ITensor& R, 
                             ProductProps& pp,
                             MatrixRefNoLink& lref, MatrixRefNoLink& rref);

    void
    directMultiply(const ITensor& other, ProductProps& pp, 
                   int& new_rn_, boost::array<Index,NMAX+1>& new_index_);

    int _ind(int i1, int i2, int i3, int i4, 
             int i5, int i6, int i7, int i8) const;

    int _ind2(const IndexVal& iv1, const IndexVal& iv2) const;

    int _ind8(const IndexVal& iv1, const IndexVal& iv2, 
              const IndexVal& iv3, const IndexVal& iv4 = IndexVal::Null(), 
              const IndexVal& iv5 = IndexVal::Null(),const IndexVal& iv6 = IndexVal::Null(),
              const IndexVal& iv7 = IndexVal::Null(),const IndexVal& iv8 = IndexVal::Null())
        const;

    friend class ITSparse;

    friend void 
    product(const ITSparse& S, const ITensor& T, ITensor& res);

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
    boost::array<int,NMAX+1> n, i;
    int ind;
    int rn_,r_;

    Counter();

    Counter(const boost::array<Index,NMAX+1>& ii,int rn,int r);

    Counter(const IndexSet& is);

    void 
    init(const boost::array<Index,NMAX+1>& ii, int rn, int r);

    void 
    init(const IndexSet& is);

    Counter& 
    operator++();

    bool 
    operator!=(const Counter& other) const;

    bool 
    operator==(const Counter& other) const;

    bool 
    notDone() const 
        { return i[1] != 0; }

    friend std::ostream& 
    operator<<(std::ostream& s, const Counter& c);

    void 
    reset(int a);

    };

//
// ITDat
//
class ITDat
    {
public:

    Vector v;

    ITDat() 
        : v(0), numref(0)
        { }

    explicit 
    ITDat(int size) 
        : v(size), numref(0)
        { assert(size > 0); v = 0; }

    explicit 
    ITDat(const VectorRef& v_) 
        : v(v_), numref(0)
        { }

    explicit 
    ITDat(Real r) 
        : v(1), numref(0)
        { v = r; }

    explicit 
    ITDat(std::istream& s) 
        : numref(0) 
        { read(s); }

    explicit 
    ITDat(const ITDat& other) 
        : v(other.v), numref(0)
        { }

    void 
    read(std::istream& s);

    void 
    write(std::ostream& s) const;
    
    void 
    print() const 
        { std::cout << "ITDat: v = " << v; }

    void* operator 
    new(size_t) throw(std::bad_alloc)
        { return allocator().alloc(); }

    void operator 
    delete(void* p) throw()
        { return allocator().dealloc(p); }

    friend class ITensor;

    ENABLE_INTRUSIVE_PTR(ITDat)

private:

    mutable unsigned int 
    numref;

    //Must be dynamically allocated:
    void operator=(const ITDat&);
    ~ITDat() { }

    static DatAllocator<ITDat>& allocator()
        {
        static DatAllocator<ITDat> allocator_;
        return allocator_;
        }

    };


class commaInit
    {
public:
    commaInit(ITensor& T_)
        : T(T_)
        { 
        if(T.isNull()) Error("Can't assign to null ITensor");
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

template <typename Callable> void ITensor::
mapElems(const Callable& f)
    {
    solo();
    scaleTo(1);
    for(int j = 1; j <= p->v.Length(); ++j)
        p->v(j) = f(p->v(j));
    }

Real 
Dot(const ITensor& x, const ITensor& y, bool doconj = true);

void 
Dot(const ITensor& x, const ITensor& y, Real& re, Real& im, 
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

#endif
