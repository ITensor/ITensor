//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITENSOR_H
#define __ITENSOR_ITENSOR_H
#include "real.h"
#include "counter.h"

namespace itensor {

//Forward declarations
struct ProductProps;
class Combiner;
class ITDat;
class ITensor;

void toMatrixProd(const ITensor& L, const ITensor& R, 
                  ProductProps& pp,
                  MatrixRefNoLink& lref, MatrixRefNoLink& rref,
                  bool& L_is_matrix, bool& R_is_matrix, bool doReshape = true);

//
// ITensor
//
class ITensor : public safe_bool<ITensor>
    {
    public:

    //
    //Accessor Methods
    //

    //Rank of this ITensor (number of indices)
    int 
    r() const { return is_.r(); }

    //false if ITensor is default constructed
    bool 
    valid() const { return bool(r_); }

    bool
    isComplex() const { return bool(i_); }

    //Sets this ITensor to its real part only
    ITensor&
    takeRealPart();

    //Sets this ITensor to its imaginary part only
    //(will be real afterward)
    ITensor&
    takeImagPart();

    //Enables looping over Indices in a Foreach
    //e.g. Foreach(const Index& I, t.index() ) { ... }
    const IndexSet<Index>&
    indices() const { return is_; }

    //Read-only access to scale factor, used internally for efficient scalar ops
    const LogNumber&
    scale() const { return scale_; }

    enum Type { Null, Dense, Diag };

    Type
    type() const { return type_; }

    //Return diagonal components
    Vector
    diag() const;

    //
    //Constructors
    //

    ITensor();

    //Construct rank 1 ITensor, all entries set to zero
    explicit 
    ITensor(const Index& i1);

    //Construct rank 2 ITensor, all entries set to zero
    ITensor(const Index& i1,const Index& i2);

    //Construct ITensor up to rank 8, entries set to zero
    ITensor(const Index& i1, const Index& i2, const Index& i3,
            const Index& i4 = Index::Null(),
            const Index& i5 = Index::Null(),
            const Index& i6 = Index::Null(),
            const Index& i7 = Index::Null(),
            const Index& i8 = Index::Null());

    //Construct rank 0 ITensor (scalar), value set to val
    explicit
    ITensor(Real val);

    //Construct rank 1 ITensor, all entries set to val
    ITensor(const Index& i1, Real val);

    //Construct rank 1 ITensor, entries set to those of V
    ITensor(const Index& i1, const VectorRef& V);

    //Construct rank 2 ITensor, entries set to those of M
    ITensor(const Index& i1, const Index& i2, const MatrixRef& M);

    //Construct rank 2 ITensor (a matrix) with 'a' on the diagonal
    ITensor(const Index& i1, const Index& i2, Real a);

    //Construct rank 2 ITensor (matrix), diagonal entries set to those of V
    ITensor(const Index& i1, const Index& i2, const VectorRef& V);

    // Construct rank 1 tensor T from IndexVal iv = (I,n)
    // (I is an Index, n an int)
    // such that T(I(n)) == 1
    explicit 
    ITensor(const IndexVal& iv);

    // Construct rank 2 tensor T from IndexVals 
    // iv1 = (I1,n1), iv2 = (I2,n2)
    // such that T(I1(n1),I2(n2)) == 1
    ITensor(const IndexVal& iv1, const IndexVal& iv2);

    // Construct tensor T from up to 8 IndexVals 
    // iv1 = (I1,n1), iv2 = (I2,n2), iv3 = (I3,n3), ...
    // such that T(I1(n1),I2(n2),I3(n3),...) == 1
    ITensor(const IndexVal& iv1, const IndexVal& iv2, 
            const IndexVal& iv3, const IndexVal& iv4 = IndexVal::Null(), 
            const IndexVal& iv5 = IndexVal::Null(), const IndexVal& iv6 = IndexVal::Null(), 
            const IndexVal& iv7 = IndexVal::Null(), const IndexVal& iv8 = IndexVal::Null());

    explicit 
    ITensor(const IndexSet<Index>& I);

    ITensor(const IndexSet<Index>& I, const Vector& V);

    ITensor(const IndexSet<Index>& I, const ITensor& other);

    ITensor(const IndexSet<Index>& I, const ITensor& other, 
            const Permutation& P);

    //Read in ITensor from binary stream s
    void 
    read(std::istream& s);

    //Write out ITensor to binary stream s
    void
    write(std::ostream& s) const;


    //
    // Operators
    //

    //Contracting product
    //All matching Index pairs automatically contracted
    //Cji = \sum_{k,l} Akjl * Blki
    ITensor& 
    operator*=(const ITensor& other);

    //Non-contracting product
    //Matching Index pairs are merged
    //Ckjli = Akjl * Blki
    ITensor& 
    operator/=(const ITensor& other);

    //Multiplication and division by scalar
    ITensor& 
    operator*=(Real fac);

    ITensor& 
    operator/=(Real fac) { scale_ /= fac; return *this; }

    ITensor
    operator-() const { ITensor T(*this); T.scale_ *= -1; return T; }

    //Multiplication by Complex scalar
    ITensor&
    operator*=(Complex z);

    //Multiplication with LogNumber (very large or very small Real)
    ITensor& 
    operator*=(const LogNumber& lgnum) { scale_ *= lgnum; return *this; }

    // Contract with IndexVal
    // If iv = (J,n), Index J is fixed to it's nth
    // value and rank decreases by 1
    // (similar to summing against a Kronecker
    // delta tensor \delta_{J,n})
    ITensor& 
    operator*=(const IndexVal& iv) { return operator*=(ITensor(iv)); } 

    //Tensor addition and subtraction
    //Summands must have same Indices, in any order
    //Cijk = Aijk + Bkij
    ITensor& 
    operator+=(const ITensor& other);

    ITensor& 
    operator-=(const ITensor& other);


    //
    //Primelevel Methods
    //

    //Set primeLevel of Indices to zero
    ITensor& 
    noprime(IndexType type = All) { is_.noprime(type); return *this; }

    //Set primeLevel of Index I to zero
    ITensor& 
    noprime(const Index& I) { is_.noprime(I); return *this; }

    //Increase primeLevel of Indices by 1 (or optional amount inc)
    ITensor& 
    prime(int inc = 1) { prime(All,inc); return *this;}

    //Increase primeLevel of Indices by 1 (or optional amount inc)
    ITensor& 
    prime(IndexType type, int inc = 1) { is_.prime(type,inc); return *this; }

    //Increase primeLevel of Index I by 1 (or optional amount inc)
    ITensor& 
    prime(const Index& I, int inc = 1) { is_.prime(I,inc); return *this; }

    //Change all Indices having primeLevel plevold to have primeLevel plevnew
    ITensor& 
    mapprime(int plevold, int plevnew, IndexType type = All)
        { is_.mapprime(plevold,plevnew,type); return *this; }

    //
    //Element Access Methods
    //

    //Get scalar value of rank 0 ITensor
    //Throws ITError if r() != 0
    Real
    toReal() const;

    //Get scalar value of rank 0 ITensor
    //Throws ITError if r() != 0
    Complex
    toComplex() const;

    // IndexVal element access
    // Given iv1 = (I1,n1), iv2 = (I2,n2), ...
    // returns component of ITensor such that
    // I1 temporarily set to n1, I2 to n2, etc.
    // Can be used to set components of ITensors
    // as well, for example, T(I1(2),I2(1)) = 3;
    Real& 
    operator()(const IndexVal& iv1);

    Real 
    operator()(const IndexVal& iv1) const;

    Real& 
    operator()(const IndexVal& iv1, const IndexVal& iv2);

    Real 
    operator()(const IndexVal& iv1, const IndexVal& iv2) const;

    Real& 
    operator()(const IndexVal& iv1, const IndexVal& iv2, const IndexVal& iv3, 
               const IndexVal& iv4 = IndexVal::Null(), 
               const IndexVal& iv5 = IndexVal::Null(),
               const IndexVal& iv6 = IndexVal::Null(),
               const IndexVal& iv7 = IndexVal::Null(),
               const IndexVal& iv8 = IndexVal::Null());

    Real 
    operator()(const IndexVal& iv1, const IndexVal& iv2, const IndexVal& iv3, 
               const IndexVal& iv4 = IndexVal::Null(), 
               const IndexVal& iv5 = IndexVal::Null(),
               const IndexVal& iv6 = IndexVal::Null(),
               const IndexVal& iv7 = IndexVal::Null(),
               const IndexVal& iv8 = IndexVal::Null()) const;

    //
    //Methods for Mapping to Other Objects
    //

    //
    // groupIndices combines a set of indices (of possibly different sizes) 
    // together, leaving only single grouped Index.
    //
    // RiJ = Ai(jk) <-- Here J represents the grouped pair of indices (jk)
    //                  If j.m() == 5 and k.m() == 7, J.m() == 5*7.
    //
    void 
    groupIndices(const array<Index,NMAX+1>& indices, int nind, 
                      const Index& grouped, ITensor& res) const;

    //
    // tieIndices locks a set of indices (of the same size) together,
    // leaving only a single tied Index.
    //
    // Rijl = Aijil <-- here we have tied the 1st and 3rd index of A
    //
    void
    tieIndices(const array<Index,NMAX>& indices, int nind,
               const Index& tied);

    void
    tieIndices(const Index& i1, const Index& i2,
               const Index& tied);


    void
    tieIndices(const Index& i1, const Index& i2,
               const Index& i3,
               const Index& tied);


    void
    tieIndices(const Index& i1, const Index& i2,
               const Index& i3, const Index& i4,
               const Index& tied);


    // The trace method sums over the given set of indices
    // (which must all have the same dimension).
    //
    // Rik = Aijkml.trace(j,l,m) = \sum_t Aitktt
    ITensor&
    trace(const Index& i1, 
          const Index& i2 = Index::Null(), 
          const Index& i3 = Index::Null(),
          const Index& i4 = Index::Null(),
          const Index& i5 = Index::Null(),
          const Index& i6 = Index::Null(),
          const Index& i7 = Index::Null(),
          const Index& i8 = Index::Null());

    ITensor&
    trace(const array<Index,NMAX>& indices, int nind = -1);


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

    //Set components of rank 2 ITensor using Matrix M as input
    void 
    fromMatrix11(const Index& i1, const Index& i2, const Matrix& M);

    //Convert rank 2 ITensor to a Matrix using given Index order
    void 
    toMatrix11(const Index& i1, const Index& i2, Matrix& res) const;

    //Convert rank 2 ITensor to a Matrix, but do not include
    //scale factor in result
    void 
    toMatrix11NoScale(const Index& i1, const Index& i2, 
                           Matrix& res) const;

    // group i1,i2; i3,i4
    void 
    toMatrix22(const Index& i1, const Index& i2, 
               const Index& i3, const Index& i4, Matrix& res) const;
    void 
    fromMatrix22(const Index& i1, const Index& i2, 
                 const Index& i3, const Index& i4, const Matrix& res);

    /*
    // group i1,i2; i3
    void toMatrix21(const Index& i1, const Index& i2, 
                    const Index& i3, Matrix& res) const;
    void fromMatrix21(const Index& i1, const Index& i2, 
                      const Index& i3, const Matrix& res);

    */

    // group i1; i2,i3
    void toMatrix12NoScale(const Index& i1, const Index& i2, 
                           const Index& i3, Matrix& res) const;

    void toMatrix12(const Index& i1, const Index& i2, 
                    const Index& i3, Matrix& res) const;

    void fromMatrix12(const Index& i1, const Index& i2, 
                      const Index& i3, const Matrix& M);



    //
    // Swap can be used for similar purposes
    // as operator=(const ITensor& other)
    // but is more efficient and has same
    // end result if other is just a temporary
    //
    void
    swap(ITensor& other);


    //Other Methods -------------------------------------------------

    const Real*
    datStart() const;

    const Real*
    imagDatStart() const;

    void 
    randomize(const OptSet& opts = Global::opts());

    void 
    conj();

    void 
    dag() { conj(); }

    Real 
    norm() const;

    LogNumber 
    normLogNum() const;

    Real 
    normNoScale() const;

    template <typename Callable> 
    ITensor&
    mapElems(const Callable& f);

    void 
    scaleOutNorm();

    void 
    scaleTo(const LogNumber& newscale);

    VectorRef 
    assignToVec() const;

    void
    pseudoInvert(Real cutoff = 0.);

    void
    replaceIndex(const Index& oind,
                 const Index& nind)
        { is_.replaceIndex(oind,nind); }

    //
    // Typedefs
    //

    typedef Index 
    IndexT;

    typedef IndexVal 
    IndexValT;

    typedef Combiner 
    CombinerT;

    //Deprecated methods --------------------------

    //Use indices().dim() instead of vecSize
    //int 
    //vecSize() const;

    //void 
    //assignFromVec(const VectorRef& v);


    private:

    //////////////
    //
    // Data Members
    //

    Type type_;

    //Pointer to ITDat containing tensor data
    shared_ptr<ITDat> r_, //real part
                      i_; //imag part

    //Indices, maximum of 8
    IndexSet<Index> is_;

    //scale_ absorbs scalar factors to avoid copying ITDat
    LogNumber scale_; 

    //
    //
    //////////////

    void 
    allocate(int dim);
    void 
    allocate();

    void 
    allocateImag(int dim);
    void 
    allocateImag();

    //Disattach self from current ITDat and create own copy instead.
    //Necessary because ITensors logically represent distinct
    //objects even though they may share data
    void 
    solo();

    void 
    soloReal();

    void 
    soloImag();

    void
    equalizeScales(ITensor& other);

    void
    reshapeDat(const Permutation& P);
    
    friend struct ProductProps;

    friend void toMatrixProd(const ITensor& L, const ITensor& R, 
                             ProductProps& pp,
                             MatrixRefNoLink& lref, MatrixRefNoLink& rref,
                             bool& L_is_matrix, bool& R_is_matrix, bool doReshape);

    int _ind2(const IndexVal& iv1, const IndexVal& iv2) const;

    int 
    _ind8(const IndexVal& iv1, const IndexVal& iv2, 
          const IndexVal& iv3, const IndexVal& iv4 = IndexVal::Null(), 
          const IndexVal& iv5 = IndexVal::Null(),const IndexVal& iv6 = IndexVal::Null(),
          const IndexVal& iv7 = IndexVal::Null(),const IndexVal& iv8 = IndexVal::Null())
        const;
    int 
    _diag_ind8(const IndexVal& iv1, const IndexVal& iv2, 
               const IndexVal& iv3, const IndexVal& iv4 = IndexVal::Null(), 
               const IndexVal& iv5 = IndexVal::Null(),const IndexVal& iv6 = IndexVal::Null(),
               const IndexVal& iv7 = IndexVal::Null(),const IndexVal& iv8 = IndexVal::Null())
        const;

    void
    convertToDense();

    friend class commaInit;

    friend void 
    contractDiagDense(const ITensor& S, const ITensor& T, ITensor& res);

    friend void 
    contractDiagDiag(const ITensor& A, const ITensor& B, ITensor& res);


    friend std::ostream& 
    operator<<(std::ostream & s, const ITensor& T);

    }; // class ITensor



class commaInit
    {
    public:

    commaInit(ITensor& T,
              const Index& i1,
              const Index& i2 = Index::Null(),
              const Index& i3 = Index::Null());

    commaInit& 
    operator=(Real r);

    commaInit& 
    operator<<(Real r) { return operator=(r); }

    commaInit& 
    operator,(Real r);

    ~commaInit();

    private:

    ITensor& T_;
    bool started_;
    Counter c_; 
    Permutation P_;

    };

//
// ITDat
//
class ITDat
    {
    public:

    Vector v;

    ITDat();

    explicit 
    ITDat(int size);

    explicit 
    ITDat(const VectorRef& v_);

    explicit 
    ITDat(Real r);

    explicit 
    ITDat(const ITDat& other);

    int
    size() const { return v.Length(); }

    void
    read(std::istream& s);

    void 
    write(std::ostream& s) const;
    
    friend class ITensor;

    private:

    //Must be dynamically allocated:
    void operator=(const ITDat&);


    };

ITensor inline
operator*(ITensor A, const ITensor& B) { A *= B; return A; }

ITensor inline
operator*(ITensor T, const IndexVal& iv) { T *= iv; return T; }

ITensor inline
operator*(const IndexVal& iv, const ITensor& t) { return (ITensor(iv) *= t); }

ITensor inline
operator*(ITensor T, Real fac) { T *= fac; return T; }

ITensor inline
operator*(Real fac, ITensor T) { T *= fac; return T; }

ITensor inline
operator/(ITensor T, Real fac) { T /= fac; return T; }

ITensor inline
operator*(ITensor T, Complex z) { T *= z; return T; }

ITensor inline
operator*(Complex z, ITensor T) { T *= z; return T; }

ITensor inline
operator*(ITensor T, LogNumber lgnum) { T *= lgnum; return T; }

ITensor inline
operator*(LogNumber lgnum, ITensor T) { T *= lgnum; return T; }

ITensor inline
operator/(ITensor A, const ITensor& B) { A /= B; return A; }

ITensor inline
operator+(ITensor A, const ITensor& B) { A += B; return A; }

ITensor inline
operator-(ITensor A, const ITensor& B) { A -= B; return A; }

template <class Tensor, class IndexT>
bool inline
operator==(const IndexSet<IndexT>& is, const Tensor& t)
    { return is == t.indices(); }

template <class Tensor, class IndexT>
bool inline
operator==(const Tensor& t, const IndexSet<IndexT>& is)
    { return t.indices() == is; }

template <class Tensor, class IndexT>
bool inline
operator<(const Tensor& t, const IndexSet<IndexT>& is)
    { return t.indices() < is; }

template <typename Callable> 
ITensor& ITensor::
mapElems(const Callable& f)
    {
    if(isComplex())
        throw ITError("mapElems only works for real ITensor");
    solo();
    scaleTo(1);
    for(int j = 1; j <= r_->size(); ++j)
        r_->v(j) = f(r_->v(j));
    return *this;
    }

Real
sumels(const ITensor& t);

ITensor inline
conj(ITensor res) { res.conj(); return res; }

ITensor inline
dag(ITensor res) { res.dag(); return res; }

//
// Computes the scalar/inner/dot product of two
// real-valued ITensors.
//
// Equivalent to the ITensor contraction x * y 
// except the result is a Real number instead
// of a rank 0 ITensor.
//
Real 
Dot(const ITensor& x, const ITensor& y);

//
// Scalar (inner) product of two
// possibly complex ITensors.
//
// Conjugates the first argument, therefore
// equivalent to the contraction dag(x) * y 
// (except it yields two real numbers, re and im,
// instead of a rank 0 ITensor).
//
Complex 
BraKet(const ITensor& x, const ITensor& y);

//
// Define product of IndexVal iv1 = (I1,n1), iv2 = (I2,n2)
// (I1, I2 are Index objects; n1,n2 are type int)
// to be an ITensor T such that T(I1(n1),I2(n2)) == 1
//
// Useful for creating MPOs
//
ITensor inline
operator*(const IndexVal& iv1, const IndexVal& iv2) 
    { ITensor t(iv1); return (t *= iv2); }

//
// Define product of IndexVal iv1 = (I1,n1) with a Real "val"
// to be an ITensor T such that T(I1(n1)) == val
//
// Useful for creating MPOs
//
ITensor inline
operator*(const IndexVal& iv1, Real val) 
    { ITensor res(iv1); res *= val; return res; }

ITensor inline
operator*(Real val, const IndexVal& iv) 
    { ITensor res(iv); res *= val; return res; }


//Find index of tensor A (of optional type t) 
//which is shared with tensor B
template<class TensorA, class TensorB> typename 
TensorA::IndexT
commonIndex(const TensorA& A, const TensorB& B, IndexType t = All)
    {
    typedef typename TensorA::IndexT
    IndexT;
    Foreach(const IndexT& I, A.indices())
        {
        if( (t == All || I.type() == t)
         && hasindex(B.indices(),I) ) 
            {
            return I;
            }
        }
    return IndexT::Null();
    }


//Find index of tensor A (of optional type t) 
//which is NOT shared by tensor B
template<class TensorA, class TensorB> typename 
TensorA::IndexT
uniqueIndex(const TensorA& A, 
            const TensorB& B, 
            IndexType t)
    {
    typedef typename TensorA::IndexT
    IndexT;
    Foreach(const IndexT& I, A.indices())
        {
        if( (t == All || I.type() == t)
         && !hasindex(B.indices(),I) ) 
            {
            return I;
            }
        }
    return IndexT::Null();
    }


template<class Tensor> typename
Tensor::IndexT const&
finddir(const Tensor& T, Arrow dir)
    {
    return finddir(T.indices(),dir);
    }

template<class Tensor> typename
Tensor::IndexT const&
findtype(const Tensor& T, IndexType type)
    {
    return findtype(T.indices(),type);
    }

template<class Tensor>
bool
hasindex(const Tensor& T, const typename Tensor::IndexT& I)
    {
    return hasindex(T.indices(),I);
    }

//
// Given Tensors which represent operator matrices
// (e.g. A(site1',site1), B(site1',site1) )
// multiply them, automatically adjusting primeLevels
// so that result is again an operator matrix C(site1',site1)
//
//              s'  t'
//  s'  t'      |   |
//  |   |       [-A-]
//  [-C-]  =    |   |
//  |   |       [-B-]
//  s   t       |   |
//              s   t
//
// (here s and t are indices of type Site)
//
template<class Tensor>
Tensor
multSiteOps(Tensor A, const Tensor& B) 
    {
    A.prime(Site);
    A *= B;
    A.mapprime(2,1,Site);
    return A;
    }


//Return copy of ITensor with primeLevel of Index I increased by 1
//(or optional amount inc)
template <class Tensor, class IndexT>
Tensor
prime(Tensor A, const IndexT& I, int inc = 1)
    { A.prime(I,inc); return A; }

//Return copy of ITensor with primeLevel of Index I set to zero
template <class Tensor, class IndexT>
Tensor
noprime(Tensor A, const IndexT& I)
    { A.noprime(I); return A; }

//
//Return copy of a tensor with primeLevels plev1 and plev2 swapped
//
//For example, if T has indices i,i' (like a matrix or a site
//operator) then swapPrime(T,0,1) will have indices i',i 
//i.e. the transpose of T.
//
template <class Tensor>
Tensor
swapPrime(Tensor T, int plev1, int plev2,
          IndexType type = All)
    { 
    const int tempLevel = 100;
#ifdef DEBUG
    Foreach(const typename Tensor::IndexT& I, T.indices())
        {
        if(I.primeLevel() == tempLevel) 
            {
            Print(tempLevel);
            Error("swapPrime fails if an index has primeLevel==tempLevel");
            }
        }
#endif
    T.mapprime(plev1,tempLevel,type);
    T.mapprime(plev2,plev1,type);
    T.mapprime(tempLevel,plev2,type);
    return T; 
    }

template <class Tensor, class IndexT>
Tensor
tieIndices(Tensor T,
           const IndexT& i1, const IndexT& i2, 
           const IndexT& tied)
    { 
    T.tieIndices(i1,i2,tied); 
    return T; 
    }

template <class Tensor, class IndexT>
Tensor
tieIndices(Tensor T,
           const IndexT& i1, const IndexT& i2, 
           const IndexT& i3, 
           const IndexT& tied)
    { 
    T.tieIndices(i1,i2,i3,tied); 
    return T; 
    }

template <class Tensor, class IndexT>
Tensor
tieIndices(Tensor T,
           const IndexT& i1, const IndexT& i2, 
           const IndexT& i3, const IndexT& i4, 
           const IndexT& tied)
    { 
    T.tieIndices(i1,i2,i3,i4,tied); 
    return T; 
    }

template <class Tensor>
Tensor
realPart(const Tensor& T)
    {
    if(!T.isComplex())
        return T;
    //else
    Tensor re(T);
    re.takeRealPart();
    return re;
    }

template <class Tensor>
Tensor
imagPart(const Tensor& T)
    {
    Tensor im(T);
    im.takeImagPart();
    return im;
    }

//Returns true if T is exactly zero.
//
//If passed the option Opt("Fast",true),
//only performs fast operations such as checking
//the scale of T, but skips computing T's norm.
//This can cause the return value to be true even
//if T is actually zero.
bool
isZero(const ITensor& T, const OptSet& opts = Global::opts());

//
// Tracing over all indices results in a Real
//
template <class Tensor>
Real
trace(Tensor T)
    {
    if(T.isComplex())
        {
        Error("ITensor is complex, use trace(T,re,im)");
        }
    if(T.indices().rn() != 0) 
        {
        T.trace(T.indices(),T.indices().rn());
        }
    return T.toReal();
    }

template<class Tensor>
void
trace(const Tensor& T, Real& re, Real& im)
    {
    if(!T.isComplex())
        {
        re = trace(T);
        im = 0;
        }
    else
        {
        re = trace(realPart(T));
        im = trace(imagPart(T));
        }
    }

template<class Tensor, class IndexT>
Tensor
trace(Tensor T, 
      const IndexT& i1,
      const IndexT& i2 = IndexT::Null(), 
      const IndexT& i3 = IndexT::Null(), 
      const IndexT& i4 = IndexT::Null(),
      const IndexT& i5 = IndexT::Null(),
      const IndexT& i6 = IndexT::Null(),
      const IndexT& i7 = IndexT::Null(),
      const IndexT& i8 = IndexT::Null())
    { 
    T.trace(i1,i2,i3,i4,i5,i6,i7,i8); 
    return T; 
    }

std::ostream& 
operator<<(std::ostream & s, const ITensor& T);


//
// Deprecated older functions.
// For backwards compatibility only.
//

//Return copy of ITensor with primeLevel of Index I increased by 1
//(or optional amount inc)
template <class Tensor, class IndexT>
Tensor
primed(Tensor A, const IndexT& I, int inc = 1)
    { A.prime(I,inc); return A; }

//Return copy of ITensor with primeLevel of Index I set to zero
template <class Tensor, class IndexT>
Tensor
deprimed(Tensor A, const IndexT& I)
    { A.noprime(I); return A; }

}; //namespace itensor

#endif
