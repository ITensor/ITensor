//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITENSOR_H
#define __ITENSOR_ITENSOR_H
#include "itdata_functions.h"

namespace itensor {

//
// ITensor
//
class ITensor
    {
    public:

    using storage_type = ITData;
    using storage_ptr = PData;
    using IndexT = Index;
    using IndexValT = IndexVal;

    //
    // Constructors
    //

    //Construct Null ITensor, ITensor will evaluate to false in boolean context
    ITensor();

    //Construct rank 1 ITensor, all entries set to zero
    explicit
    ITensor(const Index& i1);

    //Construct rank 2 ITensor, all entries set to zero
    ITensor(const Index& i1,
            const Index& i2);

    //Construct rank n ITensor, all entries set to zero
    template <typename... Indices>
    ITensor(const Index& i1, 
            const Index& i2, 
            const Index& i3, 
            const Indices&... rest);

    //Construct rank 0 ITensor (scalar), value set to val
    explicit
    ITensor(Real val);
    explicit
    ITensor(Complex val);

    //Construct ITensor with diagonal elements set to z
    //(z can be Real argument too; convertible to Complex)
    template<typename... Inds>
    ITensor(Complex z, 
            const Index& i1,
            const Inds&... inds);

    //Construct diagonal ITensor,
    //diagonal elements given by VectorRef V
    template<typename... Inds>
    ITensor(const VectorRef& V,
            const Index& i1,
            const Inds&... inds);

    //Construct rank n ITensor, all
    //entries set to zero except the single
    //entry specified by the IndexVal args
    template <typename... IVals>
    explicit
    ITensor(const IndexVal& iv1, 
            const IVals&... rest);

    //
    // Accessor Methods
    //

    //Rank of this ITensor (number of indices)
    int 
    r() const { return is_.r(); }

    //Access index set
    const IndexSet&
    inds() const { return is_; }

    //evaluate to false if ITensor is default constructed
    explicit operator bool() const { return bool(store_); }

    template <typename... IndexVals>
    Real
    real(const IndexVals&... ivs) const;

    template <typename... IndexVals>
    Complex
    cplx(const IndexVals&... ivs) const;

    template<typename... IndexVals>
    void
    set(Real val, const IndexVals&... ivs);

    template<typename... IndexVals>
    void
    set(Complex val, const IndexVals&... ivs);

    //
    // Operators
    //

    //Contracting product
    //All matching Index pairs automatically contracted
    //Cji = \sum_{k,l} Akjl * Blki
    ITensor& 
    operator*=(const ITensor& other);

    // Contract with IndexVal
    // If iv = (J,n), Index J is fixed to it's nth
    // value and rank decreases by 1
    // (similar to summing against a Kronecker
    // delta tensor \delta_{J,n})
    ITensor& 
    operator*=(const IndexVal& iv) { return operator*=(ITensor(iv)); } 

    //Multiplication and division by scalar
    ITensor& 
    operator*=(Real fac);

    ITensor& 
    operator/=(Real fac);

    ITensor& 
    operator*=(Complex z);

    ITensor& 
    operator/=(Complex z);

    ITensor
    operator-() const;


    //Tensor addition and subtraction
    //Summands must have same Indices, in any order
    //Cijk = Aijk + Bkij
    ITensor& 
    operator+=(const ITensor& other);

    ITensor& 
    operator-=(const ITensor& other);

    //
    // Index Prime Level Methods
    //

    //Set primeLevel of Indices to zero
    ITensor& 
    noprime(IndexType type = All);

    //Set primeLevel of Index I to zero
    ITensor& 
    noprime(const Index& I);

    //Increase primeLevel of Indices by 1 (or optional amount inc)
    ITensor& 
    prime(int inc = 1);

    //Increase primeLevel of Indices by 1 (or optional amount inc)
    ITensor& 
    prime(IndexType type, int inc = 1);

    //Increase primeLevel of Index I by 1 (or optional amount inc)
    ITensor& 
    prime(const Index& I, int inc = 1);

    //Change all Indices having primeLevel plevold to have primeLevel plevnew
    ITensor& 
    mapprime(int plevold, int plevnew, IndexType type = All);

    //
    // Element Transformation Methods
    //

    ITensor&
    fill(Real r);

    ITensor&
    fill(Complex z);

    template <typename Func>
    ITensor&
    generate(Func&& f);

    template <typename Func>
    ITensor&
    apply(Func&& f);

    template <typename Func>
    const ITensor&
    visit(Func&& f) const;

    private:

    //Disattach self from current ITData and create own copy instead.
    //Necessary because ITensors logically represent distinct
    //tensors even though they may share data (copy-on-write idiom)
    void 
    solo();

    void
    scaleOutNorm();

    void
    equalizeScales(ITensor& other);

    public:

    //
    // Developer / advanced methods
    //
    // The following methods should not
    // be needed for most user code.
    //

    //Construct by explicitly providing data members
    ITensor(IndexSet&& iset,
            NewData nd,
            LogNumber scale);

    //Provide indices from IndexSet
    explicit
    ITensor(const IndexSet& is);

    //Provide indices from an index set
    //and elements from a VectorRef
    ITensor(const IndexSet& is,
            const VectorRef& v);

    ITensor(const IndexSet& is,
            const ITensor& t);

    //ITensor(const IndexSet& is,
    //        const ITensor& t,
    //        const Permutation& P);

    //Scale factor, used internally for efficient scalar ops.
    //Mostly for developer use; not necessary to explicitly involve
    //scale factors in user-level ITensor operations.
    const LogNumber&
    scale() const { return scale_; }

    const ITData&
    data() const { return *store_; }

    void 
    scaleTo(const LogNumber& newscale);

    //
    // Deprecated methods
    //

    //Construct matrix-like rank 2 ITensor,
    //elements given by MatrixRef M
    ITensor(const Index& i1,
            const Index& i2,
            const MatrixRef& M);

    Real
    norm() const;

    template <typename Callable> 
    ITensor&
    mapElems(Callable&& f)
        {
        return apply(std::forward<Callable>(f));
        }

    const IndexSet&
    indices() const { return inds(); }

    private:
    ////////////////////////
    IndexSet is_;
    LogNumber scale_;
    storage_ptr store_;
    ////////////////////////
    }; // class ITensor


template <typename... Indices>
ITensor::
ITensor(const Index& i1, 
        const Index& i2,
        const Index& i3,
        const Indices&... rest)
    :
    is_(i1,i2,i3,rest...),
    scale_(1.),
    store_(make_shared<ITDense<Real>>(area(is_),0.))
	{ }

template<typename... Inds>
ITensor::
ITensor(Complex z, 
        const Index& i1,
        const Inds&... inds)
    :
    is_(i1,inds...),
    scale_(1.)
    { 
    if(z.imag() == 0)
        store_ = make_shared<ITDiag<Real>>(z.real());
    else
        store_ = make_shared<ITDiag<Complex>>(z);
    }

template<typename... Inds>
ITensor::
ITensor(const VectorRef& V, 
        const Index& i1,
        const Inds&... inds)
    :
    is_(i1,inds...),
    scale_(1.),
    store_(std::make_shared<ITDiag<Real>>(V.begin(),V.end()))
    { 
#ifdef DEBUG
    //Compute min of all index dimensions
    long minm = i1.m();
    for(const auto& ind : is_)
        if(ind.m() < minm) minm = ind.m();
    if(V.Length() != minm)
        {
        Print(minm);
        Print(V.Length());
        Error("Wrong size of data in diagonal ITensor constructor");
        }
#endif
    }


template <typename... IVals>
ITensor::
ITensor(const IndexVal& iv1, 
        const IVals&... rest)
    :
    scale_(1.)
    {
    const auto size = 1+sizeof...(rest);
    auto ivs = std::array<IndexVal,size>{{iv1,rest...}};
    auto inds = std::vector<Index>();
    inds.reserve(size);
    for(const auto& iv : ivs) inds.push_back(iv.index);
    is_ = IndexSet(std::move(inds));

    store_ = make_shared<ITDense<Real>>(area(is_),0.);

    set(1.,iv1,rest...);
    }


template <typename... IndexVals>
Complex ITensor::
cplx(const IndexVals&... ivs) const
    {
#ifdef DEBUG
    if(!*this) Error("ITensor is default constructed");
#endif
    static constexpr auto size = sizeof...(ivs);
    std::array<IndexVal,size> vals{{static_cast<IndexVal>(ivs)...}};
    std::vector<long> inds(is_.size(),0);
    detail::permute_map(is_,vals,inds,[](const IndexVal& iv) { return iv.i-1; });
    auto g = applyFunc<GetElt<Complex,size>>(store_,{is_,inds});
	try {
	    return Complex(g)*scale_.real(); 
	    }
	catch(const TooBigForReal& e)
	    {
	    println("too big for real in cplx(...), scale = ",scale());
	    throw e;
	    }
	catch(TooSmallForReal)
	    {
        println("warning: too small for real in cplx(...)");
	    return Complex(0.,0.);
	    }
    return Complex(NAN,NAN);
    }


template <typename... IndexVals>
Real ITensor::
real(const IndexVals&... ivs) const
    {
    auto z = cplx(ivs...);
    if(fabs(z.imag()) != 0)
        {
        printfln("element = (%.5E,%.5E)",z.real(),z.imag());
        Error("ITensor is Complex-valued, use .cplx(...) method");
        }
    return z.real();
    }

template <typename... IndexVals>
void ITensor::
set(Real val, const IndexVals&... ivs)
    {
    static constexpr auto size = sizeof...(ivs);
    scaleTo(1.);
    const std::array<IndexVal,size> vals = {{ static_cast<IndexVal>(ivs)...}};
    std::array<int,size> inds;
    detail::permute_map(is_,vals,inds,[](const IndexVal& iv) { return iv.i-1; });
    applyFunc<SetEltReal<size>>(store_,{val,is_,inds});
    }

template <typename... IndexVals>
void ITensor::
set(Complex val, const IndexVals&... ivs)
    {
    static constexpr auto size = sizeof...(ivs);
    scaleTo(1.);
    const std::array<IndexVal,size> vals = {{ static_cast<IndexVal>(ivs)...}};
    std::array<int,size> inds;
    detail::permute_map(is_,vals,inds,[](const IndexVal& iv) { return iv.i-1; });
    applyFunc<SetEltComplex<size>>(store_,{val,is_,inds});
    }

template <typename Func>
ITensor& ITensor::
generate(Func&& f)
    {
    solo();
    scaleTo(1);
    applyFunc<GenerateIT<decltype(f)>>(store_,{std::forward<Func>(f)});
    return *this;
    }

template <typename Func>
ITensor& ITensor::
apply(Func&& f)
    {
    solo();
    scaleTo(1);
    applyFunc<ApplyIT<decltype(f)>>(store_,{std::forward<Func>(f)});
    return *this;
    }

template <typename Func>
const ITensor& ITensor::
visit(Func&& f) const
    {
    applyFunc<VisitIT<decltype(f)>>(store_,{std::forward<Func>(f),scale()});
    return *this;
    }

std::ostream& 
operator<<(std::ostream & s, const ITensor& T);

ITensor inline
operator*(ITensor A, const ITensor& B) { A *= B; return A; }
ITensor inline
operator*(ITensor T, Real fac) { T *= fac; return T; }
ITensor inline
operator*(Real fac, ITensor T) { T *= fac; return T; }
ITensor inline
operator*(ITensor T, Complex fac) { T *= fac; return T; }
ITensor inline
operator*(Complex fac, ITensor T) { T *= fac; return T; }
ITensor inline
operator/(ITensor T, Real fac) { T /= fac; return T; }
ITensor inline
operator/(ITensor T, Complex fac) { T /= fac; return T; }
ITensor inline
operator+(ITensor A, const ITensor& B) { A += B; return A; }
ITensor inline
operator-(ITensor A, const ITensor& B) { A -= B; return A; }

ITensor inline
operator*(ITensor T, const IndexVal& iv) { T *= iv; return T; }
ITensor inline
operator*(const IndexVal& iv, const ITensor& t) { return (ITensor(iv) *= t); }

ITensor
combiner(std::vector<Index> inds);

template<typename... Inds>
ITensor
combiner(const Index& i1, const Inds&... inds)
    {
    return combiner({i1,inds...});
    }

ITensor
delta(const Index& i1, const Index& i2);

//
// Define product of IndexVal iv1 = (I1,n1), iv2 = (I2,n2)
// (I1, I2 are Index objects; n1,n2 are type int)
// to be an ITensor T such that T(I1(n1),I2(n2)) == 1
//
// Useful for creating MPOs
//
ITensor inline
operator*(const IndexVal& iv1, const IndexVal& iv2) 
    { 
    ITensor t(iv1); 
    return (t *= iv2); 
    }
//
// Define product of IndexVal iv1 = (I1,n1) with a Real "val"
// to be an ITensor T such that T(I1(n1)) == val
//
// Useful for creating MPOs
//
ITensor inline
operator*(const IndexVal& iv1, Real val) 
    { 
    ITensor res(iv1); 
    res *= val; 
    return res; 
    }
ITensor inline
operator*(Real val, const IndexVal& iv) { return operator*(iv,val); }

//Return copy of ITensor with primeLevel of Index I increased by 1
//(or optional amount inc)
template <class Tensor, class IndexT>
Tensor
prime(Tensor A, const IndexT& I, int inc = 1)
    { 
    A.prime(I,inc); 
    return A; 
    }

//Return copy of ITensor with primeLevel of Index I set to zero
template <class Tensor, class IndexT>
Tensor
noprime(Tensor A, const IndexT& I)
    { 
    A.noprime(I); 
    return A; 
    }

template<class Tensor>
bool
hasindex(const Tensor& T, const typename Tensor::IndexT& I)
    {
    return detail::contains(T.inds(),I);
    }

ITensor
randIT(ITensor T, const Args& args = Global::args());

template <typename... Indices>
ITensor
randIT(const Index& i1, const Indices&... rest)
    {
    return randIT(ITensor(i1,rest...));
    }

template <typename... Indices>
ITensor
tieIndex(const ITensor& T,
         const Index& t0,
         const Index& t1,
         const Indices&... rest);

//Compute the norm of an ITensor.
//Thinking of elements as a vector, equivalent to sqrt(v*v).
//Result is equivalent to sqrt((T*T).real()) 
//[and similar for complex case] but computed much more efficiently.
Real 
norm(const ITensor& T);

ITensor
conj(const ITensor& T);

ITensor inline
dag(const ITensor& T) { return conj(T); }

bool
isComplex(const ITensor& T);

Real
sumels(const ITensor& T);

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
          IndexType type = All);

//Find index of tensor A (of optional type t) 
//which is shared with tensor B
template<class TensorA, class TensorB> typename 
TensorA::IndexT
commonIndex(const TensorA& A, const TensorB& B, IndexType t = All);

//Find index of tensor A (of optional type t) 
//which is NOT shared by tensor B
template<class TensorA, class TensorB> typename 
TensorA::IndexT
uniqueIndex(const TensorA& A, 
            const TensorB& B, 
            IndexType t);






//
//
//
// Template Method Implementations
//
//
//

template <class Tensor>
Tensor
swapPrime(Tensor T, int plev1, int plev2,
          IndexType type)
    { 
    int tempLevel = 100;
#ifdef DEBUG
    for(const auto& I : T.inds())
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

template<class TensorA, class TensorB> typename 
TensorA::IndexT
commonIndex(const TensorA& A, const TensorB& B, IndexType t)
    {
    using IndexT = typename TensorA::IndexT;
    for(const IndexT& I : A.inds())
        {
        if( (t == All || I.type() == t)
         && hasindex(B.inds(),I) ) 
            {
            return I;
            }
        }
    return IndexT();
    }


template<class TensorA, class TensorB> typename 
TensorA::IndexT
uniqueIndex(const TensorA& A, 
            const TensorB& B, 
            IndexType t)
    {
    using IndexT = typename TensorA::IndexT;
    for(const IndexT& I : A.inds())
        {
        if( (t == All || I.type() == t)
         && !hasindex(B.inds(),I) ) 
            {
            return I;
            }
        }
    return IndexT::Null();
    }



//template <typename... Indices>
//ITensor
//tieIndex(const ITensor& T,
//         const Index& t0,
//         const Index& t1,
//         const Indices&... rest)
//    {
//    static constexpr auto size = 2 + sizeof...(rest);
//    if(size > T.r()) Error("Cannot tie more indices than ITensor rank.");
//    std::array<Index,size> totie = {{ t0, t1, static_cast<Index>(rest)...}};
//    std::array<size_t,size> I;
//    NewIndexSet<Index> new_index(T.r()-size+1);
//    size_t nt = 0;
//    for(int j = 0; j < T.r(); ++j)
//        {
//        const auto& J = T.inds()[j];
//        if(detail::contains(totie,J))
//            {
//            if(J == totie.front()) new_index.add(J);
//            I[nt++] = j;
//            }
//        else
//            {
//            new_index.add(J);
//            }
//        }
//    if(nt != totie.size())
//        Error("ITensor does not have requested Index to tie");
//
//    auto nd = T.data().clone();
//    const auto f = [&I](const btas::Range& r) { return tieIndex(r,I); };
//    applyFunc<ApplyRange<decltype(f)>>(nd,{f});
//
//    return ITensor(new_index,std::move(nd),T.scale());
//    }

}; //namespace itensor


#endif
