//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MATRIX_VEC_H_
#define __ITENSOR_MATRIX_VEC_H_

#include "itensor/tensor/vecrange.h"
#include "itensor/tensor/ten.h"

namespace itensor {

using Vector     = Ten<VecRange>;
using VectorRef  = TenRef<VecRange>;
using VectorRefc = TenRefc<VecRange>;

auto inline
stride(VectorRefc const& v) { return v.stride(0); }

auto inline
stride(Vector const& v) { return v.stride(0); }

VectorRef
operator*=(VectorRef v, Real fac);

VectorRef
operator/=(VectorRef v, Real fac);

VectorRef
operator+=(VectorRef a, VectorRefc b);

VectorRef
operator-=(VectorRef a, VectorRefc b);

//Copy data referenced by b to memory referenced by a
VectorRef
operator&=(VectorRef a, VectorRefc b);

//Dot product
Real
operator*(VectorRefc a, VectorRefc b);

inline Vector&
operator*=(Vector & v, Real fac) { makeRef(v) *= fac; return v; }

inline Vector&
operator/=(Vector & v, Real fac) { makeRef(v) /= fac; return v; }

inline Vector& 
operator+=(Vector & v, Vector const& other) { makeRef(v) += makeRef(other); return v; }

inline Vector&
operator-=(Vector & v, Vector const& other) { makeRef(v) -= makeRef(other); return v; }

inline Vector&
operator+=(Vector & v, VectorRefc other) {  makeRef(v) += other; return v; }

inline Vector&
operator-=(Vector & v, VectorRefc other) { makeRef(v) -= other; return v; }

Vector inline
operator*(Vector A, Real fac) { A *= fac; return A; }

Vector inline
operator*(Real fac, Vector A) { A *= fac; return A; }

Vector inline
operator/(Vector A, Real fac) { A /= fac; return A; }

Vector inline
operator+(VectorRefc A, VectorRefc B)
    { 
    Vector res(A);
    res += B;
    return res;
    }

Vector inline
operator+(VectorRefc A, Vector&& B) 
    { 
    Vector res(std::move(B));
    res += A;
    return res;
    }

Vector inline
operator+(Vector&& A, VectorRefc B) 
    { 
    Vector res(std::move(A));
    res += B;
    return res;
    }

Vector inline
operator-(VectorRefc A, VectorRefc B)
    { 
    Vector res(A);
    res -= B;
    return res;
    }

Vector inline
operator-(VectorRefc A, Vector&& B) 
    { 
    Vector res(std::move(B)); 
    res *= -1;
    res += A; 
    return res; 
    }

Vector inline
operator-(Vector&& A, VectorRefc B) 
    { 
    Vector res(std::move(A)); 
    res -= B; 
    return res; 
    }

template<>
Real
norm(VectorRefc const& v);

void
randomize(VectorRef v);

Vector
randomVec(long size);

Real
sumels(VectorRefc v);

void 
resize(Vector & v, size_t newsize);


//
// These versions of op-assign to VectorRef can
// work safely for temporary Vectors since
// const references extend lifetime of rvalues
//

void inline
operator&=(VectorRef a, Vector const& b) { a &= makeRef(b); }
                                 
void inline                      
operator+=(VectorRef a, Vector const& b) { a += makeRef(b); }
                                 
void inline                      
operator-=(VectorRef a, Vector const& b) { a -= makeRef(b); }


template<>
std::ostream&
operator<<(std::ostream& s, VectorRefc const& v);

template<> inline
std::ostream&
operator<<(std::ostream& s, VectorRef const& v) { return operator<<(s,makeRefc(v)); }

inline std::ostream&
operator<<(std::ostream& s, Vector const& v) { return operator<<(s,makeRefc(v)); }

//
// makeVecRef functions
//

auto inline
makeVecRef(Real* p,
           size_t size)
    {
    return VectorRef({p,size},VecRange(size));
    }

auto inline
makeVecRef(const Real* p,
           size_t size)
    {
    return VectorRefc({p,size},VecRange(size));
    }

auto inline
makeVecRefc(const Real* p,
            size_t size)
    {
    return makeVecRef(p,size);
    }

auto inline
makeVecRef(Real* p,
           size_t size,
           size_t stride)
    {
    return VectorRef({p,size*stride},VecRange(size,stride));
    }

auto inline
makeVecRef(const Real* p,
           size_t size,
           size_t stride)
    {
    return VectorRefc({p,size*stride},VecRange(size,stride));
    }

auto inline
makeVecRefc(const Real* p,
            size_t size,
            size_t stride)
    {
    return makeVecRef(p,size);
    }

//
// Vector slicing operations
//

//Return ref to elements [start,stop], inclusive, of a vector
template<typename Vec_>
auto
subVector(Vec_&& v,
          size_t start,
          size_t stop)
    {
    static_assert(!std::is_same<Vec_&&,Vector&&>::value,"Cannot pass temp/rvalue Vector to subVector");
    auto offset = start;
    return makeRef(std::forward<Vec_>(v).store()+offset,VecRange(stop-start));
    }

} //namespace itensor

#endif
