//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef __ITENSOR_MATRIX_VEC_H_
#define __ITENSOR_MATRIX_VEC_H_

#include "itensor/tensor/vecrange.h"
#include "itensor/tensor/ten.h"

namespace itensor {

template<typename V>
using VecRefc = TenRefc<VecRange,V>;
template<typename V>
using VecRef = TenRef<VecRange,V>;
template<typename V>
using Vec = Ten<VecRange,V>;

using Vector     = Vec<Real>;
using VectorRef  = VecRef<Real>;
using VectorRefc = VecRefc<Real>;

using CVector     = Vec<Cplx>;
using CVectorRef  = VecRef<Cplx>;
using CVectorRefc = VecRefc<Cplx>;

using Vector1     = Ten<VecRange1,Real>;
using VectorRef1  = TenRef<VecRange1,Real>;
using VectorRefc1 = TenRefc<VecRange1,Real>;

using CVector1     = Ten<VecRange1,Real>;
using CVectorRef1  = TenRef<VecRange1,Real>;
using CVectorRefc1 = TenRefc<VecRange1,Real>;

template<typename V>
using hasVecRange = std::is_base_of<VecRangeType,typename stdx::decay_t<V>::range_type>;

template<typename V>
auto
stride(VecRefc<V> const& v) -> decltype(v.stride(0)) { return v.stride(0); }

template<typename V>
auto
stride(Vec<V> const& v) -> decltype(v.stride(0)) { return v.stride(0); }

VectorRef
operator*=(VectorRef v, Real fac);

CVectorRef
operator*=(CVectorRef v, Real fac);

CVectorRef
operator*=(CVectorRef v, Cplx fac);

VectorRef
operator/=(VectorRef v, Real fac);

CVectorRef
operator/=(CVectorRef v, Real fac);

CVectorRef
operator/=(CVectorRef v, Cplx fac);

VectorRef
operator+=(VectorRef a, VectorRefc b);

VectorRef
operator-=(VectorRef a, VectorRefc b);

//Copy data referenced by b to data referenced by a
void
operator&=(VectorRef a, VectorRefc const& b);

void
operator&=(CVectorRef a, CVectorRefc const& b);
void
operator&=(CVectorRef a, VectorRefc const& b);

//Dot product
Real
operator*(VectorRefc a, VectorRefc b);
//Complex dot product, conjugates a: \sum_j conj(a_j)*b_j
Cplx
operator*(CVectorRefc a, CVectorRefc b);

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

CVector inline
operator*(CVector A, Cplx fac) { A *= fac; return A; }

CVector inline
operator*(Cplx fac, CVector A) { A *= fac; return A; }

CVector inline
operator/(CVector A, Cplx fac) { A /= fac; return A; }

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
operator+(Vector&& A, Vector&& B) 
    { 
    Vector a(std::move(A));
    a += B;
    return a;
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

Vector inline
operator-(Vector&& A, Vector&& B) 
    { 
    Vector a(std::move(A));
    a -= B;
    return a;
    }

template<>
Real
norm(VectorRefc const& v);

Vector
randomVec(long size);

CVector
randomCVec(long size);

CVector inline
randomVecC(long size) { return randomCVec(size); }

Real
sumels(VectorRefc v);

void 
resize(Vector & v, size_t newsize);

void 
resize(CVector & v, size_t newsize);

template<typename T>
void 
resize(VecRefc<T> const& v, size_t newsize)
    {
    if(v.size() != newsize)
        {
        auto msg = format("Vector ref has wrong size, expected=%d, actual=%d",newsize,v.size());
        throw std::runtime_error(msg);
        }
    }


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

template<typename T>
auto
makeVecRef(T * p,
           size_t size)
    -> VecRef<T>
    {
    return VecRef<T>({p,size},VecRange(size));
    }

template<typename T>
auto
makeVecRef(T const* p,
           size_t size)
    -> VecRefc<T>
    {
    return VecRefc<T>({p,size},VecRange(size));
    }

template<typename T>
auto
makeVecRefc(T const* p,
            size_t size)
    -> VecRefc<T>
    {
    return makeVecRef(p,size);
    }

template<typename T>
auto
makeVecRef(T* p,
           size_t size,
           size_t stride)
    -> VecRef<T>
    {
    return VecRef<T>({p,size*stride},VecRange(size,stride));
    }

template<typename T>
auto
makeVecRef(T const* p,
           size_t size,
           size_t stride)
    -> VecRefc<T>
    {
    return VecRefc<T>({p,size*stride},VecRange(size,stride));
    }

template<typename T>
auto
makeVecRefc(T const* p,
            size_t size,
            size_t stride)
    -> VecRefc<T>
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
    -> decltype(makeRef(std::forward<Vec_>(v).store(),VecRange{}))
    {
    static_assert(!std::is_same<Vec_&&,Vector&&>::value,"Cannot pass temp/rvalue Vector to subVector");
    auto offset = start;
    return makeRef(std::forward<Vec_>(v).store()+offset,VecRange(stop-start));
    }

} //namespace itensor

#endif
