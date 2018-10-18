//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MAT__H_
#define __ITENSOR_MAT__H_

#include "itensor/tensor/matrange.h"

namespace itensor {

template<typename V>
using MatRefc = TenRefc<MatRange,V>;
template<typename V>
using MatRef = TenRef<MatRange,V>;
template<typename V>
using Mat = Ten<MatRange,V>;

using MatrixRef = MatRef<Real>;
using MatrixRefc = MatRefc<Real>;
using Matrix = Mat<Real>;

using CMatrixRef = MatRef<Cplx>;
using CMatrixRefc = MatRefc<Cplx>;
using CMatrix = Mat<Cplx>;

using MatrixRef1 = TenRefc<MatRange1,Real>;
using MatrixRefc1 = TenRef<MatRange1,Real>;
using Matrix1 = Ten<MatRange1,Real>;

using CMatrixRef1 = TenRefc<MatRange1,Cplx>;
using CMatrixRefc1 = TenRef<MatRange1,Cplx>;
using CMatrix1 = Ten<MatRange1,Cplx>;

template<typename M>
using hasMatRange = std::is_base_of<MatRangeType,typename stdx::decay_t<M>::range_type>;

template<typename Mat_>
auto
nrows(Mat_ const& M) -> MatRange::size_type { return M.range().rn; }

template<typename Mat_>
auto
ncols(Mat_ const& M) -> MatRange::size_type { return M.range().cn; }

template<typename Mat_>
auto
rowStride(Mat_ const& M) -> MatRange::size_type { return M.range().rs; }

template<typename Mat_>
auto
colStride(Mat_ const& M) -> MatRange::size_type { return M.range().cs; }

template<typename Mat_>
bool
isTransposed(Mat_ const& M) { return isTransposed(M.range()); }

void
operator*=(MatrixRef const& M, Real fac);
void
operator*=(CMatrixRef const& M, Real fac);

void
operator/=(MatrixRef const& M, Real fac);
void
operator/=(CMatrixRef const& M, Real fac);

void
operator+=(MatrixRef const& A, MatrixRefc const& B);
void
operator+=(MatrixRef const& A, Matrix && B);
void
operator+=(CMatrixRef const& A, CMatrixRefc const& B);
void
operator+=(CMatrixRef const& A, CMatrix && B);

void
operator-=(MatrixRef const& A, MatrixRefc const& B);
void
operator-=(MatrixRef const& A, Matrix && B);
void
operator-=(CMatrixRef const& A, CMatrixRefc const& B);
void
operator-=(CMatrixRef const& A, CMatrix && B);

//Copy data referenced by B to memory referenced by A
void
operator&=(MatrixRef const& A, MatrixRefc const& B);

//Copy data of B to memory referenced by A
void inline
operator&=(MatrixRef const& A, Matrix const& B) { A &= makeRefc(B); }

// C = beta*C + alpha*A*B
template<typename VA, typename VB>
void
gemm(MatRefc<VA> A, 
     MatRefc<VB> B, 
     MatRef<common_type<VA,VB>>  C,
     Real alpha,
     Real beta);

template<typename VA, typename VB>
void
mult(MatRefc<VA> A, 
     MatRefc<VB> B, 
     MatRef<common_type<VA,VB>> C);

template<typename MatA, 
         typename MatB,
         typename MatC,
class=stdx::require<hasMatRange<MatA>,
                    hasMatRange<MatB>,
                    hasMatRange<MatC>>>
void
mult(MatA const& A, 
     MatB const& B, 
     MatC      & C)
    {
    gemm(makeRef(A),makeRef(B),makeRef(C),1.,0.);
    }


// compute matrix multiply (dgemm) A*B
// add result to memory referenced by C
void
multAdd(MatrixRefc A, 
        MatrixRefc B, 
        MatrixRef  C);

template<typename V>
void
mult(MatRefc<V> M,
     VecRefc<V> x,
     VecRef<V> y,
     bool fromleft = false);

template<typename VM, typename Vx>
Vec<common_type<VM,Vx>>
mult(MatRefc<VM> M,
     VecRefc<Vx> x,
     bool fromleft = false);

//y = y+M*x
template<typename V>
void
multAdd(MatRefc<V> M,
        VecRefc<V> x,
        VecRef<V> y,
        bool fromleft = false);

//y = y-M*x
template<typename V>
void
multSub(MatRefc<V> M,
        VecRefc<V> x,
        VecRef<V> y,
        bool fromleft = false);

//void
//mult(CMatrixRefc A,
//     CMatrixRefc B,
//     CMatrixRef  C);

//Create an NrxNc matrix
//with 1's along the diagonal
Matrix
eye(size_t Nr, size_t Nc);

//Reducing number of columns does not affect
//remaining data (column major storage)
template<typename V>
void
reduceCols(Mat<V> & M, size_t new_ncols);

template<typename V>
void
resize(Mat<V> & M, size_t nrows, size_t ncols);

template<typename T>
void 
resize(MatRefc<T> const& M, size_t nr, size_t nc)
    {
    if((nrows(M)!=nr) || (ncols(M)!=nc))
        {
        auto msg = format("Matrix ref has wrong size, expected=%dx%d, actual=%dx%d",
                          nr,nc,nrows(M),ncols(M));
        throw std::runtime_error(msg);
        }
    }


template<typename V>
Mat<V> 
operator*(MatRefc<V> const& A, Real fac);

template<typename V>
Mat<V> 
operator*(Real fac, MatRefc<V> const& A);

template<typename V>
Mat<V> 
operator*(Mat<V> && A, Real fac);

template<typename V>
Mat<V> 
operator*(Real fac, Mat<V> && A);

template<typename V>
Mat<V> 
operator/(MatRefc<V> const& A, Real fac);

template<typename V>
Mat<V> 
operator/(Mat<V> && A, Real fac);

template<typename MatA,typename MatB,class>
auto
operator+(MatA && A, MatB && B) -> Mat<common_type<MatA,MatB>>;

template<typename MatA,typename MatB,class>
auto
operator-(MatA && A, MatB && B) -> Mat<common_type<MatA,MatB>>;

template<typename MatA,
         typename MatB,
         class = stdx::require<hasMatRange<MatA>,hasMatRange<MatB>> >
Mat<common_type<MatA,MatB>> 
mult(MatA const& A,
     MatB const& B);

Vector
operator*(MatrixRefc const& A,
          VectorRefc const& b);

Vector
operator*(VectorRefc const& a,
          MatrixRefc const& B);

template<typename VA, typename VB>
auto
operator*(MatRefc<VA> const& A, MatRefc<VB> const& B) -> decltype(mult(A,B))
    { return mult(A,B); }

template<typename VA, typename VB>
auto
operator*(Mat<VA> const& A, MatRefc<VB> const& B) -> decltype(mult(makeRef(A),B))
    { return mult(makeRef(A),B); }

template<typename VA, typename VB>
auto
operator*(MatRefc<VA> const& A, Mat<VB> const& B) -> decltype(mult(A,makeRef(B)))
    { return mult(A,makeRef(B)); }

template<typename VA, typename VB>
auto
operator*(Mat<VA> const& A, Mat<VB> const& B) -> decltype(mult(makeRef(A),makeRef(B)))
    { return mult(makeRef(A),makeRef(B)); }

template<typename VA, typename VB>
Mat<VA>&
operator*=(Mat<VA> & A, MatRefc<VB> const& B) { A = mult(makeRef(A),B); return A; }

template<typename... CtrArgs>
Matrix
randomMat(CtrArgs&&... args);

template<typename... CtrArgs>
CMatrix
randomMatC(CtrArgs&&... args);

template<>
std::ostream&
operator<<(std::ostream& s, MatrixRefc const& M);

template<>
std::ostream&
operator<<(std::ostream& s, CMatrixRefc const& M);

template<typename V>
std::ostream&
operator<<(std::ostream& s, TenRef<MatRange,V> const& M) { return s << makeRefc(M); }

template<typename V>
std::ostream&
operator<<(std::ostream& s, Ten<MatRange,V> const& M) { return s << makeRefc(M); }


//
// makeMatRef functions
//

template<typename T>
auto
makeMatRef(T* p,
           size_t max_offset,
           size_t nrows,
           size_t ncols)
    -> MatRef<T>
    {
    return MatRef<T>({p,max_offset},MatRange{nrows,ncols});
    }

template<typename T>
auto
makeMatRef(T const* p,
           size_t max_offset,
           size_t nrows,
           size_t ncols)
    -> MatRefc<T>
    {
    return MatRefc<T>({p,max_offset},MatRange{nrows,ncols});
    }

template<typename T>
auto
makeMatRefc(T const* p,
            size_t max_offset,
            size_t nrows,
            size_t ncols)
    -> MatRefc<T>
    {
    return MatRefc<T>({p,max_offset},MatRange{nrows,ncols});
    }

template<typename T, 
         class = stdx::enable_if_t<not std::is_const<T>::value> >
auto
makeMatRef(DataRange<T> const& D,
           size_t nrows,
           size_t ncols)
    -> MatRef<T>
    {
    return MatRef<T>(D,MatRange{nrows,ncols});
    }

template<typename T>
auto
makeMatRef(DataRange<const T> const& D,
           size_t nrows,
           size_t ncols)
    -> MatRefc<T>
    {
    return MatRefc<T>(D,MatRange{nrows,ncols});
    }

template<typename T>
auto 
makeMatRefc(DataRange<T> const& D,
            size_t nrows,
            size_t ncols)
    -> MatRefc<stdx::remove_const_t<T>>
    {
    return MatRefc<stdx::remove_const_t<T>>(DataRange<const T>(D),MatRange{nrows,ncols});
    }

} //namespace itensor

#include "mat_impl.h"

#endif
