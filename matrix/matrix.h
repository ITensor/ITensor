// matrix.h -- S.R. White 9/94

#ifndef _matrix_h
#define _matrix_h

#include "matrixref.h"
#include <vector>


namespace itensor {

class Vector;			// Defined later 
class SparseMatrix;

// The Matrix class
class Matrix : public MatrixRef
    {
public:
// Making and resizing:
    inline Matrix (int, int);		// nrows, ncols
    inline void ReDimension(int, int);	// Change size
    void ReduceDimension(int, int);	// ReDimension without reducing storage 
    void Enlarge(int, int);		// Change size while keeping contents

// Making out of other things:
    inline Matrix(const MatrixRef &);
    inline void CopyDestroy(Matrix &);	// Assignment by destroying other matrix 
    void CopyPointer(const Matrix &);	// Reference-count copy of pointer
    inline MatrixRef& operator = (const MatrixRef &);
    inline MatrixRef& operator = (const MatrixMatrixRes &);
    inline MatrixRef& operator = (Real);
    MatrixRef& operator=(const SparseMatrix &);

// Element access: (pretty efficient)
    inline Real& operator() (int, int);	// Start at 1
    inline Real  operator() (int, int) const;
    inline Real& el(int, int);		// Start at 0 
    inline Real  el(int, int) const;

    inline Matrix operator + (Real) const;
    inline friend Matrix operator + (Real, const Matrix&);
    inline Vector vector() const;	// First column

    inline MatrixRef &operator = (const Matrix &);
    inline Matrix (const Matrix &);	// Copy constructor 
    inline Matrix ();			// Default Constructor 
    inline ~Matrix ();
    inline int Storage() const; 	// Number of Reals allocated
    inline int memory() const;		// return memory used in bytes 
    inline void MakeTemp();

    void read(std::istream& s);

    inline static int GetNumMats();	// return the total number of mats 
    inline static void ResetNumMats();	// Reset matrix counter to zero 
    inline static int GetNumCon();	// return the total number of con 

    ARRAY1H_DEFS(Matrix)
private:
    char temporary;             // 1 if current matrix is a temporary
protected:
    void makematrix(int, int);	// Real Resize/Constructor 
    void copy(const Matrix &);	// real copy 
    void copytransfer(Matrix &);// copy by grabbing storage 
    inline void init();			// Initialize null matrix 
    inline void fixref();		

    static int& nummats()
        {
        static int nummats_ = 0; // number of news - number of deletes
        return nummats_;
        }

    static int& numcon()
        {
        static int numcon_ = 0;	// number of constructor calls 
        return numcon_;
        }
    };

// Functions not a member of Matrix class

Matrix Inverse(const MatrixRef &);
void Inverse(const MatrixRef &M,Matrix& result,Real& logdet, Real& sign);
void Inversenew(const MatrixRef &M,Matrix& result,Real& logdet, Real& sign);

Matrix Exp(const MatrixRef &);

Matrix Solve(const MatrixRef &A, const MatrixRef &B);	// return A^-1 * B 

void Orthog(const MatrixRef &, int nr = 0, int numpass = 2);
void Orthog(const MatrixRef& Mre, const MatrixRef& Mim, int nr = 0, int numpass = 2);

void 
QRDecomp(const MatrixRef& M, Matrix& Q, Matrix& R);

void
SchurDecomp(const MatrixRef& M, Matrix& T, Matrix& Z);

// one argument means do all columns < rows 

void EigenValues(const MatrixRef &, Vector &, Matrix &);
void GenEigenValues(const MatrixRef&, Vector&, Vector&);
void GenEigenValues(const MatrixRef& A, Vector& Re, Vector& Im, Matrix& ReV, Matrix& ImV);
void HermitianEigenvalues(const Matrix& re, const Matrix& im, 
                          Vector& evals,
	                      Matrix& revecs, Matrix& ievecs);
void
ComplexEigenvalues(const MatrixRef& Mre, const MatrixRef& Mim,
                   Vector& revals, Vector& ievals,
                   Matrix& revecs, Matrix& ievecs);
void
GeneralizedEV(const MatrixRef& A, const MatrixRef& B, Vector& D, Matrix& Z);

//void SVD(const MatrixRef& A, Matrix& U, Vector& d, Matrix& V);
//void newSVD(const MatrixRef& A, Matrix& U, Vector& d, Matrix& V);
//void thinsvd(const Matrix& A,Matrix& U, Vector& d,  Matrix &V,Real cutsq = 1.0e-20);
void checkSVD(const Matrix& A,const Matrix &U, const Vector &d, const Matrix& V);
void newcomplexSVD(const Matrix& Are, const Matrix& Aim, Matrix& Ure, Matrix& Uim, 
	Vector& d, Matrix& Vre, Matrix& Vim);
Real check_complex_SVD(const Matrix& Mre, const Matrix& Mim, const Matrix& Ure, 
	        const Matrix& Uim, const Vector& d, const Matrix& Vre, const Matrix& Vim);

void Determinant(const MatrixRef &M,Real& logdet,Real& sign);

Real Determinant(const MatrixRef &M);

void resetev(MatrixRef &);

#ifdef __alpha
int FFT(const VectorRef& in, Vector& outre, Vector& outim);
#endif

std::ostream & operator << (std::ostream &, const Matrix &);	// Overload for I/O 

class Vector : public VectorRef
    {
public:
    inline VectorRef& operator = (const VectorRef &);
    inline VectorRef& operator = (const MatrixVectorRes &);
    inline VectorRef& operator = (const VectorVectorRes &);
    inline VectorRef& operator = (Real);
    inline VectorRef& operator = (const Vector &);

// Making and resizing:
    explicit Vector (int);
    Vector (int, Real);
    explicit Vector (const std::vector<Real>& v);
    void ReDimension(int);
    void ReduceDimension(int);
    void Enlarge(int);			// Change size while keeping contents
// Element access 
    inline Real &operator() (int);
    inline Real operator() (int) const;
    inline Real &operator[] (int);
    inline Real operator[] (int) const;
    inline Real &el(int);		// Same as []
    inline Real el(int) const;

// Making out of other things:
    inline Vector (const VectorRef &);
    inline Vector (const Vector &);
    void CopyPointer(const Vector &);	// Reference-count copy of pointer
    inline void CopyDestroy(Vector &);
    inline void MakeTemp();

    inline int Storage() const;
    inline int memory() const;		// return memory used in bytes 
    inline Vector ();
    inline ~Vector ();

    void read(std::istream& s);

    friend class SparseVector;

    ARRAY1H_DEFS(Vector)
private:
    char temporary;             // 1 if current vector is a temporary
protected:
    void makevector(int);	// Real Resize/Constructor 
    void copy(const Vector &);	// real copy 
    void copytransfer(Vector &);// copy by grabbing storage 
    inline void init();
    inline void fixref();

    static int& numvecs()
        {
        static int numvecs_ = 0; // number of news - number of deletes
        return numvecs_;
        }

    static int& numcon()
        {
        static int numcon_ = 0;	// number of constructor calls 
        return numcon_;
        }
    };

// Functions not members of Vector class

void Sort(Vector &);
void Sort(Vector &, IntArray1 &);

std::ostream& operator<<(std::ostream &, const Vector &);

// Inlines for matrix.h

inline void Matrix::init()
    { MatrixRef::init(); 
      temporary = 0; 
      //numcon++; 
    }

inline void Matrix::fixref()		
    { transpose=0; scale=1.0; rowstride=ncols; store=slink.Store(); }

inline Matrix::Matrix ()
    { init(); }

inline Matrix::Matrix(const Matrix &M)	
    { init(); copy(M); }

inline Matrix::Matrix(int s1, int s2)
    { init(); makematrix(s1, s2); }

inline Matrix::Matrix(const MatrixRef &M)
    { init(); makematrix(M.Nrows(), M.Ncols()); MatrixRef::operator=(M); }

inline Matrix::~Matrix ()
    { 
        makematrix(0, 0); 
        Matrix::numcon()--; 
    }

inline void Matrix::ReDimension(int s1, int s2)
    { makematrix(s1, s2); }

inline int Matrix::Storage() const
    { return slink.Storage(); }

inline MatrixRef & Matrix::operator = (const Matrix &M)
    {
    if (&M != this) copy(M);
    return *this;
    }

// Assignment by destroying other matrix:
inline void Matrix::CopyDestroy(Matrix &M)
    { if (&M != this) copytransfer(M); }

inline MatrixRef & Matrix::operator = (const MatrixRef &M)
    {
    makematrix(M.Nrows(), M.Ncols());
    MatrixRef::operator=(M);
    return *this;
    }

inline MatrixRef& Matrix::operator = (const MatrixMatrixRes &R)
    {
    makematrix(R.Nrows(), R.Ncols());
    MatrixRef::operator=(R);
    return *this;
    }

inline MatrixRef& Matrix::operator = (Real a)
    { MatrixRef::operator=(a); return *this; }

inline Matrix Matrix::operator + (Real a) const
    { Matrix res(*this); res += a; return res; }

inline Matrix operator + (Real a, const Matrix& M)
    { return M+a; }

inline int Matrix::GetNumMats()
    { return Matrix::nummats(); }

inline void Matrix::ResetNumMats()
    { Matrix::nummats() = 0; }

inline int Matrix::GetNumCon()
    { return Matrix::numcon(); }

inline void Matrix::MakeTemp()
    { temporary = 1; }

inline int Matrix::memory() const
    { return sizeof(MatrixRef) + slink.memory(); }

inline Real & Matrix::el(int i1, int i2)
    {
    CHECKIND0(i1,i2);
    return store[i1 * ncols + i2];
    }

inline Real Matrix::el(int i1, int i2) const
    {
    CHECKIND0(i1,i2);
    return store[i1 * ncols + i2];
    }

inline Real 
Matrix::operator() (int i1, int i2) const
    {
    CHECKIND(i1,i2);
    return store[(i1 - 1) * ncols + i2 - 1];
    }

inline Real &
Matrix::operator() (int i1, int i2)
    {
    CHECKIND(i1,i2);
    return store[(i1 - 1) * ncols + i2 - 1];
    }

inline void Vector::init()
    { 
        VectorRef::init(); 
        temporary = 0; 
        Vector::numcon()++; 
    }

inline void Vector::fixref()		
    { scale = 1.0; stride = 1; store=slink.Store(); }

inline Vector::Vector () 
    { init(); }

inline Vector::Vector (const Vector &V)
    { init(); copy(V); }

inline Vector::Vector (int s)
    { init(); makevector(s); }

inline Vector::Vector (int s, Real val)
    { init(); makevector(s); operator=(val); }

inline Vector::Vector (const std::vector<Real>& v)
    { 
    init(); 
    makevector(v.size()); 
    for(size_t j = 0; j < v.size(); ++j)
        {
        operator[](j) = v[j];
        }
    }

inline Vector::Vector (const VectorRef &V)
    { init(); makevector(V.Length()); VectorRef::operator=(V); }

inline Vector::~Vector ()
    { 
        makevector(0); 
        Vector::numcon()--; 
    }

inline int Vector::Storage() const
    { return slink.Storage(); }

inline void Vector::ReDimension(int s)
    { makevector(s); }

inline VectorRef& Vector::operator = (const Vector &V)
    { if (&V != this) copy(V); return *this; }

inline void Vector::CopyDestroy(Vector &V)
    { if (&V != this) copytransfer(V); }

inline VectorRef& Vector::operator = (const VectorRef &V)
    { makevector(V.Length()); VectorRef::operator=(V); return *this; }

inline VectorRef& Vector::operator = (const MatrixVectorRes &R)
    { makevector(R.Length()); VectorRef::operator=(R); return *this; }

inline VectorRef& Vector::operator = (const VectorVectorRes &R)
    { makevector(R.Length()); VectorRef::operator=(R); return *this; }

inline VectorRef& Vector::operator = (Real a)
    { VectorRef::operator=(a); return *this; }

inline void Vector::MakeTemp()
    { temporary = 1; }

inline int Vector::memory() const		// return memory used in bytes 
    { return sizeof(VectorRef) + slink.memory(); }

inline Real & Vector::operator() (int i)
    { CHECKINDEX(i); return store[i - 1]; }

inline Real Vector::operator() (int i) const
    {
    CHECKINDEX(i);
    return store[i - 1];
    }

inline Real & Vector::operator[] (int i)
    { CHECKINDEX0(i); return store[i]; }

inline Real Vector::operator[] (int i) const
    { CHECKINDEX0(i); return store[i]; }

inline Real & Vector::el(int i)
    { return (*this)[i]; }

inline Real Vector::el(int i) const
    { return (*this)[i]; }

inline Vector Matrix::vector() const	// Make a Vector from a matrix 
    {
    if (Ncols() != 1)
	_merror("Bad conversion from Matrix to Vector");
    return Column(1);
    }

inline MatrixMatrixRes::operator Matrix() const
    { Matrix res(Nrows(),Ncols()); res = *this; return res;  }

inline MatrixVectorRes::operator Vector() const
    { Vector res(Length()); res = *this; return res;  }
// inline Vector& vector(const MatrixVectorRes& MV) 
//     { Vector res(MV.Length()); res = MV; return res;  }

inline VectorVectorRes::operator Vector() const
    { Vector res(Length()); res = *this; return res;  }

};

#include "svd.h"


namespace itensor {
ARRAY1CC_DEFS(Matrix)
ARRAY1CC_DEFS(Vector)
};

#endif
