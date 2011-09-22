// matrix.h -- S.R. White 9/94

#ifndef _matrix_h
#define _matrix_h

#include "matrixref.h"

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

    inline static int GetNumMats();	// return the total number of mats 
    inline static void ResetNumMats();	// Reset matrix counter to zero 
    inline static int GetNumCon();	// return the total number of con 

    ARRAY1H_DEFS(Matrix)
private:
    char temporary;             // 1 if current matrix is a temporary
protected:
    //static int nummats;	// number of news - number of deletes 
    //static int numcon;	// number of constructor calls 
    void makematrix(int, int);	// Real Resize/Constructor 
    void copy(const Matrix &);	// real copy 
    void copytransfer(Matrix &);// copy by grabbing storage 
    inline void init();			// Initialize null matrix 
    inline void fixref();		
    };

// Functions not a member of Matrix class

Matrix Inverse(const MatrixRef &);
void Inverse(const MatrixRef &M,Matrix& result,Real& logdet, Real& sign);
void Inversenew(const MatrixRef &M,Matrix& result,Real& logdet, Real& sign);

Matrix Exp(const MatrixRef &);

Matrix Solve(const MatrixRef &A, const MatrixRef &B);	// return A^-1 * B 

void Orthog(const MatrixRef &, int nr = 0, int numpass = 2);

// one argument means do all columns < rows 

void EigenValues(const MatrixRef &, Vector &, Matrix &);
void GenEigenValues(const MatrixRef&, Vector&, Vector&);

void SVD(const MatrixRef& A, Matrix& U, Vector& d, Matrix& V);
void newSVD(const MatrixRef& A, Matrix& U, Vector& d, Matrix& V);
void thinsvd(const Matrix& A,Matrix& U, Vector& d,  Matrix &V,Real cutsq = 1.0e-20);
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
    inline Vector (int);
    inline void ReDimension(int);
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
    //static int numvecs;		// number of news - number of deletes 
    //static int numcon;		// number of constructor calls 
    };

// Functions not members of Vector class

void Sort(Vector &);
void Sort(Vector &, IntArray1 &);

std::ostream & operator << (std::ostream &, const Vector &);

#include "matrix.ih"

#ifdef THIS_IS_MAIN

//int Matrix::nummats = 0;        // number of new's - number of delete's
//int Matrix::numcon = 0;         // Constructor counter 
//int Vector::numvecs = 0;        // number of new's - number of delete's
//int Vector::numcon = 0;         // Constructor counter 
ARRAY1CC_DEFS(Matrix)
ARRAY1CC_DEFS(Vector)

#endif

#endif
