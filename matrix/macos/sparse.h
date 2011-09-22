// sparse.h -- A sparse Matrix package -- Reinhard M. Noack 5/3/93

#ifndef _SPARSE_H
#define _SPARSE_H

#include "bigmatrix.h"		/* includes matrix.h */

const Real default_thresh = 1e-10; // Default threshold
const Real default_efactor = 1.5;  // default expansion factor for rows
const int default_minrsize = 20;   // Default minimum row size

class SparseVector
    {
public:
    IntArray1 index;
    Vector data;
    int sorted;
    static Real efactor;		// Factor for expanding rows 
    static Real thresh;		// threshold for accepting elements 
    static int minrsize;	// static default min. row size 

    SparseVector() {sorted = 1;}

    SparseVector(int storage, int size = 0) : index(storage), data(storage)
	{ data.ReduceDimension(size); sorted = 1;}

    SparseVector(const SparseVector& other) 
	: index(other.index), data(other.data) 
	{ sorted = other.sorted; }

    SparseVector & operator = (const SparseVector& other)
	{
	data = other.data; index = other.index; sorted = other.sorted;
	return *this;
	}

    int Length() const { return data.Length(); }

    int Storage() { return data.Storage(); }

    void ReDimension(int storage, int size = 0)
	{
	int s = max(storage,size);
	data.ReDimension(s); index->ReDimension(s);
	data.ReduceDimension(size);
	}

    void ReduceDimension(int s)
	{
	if(s <= Storage()) data.ReduceDimension(s);
	else			ReDimension(s,s);
	}

    int adjustsize(int minsize,int maxsize);	// Returns new size

    void CopyDestroy(SparseVector &other)
	{
	index = other.index;
	data.CopyDestroy(other.data);
	sorted = other.sorted; other.sorted = 1;
	}

    void sort();

    void checksort() 
	{ if(Length() > 1 && !sorted) sort(); sorted = 1; }

    Real el(int i) const;

    void AddElement(int i, Real v);	// Add an element (+=) 

    void AddElement0(int i, Real v)
	{ AddElement(i+1,v); }

    SparseVector & operator *= (Real fac)
    	{ data *= fac; return *this; }

    friend Real dot(const SparseVector & S1, const SparseVector & S2);

    friend Real dot(const SparseVector & S,const VectorRef & V);

    friend void sparsesaxpy(Vector &V, Real alpha, const SparseVector & S);
	// V += alpha * S

    void read(std::istream & s);	// Read and write in binary to a stream 

    void write(std::ostream & s);	// read dims new storage, write does not 
    // delete old storage 

    void SetMinRowSize(int s)	// Set it 
	{ minrsize = s; }

    Real Threshold() const
	{ return thresh; }

    void SetThresh(Real x)	// Set threshold 
	{ thresh = x; }

    void SetExFactor(Real f)	// Set the expansion factor for a row 
	{ efactor = f; }

    Real ExFactor() const	// get the expansion factor 
	{ return efactor; }

    int MinRowSize() const	// Get the minimum row size 
	{ return minrsize; }

    int memory() const		// Return amount of memory used in bytes 
    	{ return data.memory() +
	sizeof(int)*index->Size()+sizeof(IntArray1) + sizeof(int); }
    ARRAY1H_DEFS(SparseVector)
    };

inline Real 
dot(const VectorRef & V, const SparseVector & S)
    { return dot(S,V); }

class SparseMatrix : public BigMatrix 
    {
protected:
    int nrows;			// Row and column dimension 
    int ncols;
    Vector diag;
    SparseVector *row;
    int noff;			// Number of off-diagonal elements 
    int minrsize;		// Minimum row size. Can be set. 
    char sorted;		// 1 if matrix is sorted 
    char temporary;		// 1 if marked temporary 

    void copy(const SparseMatrix &);	// copy functions 
    void copytransfer(SparseMatrix &);	// Grab storage 
    void copy(const MatrixRef &);
    void make(int, int);	// Real constructor/resize 
    void init();		// Real default constructor 

    // Make a row of a given minimum,maximum size 
    void makerow(int r, int minsize, int maxsize = 0);

public:
    SparseMatrix()		// Constructors 
	{ init(); }

    SparseMatrix(int r, int c)	// rows r, columns c 
	{ init(); make(r, c); }

    SparseMatrix(const SparseMatrix & S)
	{ init(); copy(S); }

    SparseMatrix(const MatrixRef &M)
	{ init(); copy(M); }

    ~SparseMatrix()		// Destructor 
	{ make(0, 0); }

    SparseMatrix & operator = (const SparseMatrix & S)	// Assignment 
	{
	if (this != &S)
	    copy(S);
	return *this;
	}

    SparseMatrix & operator = (const MatrixRef &M)
	{
	copy(M);
	return *this;
	}

    // Set matrix equal to real times identity. Sets rowsizes to zero, 
    // but does not delete storage 
    SparseMatrix & operator = (Real);

    void CopyDestroy(SparseMatrix & M)	// Assignment by destroying other matrix 
	{
	if (&M != this)
	    copytransfer(M);
	}

    void MakeTemp()
	{ temporary = 1; }	// Mark sparsematrix as temporary so 
    				// storage can be grabbed 

    SparseMatrix & operator *= (Real);	// Multiply each element by a Real 

    SparseMatrix operator *(Real val) const
	{
	SparseMatrix result = *this;
	result *= val;
	return result;
	}

    friend SparseMatrix operator *(Real val, const SparseMatrix & S)
	{ return S * val; }

    SparseMatrix & operator += (const SparseMatrix &);	// Addition 

    SparseMatrix operator + (const SparseMatrix & S) const
	{
	SparseMatrix result = *this;
	result += S;
	return result;
	}

    int Nrows() const
	{ return nrows; }		// Information functions 

    int Ncols() const
	{ return ncols; }

    int Size() const
	{ return ncols; }

    int Storage() const
	{ return noff + min(nrows, ncols); }

    VectorRef DiagRef() const	// Pointer to diagonal elements 
	{ return diag; }

    int Rowsize(int i) const	// Return rowsize of row i 
	{ return (i > 0 && i <= nrows) ? row[i-1].Length() : -1; }

    int Rowsize0(int i) const	// Return rowsize of row i, starting at 0 
	{ return (i >= 0 && i < nrows) ? row[i].Length() : -1; }

    Real DiagElement(int i) const	// return diagonal element i 
	{ return (i > 0 && i <= min(nrows, ncols)) ? diag(i) : 0.0; }

    Real DiagElement0(int i) const  // return diagonal element i, starting at 0
	{ return (i >= 0 && i < min(nrows, ncols)) ? diag.el(i) : 0; }

    Real ODElement(int r, int i) const    // return ith OD element in Row row 
	{ return row[r-1].data(i); }

    Real ODElement0(int r, int i) const	// indices starting at 0 
	{ return row[r].data.el(i); }

    int ODColumn(int r, int i) const	// return column index 
	{ return row[r-1].index(i); }

    int ODColumn0(int r, int i) const	// indices starting at 0 
	{ return row[r].index(i+1)-1; }   // Column index starts at 0

    Real Element(int r, int i) const	// return ith OD element in Row row 
	{ return ODElement(r,i); }

    Real Element0(int r, int i) const	// indices starting at 0 
	{ return ODElement0(r,i); }

    int Column(int r, int i) const	// return column index 
	{ return ODColumn(r,i); }

    int Column0(int r, int i) const	// indices starting at 0 
	{ return ODColumn0(r,i); }

    void ReDimension(int r, int c)	// Resize matrix 
	{ make(r, c); }

    Real el(int, int) const;	// Get an element,indices starting at zero 

    Real operator() (int i, int j) const	// indices starting at 1 
	{ return el(i - 1, j - 1); }

    void AddElement(int, int, Real);	// Add an element (+=) 

    void AddElement0(int i, int j, Real v)	// Index starting at 1 
	{ AddElement(i + 1, j + 1, v); }

    Vector operator *(const VectorRef &) const;	// Multiply 
    friend void mult(const SparseMatrix &, const VectorRef &, VectorRef &);
    // basic routine 

    void product(const VectorRef &A, VectorRef &B) const
	{ mult(*this,A,B); }

    friend void mult(const SparseMatrix &, const MatrixRef &A, MatrixRef &B);
    friend void mult(const MatrixRef &A, const SparseMatrix &, MatrixRef &B);

    Vector TransposeTimes(const VectorRef &) const;	// by Transpose 

    int Sorted() const		// Return sorted flag (1 if matrix is sorted) 
	{ return sorted; }

    void Sort();		// Sort elements, consolidating rows 

    int memory() const;		// Return amount of memory used in bytes 

    void PrintMemory(std::ostream &) const;	// Print a summary of memory usage 

    void read(std::istream & s);	// Read and write in binary to a stream 
    void write(std::ostream & s);	// read dims new storage, write does not 
    // delete old storage 

    // Output a sparse matrix 
    friend std::ostream & operator << (std::ostream &, const SparseMatrix &);

    void PrintSymmetric(std::ostream &);// Print in an efficient way, assuming the 
    // matrix is symmetric.  Prints lower half. 

    friend MatrixRef &
     Matrix::operator = (const SparseMatrix &);	// Assign to a Matrix 

    friend class MatrixRef;
    };

inline std::ostream & operator << (std::ostream &s, const SparseVector &a) {return s; }

#ifdef THIS_IS_MAIN
int  SparseVector::minrsize = default_minrsize;
Real SparseVector::efactor = default_efactor;
Real SparseVector::thresh = default_thresh;
ARRAY1CC_DEFS(SparseVector)
#endif

#endif
