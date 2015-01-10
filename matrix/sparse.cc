// sparse.cc -- Routines for sparse matrix class

#include "sparse.h"
#include <iomanip>
#include <math.h>
//#include <memory.h>
#include "indent.h"
#ifdef _CRAY
#define VECTOR #pragma _CRI ivdep
#else
#define VECTOR
#endif

namespace itensor {

using std::istream;
using std::ostream;
using std::endl;

void SparseVector::AddElement(int i, Real value) 
// SparseMatrix must make sure the storage is big enough.
    {
#ifdef	MATRIXBOUNDS
    if (i < 1 || data.length >= data.Storage())
	_merror("SparseVector::AddElement: Bad index range");
#endif

    data.length++;
    data(Length()) = value;
    index[Length()] = i;
    if(Length() > 1)sorted = 0;
    else sorted = 1;
    }

// Redimension the storage of a SparseVector appropriately

int 
SparseVector::adjustsize(int minsize,int maxsize)
    {
    int cursize = Storage();

    if (cursize >= minsize && cursize <= maxsize)
	return cursize;
    if (maxsize < minsize)
	maxsize = minsize;

    int s = std::max(minsize, int (cursize * SparseVector::efactor()));
    s = std::max(SparseVector::minrsize(), s);
    s = std::min(maxsize, s);

    if (Length() > 1 && !sorted)
	{				// sort and consolidate the row 
	sort();
	if (Length() >= minsize && Length() < int (cursize * 0.80))
	    return Length();
	}

    SparseVector newvec(s);
    if(Length() > 0)
	{
	newvec.ReduceDimension(Length());
	newvec.data = data;
	for (int i = 1; i <= Length(); i++)
	    newvec.index[i] = index(i);
	}
    CopyDestroy(newvec);
    return Length();
    }

void 
SparseVector::sort()
    {
    sorted = 1;
    void hpsortir(int, int *, Real *);	// Heap sort int and Real arrays 

    if(Length() < 2) return;

    hpsortir(Length(), (int *) index.Store(), data.Store());	
					// Sort the arrays in place 

    int kk = 1;			// Consolidate elements that are the same 
    for (int i = 2; i <= Length(); i++)
	if (index(i) == index(i-1))
	    data(kk) += data(i);
	else if(++kk != i)
	    { data(kk) = data(i); index[kk] = index(i); }
    data.ReduceDimension(kk);		// Adjust length.
    }

int 
binarysearch(const IntArray1& index,int lower, int upper, int i)
	// Assume I is sorted, low to high; search for matching index.
    {
    if(upper < lower) return -1;
    if(index(lower) == i) return lower;
    if(index(upper) == i) return upper;
    int mid,mi;
    while(upper-lower > 1)
	{
	mid = (lower + upper) / 2;
	if((mi=index(mid)) == i) return mid;
	else if(mi < i) lower = mid;
	else  upper = mid;
	}
    return -1;
    }

// Get element(i) of sparse vector

Real 
SparseVector::el(int i) const
    {
    i++;
#ifdef	MATRIXBOUNDS
    if (i < 1)
	_merror("SparseVector::el: Bad index range");
#endif

    if (sorted)
	{
	int ind = binarysearch(index,1,Length(),i);
	if(ind == -1) return 0.0;
	return data(ind);
	}
    else
	{
	Real val = 0;
	for (int k = 1; k <= Length(); k++)
	    if (index(k) == i)
		val += data(k);
	return val;
	}
    }

Real 
dot(const SparseVector & S1, const SparseVector & S2)
    {
    if(!S1.sorted || !S2.sorted)
	_merror("S1 or S2 not sorted in dot");
     int l1 = S1.Length(), l2 = S2.Length();
    if(l1 < 1 || l2 < 1) return 0.0;
    Real res = 0.0;
     int si1,si2,i1=1,i2=1;
    while(1)
	{
	if((si1=S1.index(i1)) == (si2=S2.index(i2)))
	    {
	    res += S1.data(i1++) * S2.data(i2++);
	    if (i1 > l1 || i2 > l2) break;
	    }
	else if(si1 > si2)
	    { if(++i2 > l2) break; }
	else
	    { if(++i1 > l1) break; }
	}
    return res;
    }

Real 
dot(const SparseVector & S, const VectorRef & V)
    {
/*
    Real* vp = V.Store();
VECTOR
    int i, s = S.Length();
    Real res = 0.0;
    for(i=1; i <= s; i++)
	res += S.data(i) * vp[S.index(i)];
    return scale*res;

*/
    VectorRef Vp(V);
    Real scale = Vp.Scale();
    Vp /= scale;
    VectorRefBare VB(Vp);
    Real res = 0.0;
VECTOR
    for(int i=1; i <= S.Length(); i++)
	res += S.data(i) * VB(S.index(i));
    return scale*res;
    }

void 
sparsesaxpy(Vector &V, Real alpha, const SparseVector & S)
// V += alpha * S
    {
VECTOR
    for(int i=1; i <= S.Length(); i++)
	V(S.index(i)) += S.data(i) * alpha;
    }

void SparseVector::read(istream& s)
    {
    int rsize;
    s.read((char *) &rsize, sizeof(rsize));
    ReDimension(rsize,rsize);
    if (rsize != 0)
	{
	s.read((char *) data.Store(), rsize * sizeof(Real));
	s.read((char *) index.Store(), rsize * sizeof(int));
	}
    sorted = 1;
    }

// Write a SparseVector out in binary.  Do not delete storage.
void SparseVector::write(ostream& s)
    {
    checksort();

    int rsize = Length();
    s.write((char *) &rsize, sizeof(rsize));
    if (rsize != 0)
	{
	s.write((char *) data.Store(),  rsize * sizeof(Real));
	s.write((char *) index.Store(), rsize * sizeof(int));
	}
    }

// SparseMatrix:

// Constructors/Destructors

void SparseMatrix::init()    // Real Default constructor -- sets size to 0.
    {
    nrows = 0;
    ncols = 0;
    diag.ReDimension(0);
    noff = 0;
    row = 0;
    sorted = 1;
    temporary = 0;
    }

// Real constructor/resize

void SparseMatrix::make(int r,int c) 
    {
    if (c == ncols && r == nrows)
	{				// Matrix is the right size already 
    	for (int i = 0; i < nrows; i++)
	    row[i].ReduceDimension(0);
	return;
	}

    for (int i = 0; i < nrows; i++)
	row[i].ReDimension(0);
    delete [] row;
    row = 0;

    if (r < 0 || c < 0)
	_merror("SparseMatrix::make: Bad arguments");
    diag.ReDimension(std::max(0,std::min(r,c)));
    nrows = r;
    ncols = c;
    noff = 0;
    sorted = 1;
    if (nrows > 0)
	row = new SparseVector[nrows];
    }

// Sort all the rows in the matrix.  Does not reduce storage.

void SparseMatrix::Sort()
    {
    for (int i = 0; i < nrows; i++)
	row[i].checksort();
    }

// Copy from another sparse matrix

void SparseMatrix::copy(const SparseMatrix& S)
    {
    if (S.temporary)
	{
	SparseMatrix *pointer;	// loophole in the regs 
	pointer = (SparseMatrix *) & S;
	copytransfer(*pointer);
	return;
	}
    else
	{
	make(S.nrows, S.ncols);

	noff = S.noff;

	diag = S.diag;
	for (int i = 0; i < nrows; i++)
	    row[i] = S.row[i];
	if (!S.sorted)
	    Sort();
	temporary = 0;
	}
    }

// Assign by grabbing other SparseMatrix's storages and redimming it to 0

void SparseMatrix::copytransfer(SparseMatrix& S)
    {
    nrows = S.nrows;
    ncols = S.ncols;
    diag.CopyDestroy(S.diag);
    noff = S.noff;
    row = S.row;
    sorted = S.sorted;
    temporary = 0;
    S.init();
    }

// Copy from a general matrix
                                                
void SparseMatrix::copy(const MatrixRef& M)
    {
    make(M.Nrows(), M.Ncols());	// Make a blank matrix 

    diag = M.Diagonal();

    IntArray1 where(ncols);
    int i,j;
    for (i = 1; i <= nrows; i++)
	{
	int count = 0;
	for (j = 1; j <= ncols; j++)
	    {			// Count up non-zero elements in row
	    Real x = M(i, j);
	    if (fabs(x) >= SparseVector::thresh() && i != j)
		where[++count] = j;
	    }
	row[i-1].ReDimension(count,count);
	row[i-1].sorted = 1;
	for (j = 1; j <= count; j++)
	    {
	    row[i-1].data(j) = M(i, where(j));
	    row[i-1].index[j] = where(j);
	    }
	noff += count;
	}
    sorted = 1;
    temporary = 0;
    }
// Get element(i,j) of sparse matrix

Real SparseMatrix::el(int i, int j) const
    {
#ifdef	MATRIXBOUNDS
    if (i < 0 || i >= nrows || j < 0 || j >= ncols)
	_merror("SparseMatrix::el: Bad index range");
#endif

    if (i == j) return (diag.el(i));	// element is diagonal 
    else	return row[i].el(j);
    }

// Adds element i,j to matrix.  Indices start at 0.  Always puts new element
// at end, makerow automatically sorts and collapses rows.

void SparseMatrix::AddElement(int i, int j, Real value) 
    {
    if (i < 1 || i > nrows || j < 1 || j > ncols)
	_merror("SparseMatrix::AddElement: Bad index range");

    if (i == j)
	{				// element is diagonal 
	diag(i) += value;
	return;
	}

    i--;
    if (fabs(value) < SparseVector::thresh())
	return;
    int newlen = row[i].Length() + 1;
    if(newlen > row[i].Storage())
	row[i].adjustsize(newlen,ncols);
    row[i].AddElement(j,value);
    noff++;
    sorted = 0;

    return;
    }

// Set a Sparse Matrix equal to the identity times a Real.  Does not delete
// storage.

SparseMatrix& SparseMatrix::operator=(Real value)
    {
    if (nrows == 0 || ncols == 0) return *this;

    diag = value;
    for (int i = 0; i < nrows; i++)
	row[i].ReduceDimension(0);

    noff = 0;
    sorted = 1;

    return *this;
    }

// Multiply each element by a Real

SparseMatrix& SparseMatrix::operator*=(Real value)
    {
    if (nrows == 0 || ncols == 0)
	return *this;

    diag *= value;
    int i;

    for (i = 0; i < nrows; i++)
	row[i] *= value;

    return *this;
    }

// Add two sparse matrices

SparseMatrix& SparseMatrix::operator+=(const SparseMatrix& S)
    {
    if (nrows != S.nrows || ncols != S.ncols)
	_merror("SparseMatrix::operator+=: Matrices incompatible");

    diag += S.diag;
    for (int i = 0; i < nrows; i++)
	for (int j = 0; j < S.row[i].Length(); j++)
	    AddElement0(i, S.Column0(i,j), S.row[i].data.el(j));

    return *this;
    }

// Multiply by a dense vector, returning a dense vector

Vector SparseMatrix::operator*(const VectorRef& V) const
    {
    Vector result(Nrows());
    mult(*this, V, result);
    result.MakeTemp();
    return result;
    }

//extern "C" Real SPDOT(int*,Real*,int*,Real*);  // Cray libsci library routine

void mult(const SparseMatrix& S,const VectorRef& V1,VectorRef& V2)
    {
    if (V1.Length() != S.Ncols())
	_merror("mult(S,V1,V2): Vector V1 wrong size");

    if (V2.Length() != S.Nrows())
	_merror("mult(S,V1,V2): Vector V2 wrong size");

    int ndiag = std::min(S.nrows, S.ncols);	// Multiply by diagonal elements 

    elmult(S.diag,V1.SubVector(1,ndiag),V2.SubVector(1,ndiag));
    if(ndiag < V2.Length())
	V2.SubVector(ndiag+1,V2.Length()) = 0.0;
    VectorRefBare V2bare(V2);

    for (int i = 1; i <= S.nrows; i++)
	V2bare(i) += dot(S.row[i-1],V1);
    }

Vector SparseMatrix::TransposeTimes(const VectorRef& V) const
    {
    if (V.Length() != ncols)
	_merror("SparseMatrix::TransposeTimes: Vector wrong size");

    Vector B;
    B.ReDimension(ncols);
    int ndiag = std::min(nrows, ncols);	// Multiply by diagonal elements 
    elmult(diag,V.SubVector(1,ndiag),B.SubVector(1,ndiag));
    if(ndiag < V.Length())
	B.SubVector(ndiag+1,B.Length()) = 0.0;

    int i;

    for (i = 0; i < nrows; i++)
	if(V.el(i) != 0.0)
	    sparsesaxpy(B,V.el(i),row[i]);

    return B;
    }

// do not zero out B!
void mult(const SparseMatrix &S, const MatrixRef &A, MatrixRef &B)
    {
    if (B.Nrows() != S.Nrows())
	_merror("mult(S,A,B): Matrix B wrong size");

    if (A.Nrows() != S.Ncols())
	_merror("mult(S,A,B): Matrix A wrong size");

    if (B.Ncols() != A.Ncols())
	_merror("mult(S,A,B): Matrix A wrong size");

    int ndiag = std::min(S.nrows, S.ncols);	// Multiply by diagonal elements 

    int i;
    for (i = 1; i <= ndiag; i++)
	B.Row(i) += S.diag(i) * A.Row(i);

    for (i = 1; i <= S.nrows; i++)
	{
	VectorRef Br(B.Row(i));
	int j;
	for(j = 1; j <= S.row[i-1].Length(); j++)
	    Br += S.row[i-1].data(j) * A.Row(S.row[i-1].index(j));
	}
    }

void mult(const MatrixRef &A, const SparseMatrix &S, MatrixRef &B)
    {
    if (B.Nrows() != A.Nrows())
	_merror("mult(A,S,B): Matrix B wrong size");

    if (S.Nrows() != A.Ncols())
	_merror("mult(A,S,B): Matrix A wrong size");

    if (B.Ncols() != S.Ncols())
	_merror("mult(A,S,B): Matrix A wrong size");

    int ndiag = std::min(S.nrows, S.ncols);	// Multiply by diagonal elements 

    int i;
    for (i = 1; i <= ndiag; i++)
	B.Column(i) += S.diag(i) * A.Column(i);

    for (i = 1; i <= S.nrows; i++)
	{
	VectorRef Ac(A.Column(i));
	int j;
	for(j = 1; j <= S.row[i-1].Length(); j++)
	    B.Column(S.row[i-1].index(j)) += S.row[i-1].data(j) * Ac;
	}
    }

// Return memory used in bytes

int SparseMatrix::memory() const
    {
    int mem = sizeof(SparseMatrix);	// Header 
    mem += diag.memory();

    int i;
    for (i = 0; i < nrows; i++)
	mem += row[i].memory();  // off-diagonals 

    return mem;
    }

// Print a summary of memory usage onto s

void SparseMatrix::PrintMemory(ostream& s) const
    {
    int nelements = Storage();
    int totmem = memory();
    Real efficiency = Real (totmem) / Real (nelements);

    s << "Size: " << nrows << "x" << ncols
	<< ", Nonzero elements: " << nelements
	<< ", Sparseness: " << Real (nelements) / nrows
    	<< iendl << "Bytes used: " << totmem
    	<< ", bytes/element: " << efficiency << iendl;
    }

// Read a SparseMatrix in in binary, allocating storage

void SparseMatrix::read(istream& s)
    {
    if (nrows != 0 || ncols != 0)
	_merror("SparseMatrix::read: SparseMatrix not null");

    int nr, nc;
    s.read((char *) &(nr), sizeof(nr));
    s.read((char *) &(nc), sizeof(nc));
    make(nr, nc);

    int ndiag = std::min(nrows, ncols);
    s.read((char *) diag.Store(), ndiag * sizeof(Real));

    noff = 0;
    int i;
    for (i = 0; i < nrows; i++)
	{
	row[i].read(s);
	noff += row[i].Length();
	}
    s.read((char *) &(sorted), sizeof(sorted));
    temporary = 0;
    }

// Write a SparseMatrix out in binary.  Do not delete storage.  This can be
// done from the calling routine using: S.ReDimension(0,0);

void SparseMatrix::write(ostream& s)
    {
    if (!sorted)
	Sort();			// Make sure rows are sorted first 

    s.write((char *) &(nrows), sizeof(nrows));
    s.write((char *) &(ncols), sizeof(ncols));

    int ndiag = std::min(nrows, ncols);
    s.write((char *) diag.Store(), ndiag * sizeof(Real));

    int i;
    for (i = 0; i < nrows; i++)
	row[i].write(s);

    s.write((char *) &(sorted), sizeof(sorted));
    }

// Output a sparse matrix

ostream& operator<<(ostream& s, const SparseMatrix& X)
    {
    int w = s.width();
    int nr = X.nrows;
    //long f = s.flags();
    s.setf(std::ios::fixed, std::ios::floatfield);

    s << "Sparse Matrix, dimension=" << X.nrows << "x" << X.ncols;
    s << ", diagonal Elements:" << iendl;

    int i,j;
    for (i = 0; i < std::min(X.nrows, X.ncols); i++)
	s << i << "," << i << ": " << std::setw(w) << X.diag.el(i) << iendl;

    s << iendl << "Off-diagonal Elements:" << iendl;

    for (i = 0; i < nr; i++)
	for (j = 0; j < X.Rowsize0(i); j++)
	    s << i << "," << X.Column0(i,j) << ": "
		<< std::setw(w) << X.Element0(i,j) << iendl;
    s << std::flush;
    //s.flags(f);
    return s;
    }

// Output a Sparse Matrix, assuming it is symmetric, reasonably efficiently

void SparseMatrix::PrintSymmetric(ostream& s)
    {
    //int w = s.width();
    int nr = nrows;
    //long f = s.flags();
    //s.setf(0, std::ios::floatfield);	// Make sure fixed and scientific are cleared 

    s << endl << "Sparse Matrix, dimension=" << nrows << "x" << ncols;
    s << ", diagonal Elements:" << endl;

    int i,j;
    for (i = 0; i < std::min(nrows, ncols); i++)
	s << i << " " << i << ": " << std::setw(0) << std::setprecision(7)
	    << diag.el(i) << endl;

    s << endl << "Off-diagonal Elements:" << endl;

    for (i = 0; i < nr; i++)
	for (j = 0; j < row[i].Length(); j++)
	    if (Column0(i,j) < i)
		s << i << " " << Column0(i,j) << " "
		    << std::setw(0) << std::setprecision(7) << Element0(i,j) << endl;
    s << std::flush;
    //s.flags(f);
    }

// Matrix class member -- Assign a sparse matrix to a dense one

MatrixRef &
Matrix::operator=(const SparseMatrix& S)
    {
    makematrix(S.ncols, S.nrows);

    *this = 0;			// Zero out target 

    Diagonal() = S.diag;

    int i,j;
    for (i = 0; i < nrows; i++)	// off-diagonals 
	for (j = 0; j < S.row[i].Length(); j++)
	    el(i, S.Column0(i,j)) += S.Element0(i,j);
    return *this;
    }

}; //namespace itensor

