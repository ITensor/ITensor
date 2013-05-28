// matrix.cc -- Code for Matrix class

#include "matrix.h"
#include <iostream>
#include <iomanip>
#include <memory.h>
#include "minmax.h"
#include "sparse.h"
#include "indent.h"
#ifdef _CRAY
#define VECTOR #pragma _CRI ivdep
#else
#define VECTOR
#endif

// Real resize/constructor

using namespace std;


void 
Matrix::makematrix(int s1, int s2)
    {
    if(s1 == nrows && s2 == ncols) return;
    if (s1 < 0 || s2 < 0)
	_merror("Matrix::makematrix: bad args");
    int size = s1 * s2;
    /*
    if (Store() == 0)
	{ if (size > 0) Matrix::nummats()++; }
    else
	{ if (size == 0) Matrix::nummats()--; }
    */
    slink.makestorage(size);
    nrows = s1;
    ncols = s2;
    fixref();
    }

// ReDimension, but don't reduce storage, increase if needed
void 
Matrix::ReduceDimension(int s1, int s2)
    {
    slink.increasestorage(s1*s2);
    nrows = s1;
    ncols = s2;
    fixref();
    }

// Real copy function

void 
Matrix::copy(const Matrix & M)
    {
    if (M.temporary == 1)
	{
	Matrix *pointer;	// loophole in the regs 
	pointer = (Matrix *) & M;
	copytransfer(*pointer);
	return;
	}
    else
	{
	if (temporary != 0)
	    cerr << "Matrix::copy: Warning: temporary corrupted" << endl;
	// If Dimension is already reduced, don't remake matrix 
	if (Storage() > nrows * ncols)
	    ReduceDimension(M.nrows, M.ncols);
	else
	    makematrix(M.nrows, M.ncols);
	if (M.Store() != 0)
	    memcpy((void *) Store(), (void*) M.Store(), sizeof(Real)*nrows*ncols);
	temporary = 0;
	}
    }

// Assign by grabbing other matrix's store and redimming it to 0

void 
Matrix::copytransfer(Matrix & M)
    {
    if (Store() != 0)
	makematrix(0,0);
    MatrixRef::operator<<(M);
    M.makematrix(0,0); 
    //Matrix::numcon()--; 
    M.init();
    temporary = 0;
    }

void 
Matrix::CopyPointer(const Matrix & M)
    {
    if(&M == this) return;
    if (Store() != 0)
	makematrix(0,0);
    MatrixRef::operator<<(M);
    temporary = 0;
    }

void 
Matrix::Enlarge(int nr, int nc)
    {
    int i,j, onr = Nrows(), onc = Ncols();
    int knr=min(nr,onr), knc=min(nc,onc);	// size of submatrix kept intact
						// Everything else is set to 0
    if(nr == onr && nc == onc) return;
    if(nr*nc <= Storage())			// no need to change store
	{
	MatrixRef Mref(*this);
	ReduceDimension(nr,nc);
	if(nc > onc)
	    for(i = knr; i > 1; i--)			// First row invariant
	        for(j = knc; j >= 1; j--)
		    (*this)(i,j) = Mref(i,j);
	else if(nc < onc)
	    for(i = 2; i <= knr; i++)
	        for(j = 1; j <= knc; j++)
		    (*this)(i,j) = Mref(i,j);
	if(knc < nc)
	    Columns(knc+1,nc) = 0.0;
	if(knr < nr)
	    Rows(knr+1,nr) = 0.0;
	}
    else
	{
	Matrix temp(nr,nc);
	temp = 0.0;
	if(knr >= 1 && knc >= 1)
	    temp.SubMatrix(1,knr,1,knc) = SubMatrix(1,knr,1,knc);
	copytransfer(temp);
	}
    }

ostream & operator << (ostream & s, const Matrix & M)
    {
    int w = s.width();
    //long f = s.flags();
    s.setf(std::ios::fixed, std::ios::floatfield);
    for (int i = 1; i <= M.Nrows(); i++)
	{
	for (int j = 1; j <= M.Ncols(); j++)
	    s << std::setw(w) << M(i, j) << " ";
	s << iendl;
	}
    s << iendl;
    //s.flags(f);
    return s;
    }

void 
Vector::makevector(int s)
    {
    if(s == length) return;
    /*
    if (Store() == 0)
	{ if (s > 0) Vector::numvecs()++; }
    else
	{ if (s == 0) Vector::numvecs()--; }
    */
    slink.makestorage(s);
    length = s;
    fixref();
    }

// ReDimension, but don't reduce storage, increase if needed
void 
Vector::ReduceDimension(int s)
    { slink.increasestorage(s); length = s; fixref(); }

// Real copy function
void 
Vector::copy(const Vector & V)
    {
    if (V.temporary == 1)
	{
	Vector *pointer;	// loophole in the regs 
	pointer = (Vector *) & V;
	copytransfer(*pointer);
	return;
	}
    else
	{
	if (temporary != 0)
	    cerr << "Vector::copy: Warning: temporary corrupted" << endl;
	// If Dimension is already reduced, don't remake matrix 
	if (Storage() > length)
	    ReduceDimension(V.length);
	else
	    makevector(V.length);
	if (V.Store() != 0)
	    memcpy((void *) Store(), (void*) V.Store(), sizeof(Real)*length);
	temporary = 0;
	}
    }

void 
Vector::copytransfer(Vector & V)
    {
    if (Store() != 0)
	makevector(0);
    VectorRef::operator<<(V);
    V.makevector(0); 
    //Vector::numcon()--; 
    V.init();
    temporary = 0;
    }

void 
Vector::Enlarge(int n)
    {
    int on = Length();
    int kn = min(n,on);

    if(n <= Storage())			// no need to change store
	{
        ReduceDimension(n);
	if(n > on)
	    SubVector(on+1,n) = 0.0;
	return;
	}
    Vector temp(n);
    temp = 0.0;
    if(kn >= 1) temp.SubVector(1,kn) = SubVector(1,kn);
    copytransfer(temp);
    }

void Vector::
read(std::istream& s)
    {
    int L;
    s.read((char*)&L,sizeof(L));
    ReDimension(L);
    Real val;
    for(int k = 0; k < L; ++k)
        {
        s.read((char*)&val,sizeof(val));
        el(k) = val;
        }
    }

ostream & operator << (ostream & s, const Vector & V)
    {
    int w = s.width();
    //long f = s.flags();
    s.setf(std::ios:: fixed, std::ios::floatfield);
    for (int i = 1; i <= V.Length(); i++)
	s << std::setw(w) << V(i) << " ";
    s << "\n" << iendl;
    //s.flags(f);
    return s;
    }

void Matrix::
read(std::istream& s)
    {
    int nr;
    Vector v;
    s.read((char*)&nr,sizeof(nr));
    v.read(s);
    int nc = v.Length()/nr;
    ReDimension(nr,nc);
    this->TreatAsVector() = v;
    }
