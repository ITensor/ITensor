#include <math.h>
#include "itensor/matrix/sparseref.h"
using std::ostream;
using std::istream;

namespace itensor {

SparseRef::SparseRef(const MatrixRef& M)
    : base(new SparseMatBase(M.Nrows(),M.Ncols())), scale(1.0), transpose(0)
    {
    int nrows = M.Nrows(), ncols = M.Ncols(), i,j;
    IntArray1 where(ncols);
    for (i = 1; i <= nrows; i++)
	{
	int count = 0;
	for (j = 1; j <= ncols; j++)
	    {			// Count up non-zero elements in row
	    Real x = M(i,j);
	    if (fabs(x) >= SparseVector::thresh())
		where[++count] = j;
	    }
	SparseVector& row(Row(i));
	row.ReDimension(count,count);
	row.sorted = 1;
	for (j = 1; j <= count; j++)
	    {
	    row.data(j) = M(i,where(j));
	    row.index[j] = where(j);
	    }
	}
    }

SparseRef& SparseRef::operator+=(const MatrixRef& M)
    {
    Matrix temp(M.Nrows(),M.Ncols());
    PutInMatrix(temp);
    temp += M;
    SparseRef s(temp);
    *this = s;
    return *this;
    }

void SparseRef::PutInMatrix(Matrix& mat) const
    {
    mat.ReduceDimension(Nrows(),Ncols());
    if(Nrows()*Ncols() == 0) return;
    mat = 0.0;
    int i,j;
    for (i = 1; i <= base->nrows; i++)
	{
	SparseVector& row(Row(i));
	int len = row.Length();
	if(!transpose)
	    for (j = 1; j <= len; j++)
		mat(i,row.index(j)) = row.data(j) * scale;
	else
	    for (j = 1; j <= len; j++)
		mat(row.index(j),i) = row.data(j) * scale;
	}
    }

void mult(const SparseRef &S, const MatrixRef &A,
		  MatrixRef &B,ClearFlag cf)
    {
    if (B.Nrows() != S.Nrows() || A.Nrows() != S.Ncols() || B.Ncols() != A.Ncols())
	error("mult(S,A,B): size mismatch");
    if(cf == Clear)
	B = 0.0;

    VectorRefNoLink Br,Ar;
    if(!S.transpose)
	for (int i = 1; i <= S.Nrows(); i++)
	    {
	    Br << B.Row(i);
	    SparseVector& row(S.Row(i));
	    for(int j = 1; j <= row.Length(); j++)
		Br += row.data(j) * S.scale * A.Row(row.index(j));
	    }
    else
	for (int j = 1; j <= S.Ncols(); j++)
	    {
	    Ar << A.Row(j);
	    SparseVector& row(S.Row(j));
	    for(int i = 1; i <= row.Length(); i++)
		B.Row(row.index(i)) += row.data(i) * S.scale * Ar;
	    }
    }

void mult(const MatrixRef &A, const SparseRef &S,
		  MatrixRef &B,ClearFlag cf)
    {
    if (B.Nrows() != A.Nrows() || S.Nrows() != A.Ncols() 
	    || B.Ncols() != S.Ncols())
	error("mult(A,S,B): size mismatch");
    if(cf == Clear)
	B = 0.0;

    int i,j;
    if(!S.transpose)
	for (i = 1; i <= S.Nrows(); i++)
	    {
	    VectorRef Ac(A.Column(i));
	    SparseVector& row(S.Row(i));
	    for(j = 1; j <= row.Length(); j++)
		B.Column(row.index(j)) += row.data(j) * S.scale * Ac;
	    }
    else
	for (j = 1; j <= S.Ncols(); j++)	// Really columns
	    {
	    VectorRef Bc(B.Column(j));
	    SparseVector& row(S.Row(j));
	    for(i = 1; i <= row.Length(); i++)
		Bc += row.data(i) * S.scale * A.Column(row.index(i));
	    }
    }

void SparseRef::write(ostream& s) const
    {
    int nrows = Nrows(), ncols = Ncols();
    s.write((char *) &(nrows), sizeof(nrows));
    s.write((char *) &(ncols), sizeof(ncols));
    s.write((char *) &(scale), sizeof(scale));
    s.write((char *) &(transpose), sizeof(transpose));
    for (int i = 1; i <= base->nrows; i++)
	{
	SparseVector& row(Row(i));
	row.write(s);
	}
    }

void SparseRef::read(istream& s) 
    {
    int nrows, ncols;
    s.read((char *) &(nrows), sizeof(nrows));
    s.read((char *) &(ncols), sizeof(ncols));
    SparseRef S(nrows,ncols);
    s.read((char *) &(S.scale), sizeof(S.scale));
    s.read((char *) &(S.transpose), sizeof(S.transpose));

    for (int i = 1; i <= nrows; i++)
	{
	SparseVector& row(S.Row(i));
	row.read(s);
	}
    *this = S;
    }

void SparseRef::ClearZeroes()
    {
    int i,j;
    IntArray1 ind(1000);
    Vector dat(1000);
    for (i = 1; i <= base->nrows; i++)
	{
	SparseVector& row(Row(i));
	int len = row.Length(), numzero = 0;
	for (j = 1; j <= len; j++)
	    if(row.data(j) == 0.0)
		numzero++;
	if(numzero == 0) continue;
	ind->ReduceDimension(len-numzero);
	dat.ReduceDimension(len-numzero);
	int k = 1;
	for (j = 1; j <= len; j++)
	    if(row.data(j) != 0.0)
		dat(k) = row.data(j), ind[k] = row.index(j), k++;
	len -= numzero;
	row.index->ReDimension(len);
	row.data.ReDimension(len);
	for (j = 1; j <= len; j++)
	    row.index[j] = ind(j), row.data(j) = dat(j);
	}
    }

int SparseRef::memory() const
    {
    int res = 0, i;
    for (i = 1; i <= base->nrows; i++)
	res += Row(i).memory();
    res /= NumRef();
    res += sizeof(base) + sizeof(transpose) + sizeof(scale);
    return res;
    }

};
