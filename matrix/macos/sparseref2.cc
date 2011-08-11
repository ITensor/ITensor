#include "sparseref2.h"

typedef CRowIter CI;

void SparseVec::write(ostream& s) const
    {
    int len = size();
    s.write((char *) &(len), sizeof(len));
    for(CI p = begin(); p != end(); ++p)
	{
	s.write((char *) &(p->first), sizeof(p->first));
	s.write((char *) &(p->second), sizeof(p->second));
	}
    }

void SparseVec::read(istream& s)
    {
    int len;
    s.write((char *) &(len), sizeof(len));
    for(int i = 1; i <= len; i++)
	{
	int j;
	Real x;
	s.read((char *) &j, sizeof(j));
	s.read((char *) &x, sizeof(x));
	(*this)[j] = x;
	}
    }

SparseRef::SparseRef(const MatrixRef& M)
    : base(new SparseMatBase(M.Nrows(),M.Ncols())), scale(1.0), transpose(0)
    {
    for (int i = 1; i <= M.Nrows(); i++)
	for (int j = 1; j <= M.Ncols(); j++)
	    {
	    Real x = M(i,j);
	    if (fabs(x) >= SparseVec::thresh)
		Row(i)[j] = x;
	    }
    }

SparseRef& SparseRef::operator+=(const MatrixRef& M)
    {
    OwnCopy();
    for (int i = 1; i <= M.Nrows(); i++)
	for (int j = 1; j <= M.Ncols(); j++)
	    {
	    Real x = M(i,j);
	    if (fabs(x) >= SparseVec::thresh)
		Row(i)[j] += x;
	    }
    return *this;
    }

void SparseRef::PutInMatrix(Matrix& mat) const
    {
    mat.ReduceDimension(Nrows(),Ncols());
    if(Nrows()*Ncols() == 0) return;
    mat = 0.0;
    int i;
    for (i = 1; i <= base->nrows; i++)
	{
	SparseVec& row(Row(i));
	int len = row.Length();
	if(!transpose)
	    for(CI p = begin(); p != end(); ++p)
		mat(i,p->first) = p->second * scale;
	else
	    for(CI p = begin(); p != end(); ++p)
		mat(p->first,i) = p->second * scale;
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
	    SparseVec& row(S.Row(i));
	    for(CI p = row.begin(); p != row.end(); ++p)
		Br += p->second * S.scale * A.Row(p->first);
	    //for(int j = 1; j <= row.Length(); j++)
	//	Br += row.data(j) * S.scale * A.Row(row.index(j));
	    }
    else
	for (int j = 1; j <= S.Ncols(); j++)
	    {
	    Ar << A.Row(j);
	    SparseVec& row(S.Row(j));
	    for(CI p = row.begin(); p != row.end(); ++p)
		B.Row(p->first) += p->second * S.scale * Ar;
	    //for(int i = 1; i <= row.Length(); i++)
		//B.Row(row.index(i)) += row.data(i) * S.scale * Ar;
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
	    SparseVec& row(S.Row(i));
	    for(CI p = row.begin(); p != row.end(); ++p)
		B.Column(p->first) += p->second * S.scale * Ac;
	    //for(j = 1; j <= row.Length(); j++)
	//	B.Column(row.index(j)) += row.data(j) * S.scale * Ac;
	    }
    else
	for (j = 1; j <= S.Ncols(); j++)	// Really columns
	    {
	    VectorRef Bc(B.Column(j));
	    SparseVec& row(S.Row(j));
	    for(CI p = row.begin(); p != row.end(); ++p)
		Bc += p->second * S.scale * A.Column(p->first);
	 //   for(i = 1; i <= row.Length(); i++)
	//	Bc += row.data(i) * S.scale * A.Column(row.index(i));
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
	SparseVec& row(Row(i));
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
	SparseVec& row(S.Row(i));
	row.read(s);
	}
    *this = S;
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
