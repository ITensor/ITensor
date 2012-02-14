#include "sparseitensor.h"
using namespace std;
using boost::format;
using boost::array;

SparseITensor::
SparseITensor()
    :
    r_(0),
    rn_(0),
    scale_(0),
    ur(0)
    { }

SparseITensor::
SparseITensor(const Index& i1)
    :
    r_(1),
    rn_(0),
    ur(0)
    { 
    _construct1(i1);
    }

SparseITensor::
SparseITensor(const Index& i1, Real d)
    :
    r_(1),
    rn_(0),
    scale_(d),
    ur(0)
    { 
    _construct1(i1);
    }

SparseITensor::
SparseITensor(const Index& i1, const Vector& diag)
    :
    diag_(diag),
    r_(1),
    rn_(0),
    ur(0)
    { 
    _construct1(i1);
    }

SparseITensor::
SparseITensor(const Index& i1, const Index& i2, Real d)
    :
    r_(2),
    rn_(0),
    scale_(d),
    ur(0)
    { 
    _construct2(i1,i2);
    }

void SparseITensor::
_construct1(const Index& i1)
	{
	assert(r_ == 1);
	if(i1.m() != 1) 
	    rn_ = 1; 
	index_[1] = i1;
	setUniqueReal();
	}

void SparseITensor::
_construct2(const Index& i1, const Index& i2)
	{
	assert(r_ == 2);
	if(i1.m()==1) 
	    {
	    index_[1] = i2; 
        index_[2] = i1; 
	    rn_ = (i2.m() == 1 ? 0 : 1);
	    }
	else 
	    { 
	    index_[1] = i1; 
        index_[2] = i2; 
	    rn_ = (i2.m() == 1 ? 1 : 2); 
	    }
	setUniqueReal();
	}

void
product(const SparseTensor& S, const ITensor& T, ITensor& res)
    {
    const bool diag_allsame = (diag_.Length()==0);

    res.rn_ = 0;
    res.r_ = 0;
    int alloc_size = 1;

    //
    // tcon[j] = i means that the 
    // jth Index of T is contracted
    // with the ith Index of S
    //
    // tcon[j] = 0 means not contracted
    //
    // (scon is similar but for S)
    //
    array<Index,NMAX+1> tcon,
                        scon,
                        tmap,
                        smap;
    tcon.assign(0);
    scon.assign(0);
    tmap.assign(0);
    smap.assign(0);
    int ncon = 0; //number contracted

    //Analyze contracted Indices
    for(int i = 1; i <= S.r_; ++i)
    for(int j = 1; j <= T.r_; ++j)
        if(S.index_[i] == T.index_[j])
            {
            scon[i] = j;
            tcon[j] = i;

            ++ncon;
            smap[i] = ncon;
            tmap[j] = ncon;
            }

    //Put uncontracted m != 1 Indices
    //of S into res
    for(int i = 1; i <= S.rn_; ++i)
        if(scon[i] != 0)
            {
            res.index_[++res.rn_] = S.index_[i];
            ++res.r_;
            alloc_size *= S.index_[i].m();
            }

    //Put uncontracted m != 1 Indices
    //of T into res
    for(int i = 1; i <= T.rn_; ++i)
        if(tcon[i] != 0)
            {
            res.index_[++res.rn_] = T.index_[i];
            ++res.r_;
            alloc_size *= T.index_[i].m();
            }

    //Put uncontracted m == 1 Indices
    //of S into res
    for(int i = S.rn_+1; i <= S.r_; ++i)
        if(scon[i] != 0)
            {
            res.index_[++res.r_] = S.index_[i];
            }

    //Put uncontracted m == 1 Indices
    //of T into res
    for(int i = T.rn_+1; i <= T.r_; ++i)
        if(tcon[i] != 0)
            {
            res.index_[++res.r_] = T.index_[i];
            }
#ifdef DEBUG
    if(res.r_ != (S.r_+T.r_ - ncon))
        {
        cout << format("res.r_ = %d != (S.r_+T.r_-ncon) = %d")
            % res.r _% (S.r_+T.r_-ncon) << endl;
        Error("Incorrect rank");
        }
#endif

    res.setUniqueReal();

    res.scale_ = S.scale_ * T.scale_;

    //If S has dimension 1
    //it is just a scalar.
    //res may have different m==1 
    //Indices than T, though.
    if(S.rn_ == 0)
        {
        res.p = T.p;
        if(!diag_allsame)
            res *= S.diag_(1);
        return;
        }

    //Allocate a new dat for res if necessary
    if(res.p->count() != 1) { res.p = new ITDat(); }
    res.p->v.ReDimension(alloc_size);

    if(diag_allsame)
        {
        //Determine the size of the diagonal
        int dsize = S.index_[1].m();
        for(int i = 2; i <= S.rn_; ++i)
            if(scon[i] != 0)
                {
                dsize = min(dsize,S.index_[i].m());
                }
        }
    else
        {
        Error("Not implemented");
        }
    }

