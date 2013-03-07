//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#include "itensor.h"
using namespace std;
using boost::format;
using boost::array;
using boost::shared_ptr;
using boost::make_shared;

#ifdef DEBUG
#define ITENSOR_CHECK_NULL if(!p) Error("ITensor is null");
#else
#define ITENSOR_CHECK_NULL
#endif

//
// ITensor
//

//
// ITensor Constructors
//

ITensor::
ITensor()  
    : 
    scale_(1)
    { }


ITensor::
ITensor(Real val) 
    :
    scale_(1)
    { 
    allocate(1);
    p->v = val;
    }

ITensor::
ITensor(const Index& i1) 
    :
    is_(i1),
    scale_(1)
	{ 
    allocate(i1.m());
    }

ITensor::
ITensor(const Index& i1, Real val) 
    :
    is_(i1),
    scale_(1)
	{ 
    allocate(i1.m());
    p->v = val; 
    }

ITensor::
ITensor(const Index& i1, const VectorRef& V) 
    : 
    p(make_shared<ITDat>(V)),
    is_(i1),
    scale_(1)
	{ 
	if(i1.m() != V.Length()) 
	    Error("Mismatch of Index and Vector sizes.");
	}

ITensor::
ITensor(const Index& i1,const Index& i2) 
    :
    is_(i1,i2),
    scale_(1)
	{ 
    allocate(i1.m()*i2.m());
    }
    

ITensor::
ITensor(const Index& i1,const Index& i2,Real a) 
    :
    is_(i1,i2),
    scale_(1)
	{
    allocate(i1.m()*i2.m());
	if(is_.rn() == 2) //then index order is i1, i2
	    {
	    const int nn = min(i1.m(),i2.m());
	    for(int i = 1; i <= nn; ++i) 
		p->v((i-1)*i1.m()+i) = a;
	    }
	else 
	    p->v(1) = a;
	}

ITensor::
ITensor(const Index& i1,const Index& i2,const MatrixRef& M) 
    :
    is_(i1,i2),
    scale_(1)
	{
    allocate(i1.m()*i2.m());
	if(i1.m() != M.Nrows() || i2.m() != M.Ncols()) 
	    Error("Mismatch of Index sizes and matrix.");
	MatrixRef dref; 
	p->v.TreatAsMatrix(dref,i2.m(),i1.m()); 
	dref = M.t();
	}

ITensor::
ITensor(const Index& i1, const Index& i2, const Index& i3,
        const Index& i4, const Index& i5, const Index& i6,
        const Index& i7, const Index& i8)
    :
    scale_(1)
	{
#ifdef DEBUG
    if(i1 == Index::Null())
        Error("i1 is null");
    if(i2 == Index::Null())
        Error("i2 is null");
    if(i3 == Index::Null())
        Error("i3 is null");
#endif
	array<Index,NMAX> ii = {{ i1, i2, i3, i4, i5, i6, i7, i8 }};
	int size = 3;
	while(ii[size] != Index::Null()) ++size;
	int alloc_size = -1; 
    is_ = IndexSet<Index>(ii,size,alloc_size,0);
	allocate(alloc_size);
	}

ITensor::
ITensor(const IndexVal& iv, Real fac) 
    :
    is_(iv.ind),
    scale_(1)
	{ 
    allocate(iv.ind.m());
	p->v(iv.i) = fac; 
	}

ITensor::
ITensor(const IndexVal& iv1, const IndexVal& iv2) 
    :
    is_(iv1.ind,iv2.ind),
    scale_(1)
	{ 
    allocate(iv1.ind.m()*iv2.ind.m());
	p->v((iv2.i-1)*iv1.ind.m()+iv1.i) = 1; 
	}

ITensor::
ITensor(const IndexVal& iv1, const IndexVal& iv2, 
        const IndexVal& iv3, const IndexVal& iv4, 
        const IndexVal& iv5, const IndexVal& iv6, 
        const IndexVal& iv7, const IndexVal& iv8)
    :
    scale_(1)
	{
    //Construct ITensor
    array<Index,NMAX+1> ii = 
        {{ iv1.ind, iv2.ind, iv3.ind, iv4.ind, iv5.ind, 
           iv6.ind, iv7.ind, iv8.ind }};
    int size = 3; 
    while(size < NMAX && ii[size+1] != IndexVal::Null().ind) ++size;
    int alloc_size = -1;
    is_ = IndexSet<Index>(ii,size,alloc_size,0);
    allocate(alloc_size);

    //Assign specified element to 1
    array<int,NMAX+1> iv = 
        {{ iv1.i, iv2.i, iv3.i, iv4.i, iv5.i, iv6.i, iv7.i, iv8.i }};
    array<int,NMAX+1> ja; ja.assign(1);
    for(int k = 1; k <= is_.rn(); ++k) //loop over indices of this ITensor
        for(int j = 0; j < size; ++j)  // loop over the given indices
    if(is_.index(k) == ii[j]) 
        { ja[k] = iv[j]; break; }
    p->v(_ind(is_,ja[1],ja[2],ja[3],ja[4],ja[5],ja[6],ja[7],ja[8])) = 1;
    }

ITensor::
ITensor(const IndexSet<Index>& I) 
    :
    is_(I),
    scale_(1)
	{
	allocate(is_.dim());
	}

ITensor::
ITensor(const std::vector<Index>& I) 
    :
    scale_(1)
	{
    int alloc_size = -1;
    is_ = IndexSet<Index>(I,I.size(),alloc_size,0);
	allocate(alloc_size);
	}

ITensor::
ITensor(const std::vector<Index>& I, const Vector& V) 
    : 
    p(make_shared<ITDat>(V)),
    scale_(1)
	{
    int alloc_size = -1;
    is_ = IndexSet<Index>(I,I.size(),alloc_size,0);
	if(alloc_size != V.Length()) 
	    { Error("incompatible Index and Vector sizes"); }
	}


ITensor::
ITensor(const std::vector<Index>& I, const ITensor& other) 
    : 
    p(other.p), 
    scale_(other.scale_)
	{
    int alloc_size = -1;
    is_ = IndexSet<Index>(I,I.size(),alloc_size,0);
	if(alloc_size != other.vecSize()) 
	    { Error("incompatible Index and ITensor sizes"); }
	}

ITensor::
ITensor(const IndexSet<Index>& I, const ITensor& other) 
    : 
    p(other.p), 
    is_(I),
    scale_(other.scale_)
	{ }

ITensor::
ITensor(const std::vector<Index>& I, const ITensor& other, Permutation P) 
    : 
    scale_(other.scale_)
    {
    int alloc_size = -1;
    is_ = IndexSet<Index>(I,I.size(),alloc_size,0);
    if(alloc_size != other.vecSize()) 
        { Error("incompatible Index and ITensor sizes"); }
    if(P.is_trivial()) { p = other.p; }
    else               { allocate(); other.reshapeDat(P,p->v); }
    }

const ITensor& ITensor::
Complex_1()
    {
    static const ITensor Complex_1_(Index::IndReIm()(1));
    return Complex_1_;
    }

const ITensor& ITensor::
Complex_i()
    {
    static const ITensor Complex_i_(Index::IndReIm()(2));
    return Complex_i_;
    }

ITensor
makeConjTensor()
    {
    ITensor ct(Index::IndReImP());
    ct(Index::IndReImP()(1)) =  1;
    ct(Index::IndReImP()(2)) = -1;
    return ct;
    }

const ITensor& ITensor::
ConjTensor()
    {
    static const ITensor ConjTensor_(makeConjTensor());
    return ConjTensor_;
    }

ITensor
makeComplexProd()
    {
    ITensor prod(Index::IndReIm(),
                 Index::IndReImP(),
                 Index::IndReImPP());

    IndexVal iv0(Index::IndReIm(),1), 
             iv1(Index::IndReImP(),1), 
             iv2(Index::IndReImPP(),1);

    iv0.i = 1; iv1.i = 1; iv2.i = 1; prod(iv0,iv1,iv2) =  1;
    iv0.i = 1; iv1.i = 2; iv2.i = 2; prod(iv0,iv1,iv2) = -1;
    iv0.i = 2; iv1.i = 2; iv2.i = 1; prod(iv0,iv1,iv2) =  1;
    iv0.i = 2; iv1.i = 1; iv2.i = 2; prod(iv0,iv1,iv2) =  1;
    return prod;
    }

const ITensor& ITensor::
ComplexProd()
    {
    static const ITensor ComplexProd_(makeComplexProd());
    return ComplexProd_;
    }

void ITensor::
read(std::istream& s)
    { 
    bool isNull_;
    s.read((char*) &isNull_,sizeof(isNull_));
    if(isNull_) { *this = ITensor(); return; }

    is_.read(s);
    scale_.read(s);
    p = make_shared<ITDat>();
    p->read(s);
    }

void ITensor::
write(std::ostream& s) const 
    { 
    bool isNull_ = isNull();
    s.write((char*) &isNull_,sizeof(isNull_));
    if(isNull_) return;

    is_.write(s);
    scale_.write(s);
    p->write(s);
    }


Real ITensor::
toReal() const 
	{ 
#ifdef DEBUG
    if(this->isNull())
        Error("ITensor is null");
#endif
    if(is_.rn() != 0)
        {
        Print(*this);
        Error("ITensor is not a real scalar");
        }

	try {
	    return p->v(1)*scale_.real(); 
	    }
	catch(const TooBigForReal& e)
	    {
	    cout << "too big for real() in toReal" << endl;
	    cerr << "too big for real() in toReal" << endl;
	    cout << "p->v(1) is " << p->v(1) << endl;
	    cout << "scale is " << scale() << endl;
	    cout << "rethrowing" << endl;
	    throw e;
	    }
	catch(TooSmallForReal)
	    {
	    cout << "warning: too small for real() in toReal" << endl;
	    cerr << "warning: too small for real() in toReal" << endl;
	    cout << "p->v(1) is " << p->v(1) << endl;
	    cout << "scale is " << scale() << endl;
	    return 0.0;
	    }
	return NAN; //shouldn't reach this line
	}

void ITensor::
toComplex(Real& re, Real& im) const 
	{ 
    if(isComplex(*this))
        {
        re = operator()(Index::IndReIm()(1));
        im = operator()(Index::IndReIm()(2));
        }
    else
        {
        re = toReal();
        im = 0;
        }
    }

Real& ITensor::
operator()()
	{ 
    if(is_.rn() != 0)
        {
        std::cerr << format("# given = 0, rn_ = %d\n")%is_.rn();
        Error("Not enough indices (requires all having m!=1)");
        }
    solo(); 
    scaleTo(1);
    return p->v(1);
    }

Real ITensor::
operator()() const
	{ 
    ITENSOR_CHECK_NULL
    if(is_.rn() != 0)
        {
        std::cerr << format("# given = 0, rn_ = %d\n")%is_.rn();
        Error("Not enough indices (requires all having m!=1)");
        }
    return scale_.real()*p->v(1);
    }

Real& ITensor::
operator()(const IndexVal& iv1)
	{
#ifdef DEBUG
    if(is_.rn() > 1) 
        {
        std::cerr << format("# given = 1, rn_ = %d\n")%is_.rn();
        Error("Not enough m!=1 indices provided");
        }
    if(is_[0] != iv1.ind)
        {
        Print(*this);
        Print(iv1);
        Error("Incorrect IndexVal argument to ITensor");
        }
#endif
    solo(); 
    scaleTo(1);
    return p->v(iv1.i);
	}

Real ITensor::
operator()(const IndexVal& iv1) const
	{
    ITENSOR_CHECK_NULL
#ifdef DEBUG
    if(is_.rn() > 1) 
        {
        std::cerr << format("# given = 1, rn() = %d\n")%is_.rn();
        Error("Not enough m!=1 indices provided");
        }
    if(is_[0] != iv1.ind)
        {
        Print(*this);
        Print(iv1);
        Error("Incorrect IndexVal argument to ITensor");
        }
#endif
    return scale_.real()*p->v(iv1.i);
	}

Real& ITensor::
operator()(const IndexVal& iv1, const IndexVal& iv2) 
    {
    solo(); 
    scaleTo(1);
    return p->v(_ind2(iv1,iv2));
    }

Real ITensor::
operator()(const IndexVal& iv1, const IndexVal& iv2) const
    {
    ITENSOR_CHECK_NULL
    return scale_.real()*p->v(_ind2(iv1,iv2));
    }

Real& ITensor::
operator()(const IndexVal& iv1, const IndexVal& iv2, 
           const IndexVal& iv3, const IndexVal& iv4, 
           const IndexVal& iv5,const IndexVal& iv6,
           const IndexVal& iv7,const IndexVal& iv8)
    {
    solo(); 
    scaleTo(1);
    return p->v(_ind8(iv1,iv2,iv3,iv4,iv5,iv6,iv7,iv8));
    }

Real ITensor::
operator()(const IndexVal& iv1, const IndexVal& iv2, 
           const IndexVal& iv3, const IndexVal& iv4,
           const IndexVal& iv5,const IndexVal& iv6,
           const IndexVal& iv7,const IndexVal& iv8) const
    {
    ITENSOR_CHECK_NULL
    return scale_.real()*p->v(_ind8(iv1,iv2,iv3,iv4,iv5,iv6,iv7,iv8));
    }

//#define DO_REWRITE_ASSIGN

/*
void ITensor::
assignFrom(const ITensor& other)
    {
    if(this == &other) return;
    if(fabs(other.is_.uniqueReal() - is_.uniqueReal()) > 1E-12)
        {
        Print(*this); Print(other);
        Error("assignFrom: unique Real not the same"); 
        }
#ifdef DO_REWRITE_ASSIGN
    is_ = other.is_;
    scale_ = other.scale_;
    p = other.p;
#else
    Permutation P; 
    is_.getperm(other.is_,P);
    scale_ = other.scale_;
    if(!p.unique())
        { 
        p = make_shared<ITDat>(); 
        }
    other.reshapeDat(P,p->v);
    DO_IF_PS(++Prodstats::stats().c1;)
#endif
    }
    */


void ITensor::
groupIndices(const array<Index,NMAX+1>& indices, int nind, 
             const Index& grouped, ITensor& res) const
    {
    array<int,NMAX+1> isReplaced; 
    isReplaced.assign(0);

    //Print(*this);

    int tot_m = 1;
    int nn = 0; //number of m != 1 indices
    for(int j = 1; j <= nind; ++j) 
        {
        //cerr << format("indices[%d] = ") % j << indices[j] << "\n";
        const Index& J = indices[j];
        if(J.m() != 1) ++nn;
        tot_m *= J.m();

        bool foundit = false;
        for(int k = 1; k <= r(); ++k) 
            { 
            if(is_.index(k) == J) 
                {
                isReplaced[k] = (J.m() == 1 ? -1 : nn);
                //cerr << format("setting isReplaced[%d] = %d\n ") % k % isReplaced[k];
                foundit = true; 
                break; 
                }
            }
        if(!foundit)
            {
            Print(*this);
            cerr << "Couldn't find Index " << J << " in ITensor.\n";
            Error("bad request for Index to replace");
            }
        }
    if(tot_m != grouped.m()) Error("ITensor::groupAndReplace: \
                                    mismatched index sizes.");

    //Compute rn_ of res
    const int res_rn_ = is_.rn() - nn + (nn == 0 ? 0 : 1);

    vector<Index> nindices; 
    nindices.reserve(r()-nind+1);
    Permutation P;
    int nkept = 0; 
    for(int j = 1; j <= is_.rn(); ++j)
        {
        //cerr << format("isReplaced[%d] = %d\n") % j % isReplaced[j];
        if(isReplaced[j] == 0)
            {
            //cerr << format("Kept index, setting P.from_to(%d,%d)\n") % j % (nkept+1);
            P.from_to(j,++nkept);
            nindices.push_back(is_.index(j)); 
            }
        else
            {
            //cerr << format("Replaced index, setting P.from_to(%d,%d)\n") % j % (res_rn_+isReplaced[j]-1);
            P.from_to(j,res_rn_+isReplaced[j]-1);
            }
        }

    nindices.push_back(grouped);

    for(int j = is_.rn()+1; j <= r(); ++j) 
        if(isReplaced[j] == 0) nindices.push_back(is_.index(j));

    if(nn == 0) 
        res = ITensor(nindices,*this);
    else        
        res = ITensor(nindices,*this,P); 
    }

void ITensor::
tieIndices(const array<Index,NMAX>& indices, int nind,
           const Index& tied)
    {
    if(nind == 0) Error("No indices given");

    const int tm = tied.m();
    
    array<Index,NMAX+1> new_index_;
    new_index_[1] = tied;
    //will count these up below
    int new_r_ = 1;
    int alloc_size = tm;

    array<bool,NMAX+1> is_tied;
    is_tied.assign(false);

    int nmatched = 0;
    for(int k = 1; k <= r(); ++k)
        {
        const Index& K = is_.index(k);
        for(int j = 0; j < nind; ++j)
        if(K == indices[j]) 
            { 
            if(indices[j].m() != tm)
                Error("Tied indices must have matching m's");
            is_tied[k] = true;

            ++nmatched;

            break;
            }

        if(!is_tied[k])
            {
            new_index_[++new_r_] = K;
            alloc_size *= K.m();
            }
        }

    //Check that all indices were found
    if(nmatched != nind)
        {
        Print(*this);
        cout << "indices = " << endl;
        for(int j = 0; j < nind; ++j)
            cout << indices[j] << endl;
        Error("Couldn't find Index to tie");
        }

    IndexSet<Index> new_is_(new_index_,new_r_,alloc_size,1);

    //If tied indices have m==1, no work
    //to do; just replace indices
    if(tm == 1)
        {
        is_.swap(new_is_);
        return;
        }

    Counter nc(new_is_);

    //Set up ii pointers to link
    //elements of res to appropriate
    //elements of *this
    array<int*,NMAX+1> ii;
    int n = 2;
    for(int j = 1; j <= r(); ++j)
        {
        if(is_tied[j])
            ii[j] = &(nc.i[1]);
        else
            ii[j] = &(nc.i[n++]);
        }

    int one = 1;
    for(int j = r()+1; j <= NMAX; ++j)
        ii[j] = &one;
    
    //Create the new dat
    shared_ptr<ITDat> np = make_shared<ITDat>(alloc_size);
    Vector& resdat = np->v;

    const Vector& thisdat = p->v;
    for(; nc.notDone(); ++nc)
        {
        resdat(nc.ind) =
        thisdat(_ind(is_,*ii[1],*ii[2],
                         *ii[3],*ii[4],
                         *ii[5],*ii[6],
                         *ii[7],*ii[8]));
        }

    is_.swap(new_is_);
    p.swap(np);

    } //ITensor::tieIndices

void ITensor::
tieIndices(const Index& i1, const Index& i2,
           const Index& tied)
    {
    array<Index,NMAX> inds =
        {{ i1, i2, 
           Index::Null(), Index::Null(), 
           Index::Null(), Index::Null(), 
           Index::Null(), Index::Null() }};

    tieIndices(inds,2,tied);
    }

void ITensor::
tieIndices(const Index& i1, const Index& i2,
           const Index& i3,
           const Index& tied)
    {
    array<Index,NMAX> inds =
        {{ i1, i2, i3,
           Index::Null(), Index::Null(), 
           Index::Null(), Index::Null(), Index::Null() }};

    tieIndices(inds,3,tied);
    }

void ITensor::
tieIndices(const Index& i1, const Index& i2,
           const Index& i3, const Index& i4,
           const Index& tied)
    {
    array<Index,NMAX> inds =
        {{ i1, i2, i3, i4,
           Index::Null(), Index::Null(), 
           Index::Null(), Index::Null() }};

    tieIndices(inds,4,tied);
    }

void ITensor::
trace(const array<Index,NMAX>& indices, int nind)
    {
    if(nind == 0) Error("No indices given");

    const int tm = indices[0].m();
    
    array<Index,NMAX+1> new_index_;

    //will count these up below
    int new_r_ = 0;
    int alloc_size = 1;

    array<bool,NMAX+1> traced;
    traced.assign(false);

    int nmatched = 0;
    for(int k = 1; k <= r(); ++k)
        {
        const Index& K = is_.index(k);
        for(int j = 0; j < nind; ++j)
        if(K == indices[j]) 
            { 
#ifdef DEBUG
            if(indices[j].m() != tm)
                {
                Print((*this));
                Print(K);
                Print(indices[j]);
                Error("Traced indices must all have same m's");
                }
#endif
            traced[k] = true;

            ++nmatched;

            break;
            }

        if(!traced[k])
            {
            new_index_[++new_r_] = K;
            alloc_size *= K.m();
            }
        }

    //Check that all indices were found
    if(nmatched != nind)
        {
        Print(*this);
        cout << "indices = " << endl;
        for(int j = 0; j < nind; ++j)
            cout << indices[j] << endl;
        Error("Couldn't find Index to trace");
        }

    IndexSet<Index> new_is_(new_index_,new_r_,alloc_size,1);

    //If traced indices have m==1, no work
    //to do; just replace indices
    if(tm == 1)
        {
        is_.swap(new_is_);
        return;
        }

    Counter nc(new_is_);

    //Set up ii pointers to link
    //elements of res to appropriate
    //elements of *this
    int trace_ind = 0;
    array<int*,NMAX+1> ii;
    int n = 1;
    for(int j = 1; j <= r(); ++j)
        {
        if(traced[j])
            ii[j] = &(trace_ind);
        else
            ii[j] = &(nc.i[n++]);
        }

    int one = 1;
    for(int j = r()+1; j <= NMAX; ++j)
        ii[j] = &one;
    
    //Create the new dat
    shared_ptr<ITDat> np = make_shared<ITDat>(alloc_size);
    Vector& resdat = np->v;

    const Vector& thisdat = p->v;
    for(; nc.notDone(); ++nc)
        {
        Real newval = 0;
        for(trace_ind = 1; trace_ind <= tm; ++trace_ind)
            {
            newval += 
            thisdat(_ind(is_,*ii[1],*ii[2],
                             *ii[3],*ii[4],
                             *ii[5],*ii[6],
                             *ii[7],*ii[8]));
            }
        resdat(nc.ind) = newval;
        }

    is_.swap(new_is_);
    p.swap(np);

    } //ITensor::trace

void ITensor::
trace(const Index& i1)
    {
    array<Index,NMAX> inds =
        {{ i1, Index::Null(), 
           Index::Null(), Index::Null(), 
           Index::Null(), Index::Null(), 
           Index::Null(), Index::Null() }};

    trace(inds,1);
    }

void ITensor::
trace(const Index& i1, const Index& i2)
    {
    array<Index,NMAX> inds =
        {{ i1, i2, 
           Index::Null(), Index::Null(), 
           Index::Null(), Index::Null(), 
           Index::Null(), Index::Null() }};

    trace(inds,2);
    }

void ITensor::
trace(const Index& i1, const Index& i2, const Index& i3)
    {
    array<Index,NMAX> inds =
        {{ i1, i2, i3,
           Index::Null(), Index::Null(), 
           Index::Null(), Index::Null(), Index::Null() }};

    trace(inds,3);
    }

void ITensor::
trace(const Index& i1, const Index& i2,
      const Index& i3, const Index& i4)
    {
    array<Index,NMAX> inds =
        {{ i1, i2, i3, i4,
           Index::Null(), Index::Null(), 
           Index::Null(), Index::Null() }};

    trace(inds,4);
    }

void ITensor::
expandIndex(const Index& small, const Index& big, int start)
    {
    assert(small.m() <= big.m());
    assert(start < big.m());

    vector<Index> indices; 
    indices.reserve(is_.r());
    int w = -1;
    for(int j = 1; j <= r(); ++j)
        {
        if(is_.index(j) == small)
            {
            w = j;
            indices.push_back(big);
            }
        else 
            {
            indices.push_back(is_.index(j));
            }
        }

    if(w == -1)
        {
        Print(*this);
        Print(small);
        Error("couldn't find index");
        }

    ITensor res(indices);
    res.scale_ = scale_;

    //Big Index not guaranteed
    //to remain at position w
    //e.g. if some m==1 Indices
    //get moved to the back
    w = res.is_.findindex(big);

    array<int,NMAX+1> inc;
    //Make sure all other inc's are zero
    inc.assign(0);
    inc.at(w) = start;

    Counter c(is_);

    const Vector& thisdat = p->v;
    Vector& resdat = res.p->v;
    for(; c.notDone(); ++c)
        {
        resdat(_ind(res.is_,c.i[1]+inc[1],c.i[2]+inc[2],
                            c.i[3]+inc[3],c.i[4]+inc[4],
                            c.i[5]+inc[5],c.i[6]+inc[6],
                            c.i[7]+inc[7],c.i[8]+inc[8]))
        = thisdat(c.ind);
        }

    this->swap(res);
    }

int ITensor::
vecSize() const 
    { 
    return (p.get() == 0 ? 0 : p->v.Length()); 
    }

void ITensor::
assignToVec(VectorRef v) const
    {
    if(p->v.Length() != v.Length()) 
        Error("ITensor::assignToVec bad size");
    if(scale_.isRealZero()) 
        {
        v *= 0;
        return;
        }
    ITENSOR_CHECK_NULL
    v = p->v;
    v *= scale_.real();
    }

void ITensor::
assignFromVec(const VectorRef& v)
    {
    ITENSOR_CHECK_NULL
    if(p->v.Length() != v.Length()) 
	Error("ITensor::assignToVec bad size");
    scale_ = 1;
    if(!p.unique())
        { 
        p = make_shared<ITDat>();
        }
    p->v = v;
    }

void ITensor::
reshapeDat(const Permutation& P, Vector& rdat) const
    {
    ITENSOR_CHECK_NULL

    const Vector& thisdat = p->v;

    if(P.is_trivial())
        {
        rdat = thisdat;
        return;
        }

    rdat.ReDimension(thisdat.Length());
    rdat = 0;

    const Permutation::int9& ind = P.ind();

    //Make a counter for thisdat
    Counter c(is_);
    array<int,NMAX+1> n;
    for(int j = 1; j <= c.rn; ++j) n[ind[j]] = c.n[j];

    //Special case loops
#define Loop6(q,z,w,k,y,s) {for(int i1 = 1; i1 <= n[1]; ++i1) for(int i2 = 1; i2 <= n[2]; ++i2)\
	for(int i3 = 1; i3 <= n[3]; ++i3) for(int i4 = 1; i4 <= n[4]; ++i4) for(int i5 = 1; i5 <= n[5]; ++i5)\
    for(int i6 = 1; i6 <= n[6]; ++i6)\
    rdat( (((((i6-1)*n[5]+i5-1)*n[4]+i4-1)*n[3]+i3-1)*n[2]+i2-1)*n[1]+i1 ) =\
    thisdat( (((((s-1)*c.n[5]+y-1)*c.n[4]+k-1)*c.n[3]+w-1)*c.n[2]+z-1)*c.n[1]+q ); return; }

#define Loop5(q,z,w,k,y) {for(int i1 = 1; i1 <= n[1]; ++i1) for(int i2 = 1; i2 <= n[2]; ++i2)\
	for(int i3 = 1; i3 <= n[3]; ++i3) for(int i4 = 1; i4 <= n[4]; ++i4) for(int i5 = 1; i5 <= n[5]; ++i5)\
    rdat( ((((i5-1)*n[4]+i4-1)*n[3]+i3-1)*n[2]+i2-1)*n[1]+i1 ) = thisdat( ((((y-1)*c.n[4]+k-1)*c.n[3]+w-1)*c.n[2]+z-1)*c.n[1]+q ); return; }

#define Loop4(q,z,w,k) {for(int i1 = 1; i1 <= n[1]; ++i1)  for(int i2 = 1; i2 <= n[2]; ++i2)\
	for(int i3 = 1; i3 <= n[3]; ++i3) for(int i4 = 1; i4 <= n[4]; ++i4)\
	rdat( (((i4-1)*n[3]+i3-1)*n[2]+i2-1)*n[1]+i1 ) = thisdat( (((k-1)*c.n[3]+w-1)*c.n[2]+z-1)*c.n[1]+q ); return; }

#define Loop3(q,z,w) {for(int i1 = 1; i1 <= n[1]; ++i1)  for(int i2 = 1; i2 <= n[2]; ++i2)\
	for(int i3 = 1; i3 <= n[3]; ++i3) rdat( ((i3-1)*n[2]+i2-1)*n[1]+i1 ) = thisdat( ((w-1)*c.n[2]+z-1)*c.n[1]+q ); return; }

#define Bif3(a,b,c) if(ind[1] == a && ind[2] == b && ind[3] == c)

#define Bif4(a,b,c,d) if(ind[1] == a && ind[2] == b && ind[3] == c && ind[4] == d)

#define Bif5(a,b,c,d,e) if(ind[1] == a && ind[2] == b && ind[3] == c && ind[4]==d && ind[5] == e)

#define Bif6(a,b,c,d,e,g) if(ind[1] == a && ind[2] == b && ind[3] == c && ind[4]==d && ind[5] == e && ind[6] == g)

    if(is_.rn() == 2 && ind[1] == 2 && ind[2] == 1)
        {
        MatrixRef xref; 
        thisdat.TreatAsMatrix(xref,c.n[2],c.n[1]);
        rdat = Matrix(xref.t()).TreatAsVector();
        return; 
        }
    else if(is_.rn() == 3)
        {
        DO_IF_PS(int idx = ((ind[1]-1)*3+ind[2]-1)*3+ind[3]; Prodstats::stats().perms_of_3[idx] += 1; )
        //Arranged loosely in order of frequency of occurrence
        Bif3(2,1,3) Loop3(i2,i1,i3)
        Bif3(2,3,1) Loop3(i2,i3,i1) //cyclic
        Bif3(3,1,2) Loop3(i3,i1,i2)
        //Bif3(1,3,2) Loop3(i1,i3,i2)
        //Bif3(3,2,1) Loop3(i3,i2,i1)
        }
    else if(is_.rn() == 4)
        {
        DO_IF_PS(int idx = (((ind[1]-1)*4+ind[2]-1)*4+ind[3]-1)*4+ind[4]; Prodstats::stats().perms_of_4[idx] += 1; )
        //Arranged loosely in order of frequency of occurrence
        Bif4(1,2,4,3) Loop4(i1,i2,i4,i3)
        Bif4(1,3,2,4) Loop4(i1,i3,i2,i4)
        Bif4(2,3,1,4) Loop4(i2,i3,i1,i4)
        Bif4(2,3,4,1) Loop4(i2,i3,i4,i1) //cyclic
        Bif4(1,4,2,3) Loop4(i1,i4,i2,i3)
        Bif4(2,1,3,4) Loop4(i2,i1,i3,i4)
        Bif4(2,1,4,3) Loop4(i2,i1,i4,i3)
        Bif4(3,4,1,2) Loop4(i3,i4,i1,i2)
        }
    else if(is_.rn() == 5)
        {
        DO_IF_PS(int idx = ((((ind[1]-1)*5+ind[2]-1)*5+ind[3]-1)*5+ind[4]-1)*5+ind[5]; Prodstats::stats().perms_of_5[idx] += 1; )
        //Arranged loosely in order of frequency of occurrence
        Bif5(3,1,4,5,2) Loop5(i3,i1,i4,i5,i2)
        Bif5(1,4,2,5,3) Loop5(i1,i4,i2,i5,i3)
        Bif5(1,4,2,3,5) Loop5(i1,i4,i2,i3,i5)
        Bif5(3,1,4,2,5) Loop5(i3,i1,i4,i2,i5)
        Bif5(2,4,1,3,5) Loop5(i2,i4,i1,i3,i5)
        Bif5(2,4,3,5,1) Loop5(i2,i4,i3,i5,i1)
        Bif5(3,1,4,5,2) Loop5(i3,i1,i4,i5,i2)
        Bif5(3,4,1,2,5) Loop5(i3,i4,i1,i2,i5)
        Bif5(2,1,3,4,5) Loop5(i2,i1,i3,i4,i5)
        Bif5(2,3,4,5,1) Loop5(i2,i3,i4,i5,i1)
        Bif5(2,3,4,1,5) Loop5(i2,i3,i4,i1,i5)
        Bif5(2,3,1,4,5) Loop5(i2,i3,i1,i4,i5)
        Bif5(2,3,4,1,5) Loop5(i2,i3,i4,i1,i5)
        Bif5(3,4,1,5,2) Loop5(i3,i4,i1,i5,i2)
        Bif5(5,1,4,2,3) Loop5(i5,i1,i4,i2,i3)
        }
    else if(is_.rn() == 6)
        {
        DO_IF_PS(int idx = (((((ind[1]-1)*6+ind[2]-1)*6+ind[3]-1)*6+ind[4]-1)*6+ind[5]-1)*6+ind[6]; Prodstats::stats().perms_of_6[idx] += 1; )
        //Arranged loosely in order of frequency of occurrence
        Bif6(2,4,1,3,5,6) Loop6(i2,i4,i1,i3,i5,i6)
        Bif6(1,4,2,3,5,6) Loop6(i1,i4,i2,i3,i5,i6)
        Bif6(2,4,1,5,3,6) Loop6(i2,i4,i1,i5,i3,i6)
        Bif6(1,2,4,5,3,6) Loop6(i1,i2,i4,i5,i3,i6)
        Bif6(3,4,1,5,6,2) Loop6(i3,i4,i1,i5,i6,i2)
        }
    DO_IF_PS(Prodstats::stats().c4 += 1;)

    //The j's are pointers to the i's of xdat's Counter,
    //but reordered in a way appropriate for rdat
    array<int*,NMAX+1> j;
    for(int k = 1; k <= NMAX; ++k) { j[ind[k]] = &(c.i[k]); }

    //Catch-all loops that work for any tensor
    switch(c.rn)
    {
    case 2:
        for(; c.notDone(); ++c)
            {
            rdat((*j[2]-1)*n[1]+*j[1])
                = thisdat(c.ind);
            }
        return;
    case 3:
        for(; c.notDone(); ++c)
            {
            rdat(((*j[3]-1)*n[2]+*j[2]-1)*n[1]+*j[1])
                = thisdat(c.ind);
            }
        return;
    case 4:
        for(; c.notDone(); ++c)
            {
            rdat((((*j[4]-1)*n[3]+*j[3]-1)*n[2]+*j[2]-1)*n[1]+*j[1])
                = thisdat(c.ind);
            }
        return;
    case 5:
        for(; c.notDone(); ++c)
            {
            rdat(((((*j[5]-1)*n[4]+*j[4]-1)*n[3]+*j[3]-1)*n[2]+*j[2]-1)*n[1]+*j[1])
                = thisdat(c.ind);
            }
        return;
    case 6:
        for(; c.notDone(); ++c)
            {
            rdat((((((*j[6]-1)*n[5]+*j[5]-1)*n[4]+*j[4]-1)*n[3]+*j[3]-1)*n[2]+*j[2]-1)*n[1]+*j[1])
                = thisdat(c.ind);
            }
        return;
    case 7:
        for(; c.notDone(); ++c)
            {
            rdat(((((((*j[7]-1)*n[6]+*j[6]-1)*n[5]+*j[5]-1)*n[4]+*j[4]-1)*n[3]+*j[3]-1)*n[2]+*j[2]-1)*n[1]+*j[1])
                = thisdat(c.ind);
            }
        return;
    default:
        for(; c.notDone(); ++c)
            {
            rdat((((((((*j[8]-1)*n[7]+*j[7]-1)*n[6]+*j[6]-1)*n[5]+*j[5]-1)*n[4]+*j[4]-1)*n[3]+*j[3]-1)*n[2]+*j[2]-1)*n[1]+*j[1])
                = thisdat(c.ind);
            }
        return;
    } //switch(c.rn)

    } // ITensor::reshapeDat

void ITensor::
reshapeDat(const Permutation& P)
    {
    if(P.is_trivial()) return;
    solo();
    Vector newdat;
    this->reshapeDat(P,newdat);
    p->v = newdat;
    }


void ITensor::
swap(ITensor& other)
    {
    p.swap(other.p);
    is_.swap(other.is_);
    scale_.swap(other.scale_);
    }

Real* ITensor::
datStart() const
    {
    if(!p) Error("ITensor is null");
    return p->v.First();
    }

void ITensor::
randomize() 
    { 
    solo(); 
    p->v.Randomize(); 
    }

void ITensor::
conj() 
    { 
    if(!isComplex(*this)) return; 
    prime(ReIm,1);
    operator/=(ITensor::ConjTensor()); 
    prime(ReIm,-1);
    }

Real ITensor::
sumels() const 
    { return p->v.sumels() * scale_.real0(); }

Real ITensor::
norm() const 
    { 
    if(scale_.isTooBigForReal())
        {
        throw TooBigForReal("Scale too large for real in ITensor::norm()");
        }
    //If scale_ is too small to be converted to Real,
    //real0 method will return 0.0
    return fabs(Norm(p->v) * scale_.real0()); 
    }

void ITensor::
pseudoInvert(Real cutoff)
    {
    solo();
    //Invert scale_
    scale_.pow(-1);

    //Invert elems
    for(int j = 1; j <= p->v.Length(); ++j)
        {
        Real elem = p->v(j);
        p->v(j) = (fabs(elem) <= cutoff ? 0 : 1./elem);
        }
    }

void ITensor::
scaleOutNorm()
    {
    Real f = Norm(p->v);
    //If norm already 1 return so
    //we don't have to call solo()
    if(fabs(f-1) < 1E-12) return;

    if(f != 0) 
        { 
        solo();
        p->v *= 1./f; 
        scale_ *= f; 
        }
    else
        {
        scale_ = LogNumber(0.0);
        }
    }

void ITensor::
scaleTo(const LogNumber& newscale)
    {
    if(newscale.sign() == 0) 
        Error("Trying to scale an ITensor to a 0 scale");
    if(scale_ == newscale) return;
    solo();
    scale_ /= newscale;
    p->v *= scale_.real0();
    scale_ = newscale;
    }

void ITensor::
print(std::string name,Printdat pdat) const 
    { 
    Global::printdat() = (pdat==ShowData); 
    std::cerr << "\n" << name << " =\n" << *this << "\n"; 
    Global::printdat() = false; 
    }

void ITensor::
allocate(int dim) 
    { 
    p = make_shared<ITDat>(dim); 
    }

void ITensor::
allocate() 
    { 
    p = make_shared<ITDat>(); 
    }

void ITensor::
solo()
	{
    ITENSOR_CHECK_NULL
    if(!p.unique())
        { 
        VectorRef oldv(p->v);
        p = make_shared<ITDat>();
        p->v = oldv;
        }
	}

int
_ind(const IndexSet<Index>& is,
     int i1, int i2, int i3, int i4, 
     int i5, int i6, int i7, int i8)
    {
    switch(is.rn())
    {
    case 0:
        return (1);
    case 1:
        return (i1);
    case 2:
        return ((i2-1)*is[0].m()+i1);
    case 3:
        return (((i3-1)*is[1].m()+i2-1)*is[0].m()+i1);
    case 4:
        return ((((i4-1)*is[2].m()+i3-1)*is[1].m()+i2-1)*is[0].m()+i1);
    case 5:
        return (((((i5-1)*is[3].m()+i4-1)*is[2].m()+i3-1)*is[1].m()+i2-1)
                        *is[0].m()+i1);
    case 6:
        return ((((((i6-1)*is[4].m()+i5-1)*is[3].m()+i4-1)*is[2].m()+i3-1)
                        *is[1].m()+i2-1)*is[0].m()+i1);
    case 7:
        return (((((((i7-1)*is[5].m()+i6-1)*is[4].m()+i5-1)*is[3].m()+i4-1)
                        *is[2].m()+i3-1)*is[1].m()+i2-1)*is[0].m()+i1);
    case 8:
        return ((((((((i8-1)*is[6].m()+i7-1)*is[5].m()+i6-1)*is[4].m()+i5-1)
                        *is[3].m()+i4-1)*is[2].m()+i3-1)*is[1].m()+i2-1)*is[0].m()+i1);
    } //switch
    Error("_ind: Failed switch case");
    return 1;
    }


int ITensor::
_ind2(const IndexVal& iv1, const IndexVal& iv2) const
    {
    if(is_.rn() > 2) 
        {
        std::cerr << format("# given = 2, rn_ = %d\n")%is_.rn();
        Error("Not enough m!=1 indices provided");
        }
    if(is_[0] == iv1.ind && is_[1] == iv2.ind)
        return ((iv2.i-1)*is_[0].m()+iv1.i);
    else if(is_[0] == iv2.ind && is_[1] == iv1.ind)
        return ((iv1.i-1)*is_[0].m()+iv2.i);
    else
        {
        Print(*this);
        Print(iv1);
        Print(iv2);
        Error("Incorrect IndexVal argument to ITensor");
        return 1;
        }
    }

int ITensor::
_ind8(const IndexVal& iv1, const IndexVal& iv2, 
      const IndexVal& iv3, const IndexVal& iv4,
      const IndexVal& iv5,const IndexVal& iv6,
      const IndexVal& iv7,const IndexVal& iv8) const
    {
    array<const IndexVal*,NMAX> iv = 
        {{ &iv1, &iv2, &iv3, &iv4, &iv5, &iv6, &iv7, &iv8 }};
    array<int,NMAX> ja; ja.assign(1);
    //Loop over the given IndexVals
    int nn = 0;
    for(int j = 0; j < is_.r(); ++j)
        {
        const IndexVal& J = *iv[j];
        if(J.ind == Index::Null()) break;
        if(J.ind.m() != 1) ++nn;
        if(J.i == 1) continue;
        //Loop over indices of this ITensor
        for(int k = 0; k < is_.r(); ++k)
            {
            if(is_[k] == J.ind)
                {
                ja[k] = J.i;
                goto next;
                }
            }
        //Either didn't find index
        Print((*this));
        Print(J);
        Error("Extra/incorrect IndexVal argument to ITensor");
        //Or did find it
        next:;
        }

    if(nn != is_.rn())
        {
        Error("Too few m > 1 indices provided");
        }

    return _ind(is_,ja[0],ja[1],ja[2],ja[3],ja[4],ja[5],ja[6],ja[7]);
    }


//
// Analyzes two ITensors to determine
// how they should be multiplied:
// how many indices do they share?
// which indices are common? etc.
//
struct ProductProps
    {
    ProductProps(const ITensor& L, const ITensor& R);

    //arrays specifying which indices match
    array<bool,NMAX+1> contractedL, contractedR; 

    int nsamen, //number of m !=1 indices that match
        cdim,   //total dimension of contracted inds
        odimL,  //outer (total uncontracted) dim of L
        odimR,  //outer (total uncontracted) dim of R
        lcstart, //where L's contracted inds start
        rcstart; //where R's contracted inds start

    //Permutations that move all matching m!=1
    //indices pairwise to the front 
    Permutation pl, pr;

    //Permutation which, if applied, will make
    //contracted indices of R match order of L
    //Permutation matchL;

    };

ProductProps::
ProductProps(const ITensor& L, const ITensor& R) 
    :
    nsamen(0), 
    cdim(1), 
    odimL(-1), 
    odimR(-1),
    lcstart(100), 
    rcstart(100)
    {
    for(int j = 1; j <= NMAX; ++j) 
        contractedL[j] = contractedR[j] = false;

    for(int j = 1; j <= L.is_.rn(); ++j)
	for(int k = 1; k <= R.is_.rn(); ++k)
	    if(L.is_.index(j) == R.is_.index(k))
		{
		if(j < lcstart) lcstart = j;
        if(k < rcstart) rcstart = k;

		++nsamen;
		pl.from_to(j,nsamen);
		pr.from_to(k,nsamen);

		contractedL[j] = contractedR[k] = true;

        cdim *= L.is_.index(j).m();

        //matchL.from_to(k,j-lcstart+1);
		}
    //Finish making pl
    int q = nsamen;
    for(int j = 1; j <= L.is_.rn(); ++j)
        if(!contractedL[j]) pl.from_to(j,++q);
    //Finish making pr and matchL
    q = nsamen;
    for(int j = 1; j <= R.is_.rn(); ++j)
        if(!contractedR[j]) 
            {
            ++q;
            pr.from_to(j,q);
            //matchL.from_to(j,q);
            }

    odimL = L.p->v.Length()/cdim;
    odimR = R.p->v.Length()/cdim;
    }

//Converts ITensor dats into MatrixRef's that can be multiplied as rref*lref
//contractedL/R[j] == true if L/R.indexn(j) contracted
void 
toMatrixProd(const ITensor& L, const ITensor& R, ProductProps& props,
             MatrixRefNoLink& lref, MatrixRefNoLink& rref, 
             bool& L_is_matrix, bool& R_is_matrix, bool doReshape)
    {
    assert(L.p != 0);
    assert(R.p != 0);
    const Vector &Ldat = L.p->v, &Rdat = R.p->v;

    //Initially, assume both L & R are matrix-like
    L_is_matrix = true, 
    R_is_matrix = true;

    if(props.nsamen != 0)
        {
        //Check that contracted inds are contiguous
        //and in the same order
        for(int i = 0; i < props.nsamen; ++i) 
            {
            if(!props.contractedL[props.lcstart+i] ||
                props.pl.dest(props.lcstart+i) != (i+1)) 
                {
                L_is_matrix = false;
                }
            if(!props.contractedR[props.rcstart+i] ||
                props.pr.dest(props.rcstart+i) != (i+1)) 
                { 
                R_is_matrix = false; 
                }
            }
        //Check that contracted inds are all at beginning or end of _indexn
        if(!(props.contractedL[1] || props.contractedL[L.is_.rn()])) 
            {
            L_is_matrix = false; 
            }
        if(!(props.contractedR[1] || props.contractedR[R.is_.rn()]))
            {
            R_is_matrix = false; 
            }
        }

    if(!doReshape && (!L_is_matrix || !R_is_matrix))
        {
        return;
        }

    if(L_is_matrix)  
        {
        if(props.contractedL[1]) 
            { 
            Ldat.TreatAsMatrix(lref,props.odimL,props.cdim); 
            lref.ApplyTrans(); 
            }
        else 
            { 
            Ldat.TreatAsMatrix(lref,props.cdim,props.odimL); 
            }
        }
    else //L not matrix, need to reshape to make lref
        {
        Vector lv; L.reshapeDat(props.pl,lv);
        lv.TreatAsMatrix(lref,props.odimL,props.cdim); lref.ApplyTrans();
        }

    if(R_is_matrix) 
        {
        if(props.contractedR[1]) 
            { Rdat.TreatAsMatrix(rref,props.odimR,props.cdim); }
        else                    
            { Rdat.TreatAsMatrix(rref,props.cdim,props.odimR); rref.ApplyTrans(); }
        }
    else //R not matrix, need to reshape to make rref
        {
        Vector rv; R.reshapeDat(props.pr,rv);
        rv.TreatAsMatrix(rref,props.odimR,props.cdim);
        }

#ifdef COLLECT_PRODSTATS
    if(L.is_.rn() > R.is_.rn()) 
        {
        ++(Prodstats::stats().global[std::make_pair(L.is_.rn(),R.is_.rn())]);
        }
    else 
        {
        ++(Prodstats::stats().global[std::make_pair(R.is_.rn(),L.is_.rn())]);
        }
    ++Prodstats::stats().total;
    if(L_is_matrix) ++Prodstats::stats().did_matrix;
    if(R_is_matrix) ++Prodstats::stats().did_matrix;
#endif
    }


//Non-contracting product: Cikj = Aij Bkj (no sum over j)
ITensor& ITensor::
operator/=(const ITensor& other)
    {
    if(this == &other)
        {
        ITensor cp_oth(other);
        return operator/=(cp_oth);
        }

    if(hasindexn(Index::IndReIm())   && other.hasindexn(Index::IndReIm()) &&
      !hasindexn(Index::IndReImP())  && !other.hasindexn(Index::IndReImP()) &&
      !hasindexn(Index::IndReImPP()) && !other.hasindexn(Index::IndReImPP()))
        {
        prime(ReIm,1);
        operator/=(primed(other,ReIm,2));
        operator*=(ComplexProd());
        return *this;
        }

    //------------------------------------------------------------------
    //Handle m==1 Indices: set union

    if(other.is_.rn() == 0)
        {
        scale_ *= other.scale_;
        scale_ *= other.p->v(1);
        for(int j = other.is_.rn()+1; j <= other.r(); ++j)
            {
            const Index& J = other.is_.index(j);
            if(!is_.hasindex1(J))
                is_.addindex1(J);
            }
        return *this;
        }
    else if(is_.rn() == 0)
        {
        scale_ *= other.scale_;
        scale_ *= p->v(1);
        p = other.p;
        IndexSet<Index> new_is(other.is_);
        for(int j = 1; j <= r(); ++j) 
            {
            const Index& J = is_.index(j);
            if(!new_is.hasindex1(J))
                new_is.addindex1(J);
            }
        is_.swap(new_is);
        return *this;
        }

    ProductProps props(*this,other);
    MatrixRefNoLink lref, rref;
    bool L_is_matrix,R_is_matrix;
    toMatrixProd(*this,other,props,lref,rref,L_is_matrix,R_is_matrix);

    if(!p.unique()) allocate();

    Vector& thisdat = p->v; 
    
    const int ni = lref.Ncols(), nj = lref.Nrows(), nk = rref.Nrows();
    Matrix L(lref), R(rref);
    thisdat.ReDimension(ni*nj*nk);
    
    for(int j = 1; j <= nj; ++j) 
    for(int k = 1; k <= nk; ++k) 
    for(int i = 1; i <= ni; ++i)
        { thisdat(((j-1)*nk+k-1)*ni+i) =  R(k,j) * L(j,i); }


    IndexSet<Index> new_index;

    //Handle m!=1 indices

    for(int j = 0; j < is_.rn(); ++j)
        if(!props.contractedL[j+1]) 
            new_index.addindexn(this->is_[j]);

    for(int j = 0; j < other.is_.rn(); ++j)
        if(!props.contractedR[j+1]) 
            new_index.addindexn(other.is_[j]);

    for(int j = 0; j < is_.rn(); ++j)
        if(props.contractedL[j+1])  
            new_index.addindexn(this->is_[j]);

    //Handle m==1 indices

    for(int j = is_.rn(); j < r(); ++j) 
        new_index.addindex1(is_[j]);

    for(int j = other.is_.rn()+1; j <= other.r(); ++j)
        {
        const Index& J = other.is_.index(j);
        if(!is_.hasindex1(J)) 
            new_index.addindex1(J);
        }

    is_.swap(new_index);
    
    scale_ *= other.scale_;

    scaleOutNorm();

    return *this;
    }


#define OLD_DIRECT_MULTIPLY

#ifdef OLD_DIRECT_MULTIPLY

int inline
ind4(int i4, int m3, int i3, int m2, int i2, int m1, int i1)
    {
    return (((i4-1)*m3+i3-1)*m2+i2-1)*m1+i1;
    }

void ITensor::
directMultiply(const ITensor& other, ProductProps& props, 
               IndexSet<Index>& new_index)
    {
    int am[NMAX+1], bm[NMAX+1], mcon[NMAX+1], mnew[NMAX+1];
    int *pa[NMAX+1], *pb[NMAX+1];
    int one = 1;
    for(int j = 1; j <= NMAX; ++j)
        {
        //Set *pa[j],*pb[j],mcon[j] and mnew[j]
        // to 1 unless set otherwise below
        pa[j] = pb[j] = &one, mcon[j] = mnew[j] = 1;
        am[j] = is_[j-1].m();
        bm[j] = other.is_[j-1].m();
        }
    int new_rn = 0;
    int icon[NMAX+1], inew[NMAX+1];
    for(int j = 0; j < this->is_.rn(); ++j)
        if(!props.contractedL[j+1]) 
            {
            new_index.addindexn(is_[j]);
            ++new_rn;
            mnew[new_rn] = am[j+1];
            pa[j+1] = inew + new_rn;
            }
        else
            {
            mcon[props.pl.dest(j+1)] = am[j+1];
            pa[j+1] = icon + props.pl.dest(j+1);
            }

    for(int j = 0; j < other.is_.rn(); ++j)
        if(!props.contractedR[j+1]) 
            {
            new_index.addindexn(other.is_[j]);
            ++new_rn;
            mnew[new_rn] = bm[j+1];
            pb[j+1] = inew + new_rn;
            }
        else
            {
            mcon[props.pr.dest(j+1)] = bm[j+1];
            pb[j+1] = icon + props.pr.dest(j+1);
            }

    if(new_rn > 4) 
        {
        Print((*this));
        Print(other);
        cout << "new_rn is " << new_rn << endl;
        Error("new_rn too big for this part!");
        }
    if(props.nsamen > 4) Error("nsamen too big for this part!");

    static Vector newdat;
    newdat.ReduceDimension(props.odimL*props.odimR);

    icon[1] = icon[2] = icon[3] = icon[4] = 1;
    inew[1] = inew[2] = inew[3] = inew[4] = 1;
    int basea = ind4(*pa[4],am[3],*pa[3],am[2],*pa[2],am[1],*pa[1]); 
    int baseb = ind4(*pb[4],bm[3],*pb[3],bm[2],*pb[2],bm[1],*pb[1]); 
    icon[1] = 2;
    int inca1 = ind4(*pa[4],am[3],*pa[3],am[2],*pa[2],am[1],*pa[1]) - basea; 
    int incb1 = ind4(*pb[4],bm[3],*pb[3],bm[2],*pb[2],bm[1],*pb[1]) - baseb; 
    icon[1] = 1;
    int inca2=0,incb2=0;
    if(props.nsamen == 2)
        {
        icon[2] = 2;
        inca2 = ind4(*pa[4],am[3],*pa[3],am[2],*pa[2],am[1],*pa[1]) - basea; 
        incb2 = ind4(*pb[4],bm[3],*pb[3],bm[2],*pb[2],bm[1],*pb[1]) - baseb; 
        icon[2] = 1;
        }

    Real *pv = p->v.Store()-1, *opv = other.p->v.Store()-1;
    for(inew[4] = 1; inew[4] <= mnew[4]; ++inew[4])
    for(inew[3] = 1; inew[3] <= mnew[3]; ++inew[3])
    for(inew[2] = 1; inew[2] <= mnew[2]; ++inew[2])
    for(inew[1] = 1; inew[1] <= mnew[1]; ++inew[1])
            {
            Real d = 0.0;
            if(props.nsamen == 1)
                {
                icon[1] = 1;
                int inda = ind4(*pa[4],am[3],*pa[3],am[2],*pa[2],am[1],*pa[1]);
                int indb = ind4(*pb[4],bm[3],*pb[3],bm[2],*pb[2],bm[1],*pb[1]);
                for(icon[1] = 1; icon[1] <= mcon[1]; icon[1]++, inda += inca1, indb += incb1)
                    d += pv[inda] * opv[indb];
                }
            else if(props.nsamen == 2)
                {
                icon[2] = icon[1] = 1;
                int inda = ind4(*pa[4],am[3],*pa[3],am[2],*pa[2],am[1],*pa[1]);
                int indb = ind4(*pb[4],bm[3],*pb[3],bm[2],*pb[2],bm[1],*pb[1]);
                int indaa = inda, indbb = indb;
                for(icon[2] = 1; icon[2] <= mcon[2]; icon[2]++, indaa += inca2, indbb += incb2)
                    {
                    inda = indaa; indb = indbb;
                    for(icon[1] = 1; icon[1] <= mcon[1]; ++icon[1], inda += inca1, indb += incb1)
                    d += pv[inda] * opv[indb];
                    }
                }
            else
                {
                for(icon[4] = 1; icon[4] <= mcon[4]; ++icon[4])
                for(icon[3] = 1; icon[3] <= mcon[3]; ++icon[3])
                for(icon[2] = 1; icon[2] <= mcon[2]; ++icon[2])
                for(icon[1] = 1; icon[1] <= mcon[1]; ++icon[1])
                    d +=      p->v(ind4(*pa[4],am[3],*pa[3],am[2],*pa[2],am[1],*pa[1])) 
                      * other.p->v(ind4(*pb[4],bm[3],*pb[3],bm[2],*pb[2],bm[1],*pb[1]));
                }

            newdat(ind4(inew[4],mnew[3],inew[3],mnew[2],inew[2],mnew[1],inew[1])) = d;
            }

    if(!p.unique()) allocate();

    p->v = newdat;


    } // directMultiply

#else

int
_ind(int rn, const array<int,NMAX>& n,
     int i1, int i2, int i3, int i4, 
     int i5, int i6, int i7, int i8)
    {
    switch(rn)
    {
    case 0:
        return (1);
    case 1:
        return (i1);
    case 2:
        return ((i2-1)*n[0]+i1);
    case 3:
        return (((i3-1)*n[1]+i2-1)*n[0]+i1);
    case 4:
        return ((((i4-1)*n[2]+i3-1)*n[1]+i2-1)*n[0]+i1);
    case 5:
        return (((((i5-1)*n[3]+i4-1)*n[2]+i3-1)*n[1]+i2-1)
                        *n[0]+i1);
    case 6:
        return ((((((i6-1)*n[4]+i5-1)*n[3]+i4-1)*n[2]+i3-1)
                        *n[1]+i2-1)*n[0]+i1);
    case 7:
        return (((((((i7-1)*n[5]+i6-1)*n[4]+i5-1)*n[3]+i4-1)
                        *n[2]+i3-1)*n[1]+i2-1)*n[0]+i1);
    case 8:
        return ((((((((i8-1)*n[6]+i7-1)*n[5]+i6-1)*n[4]+i5-1)
                        *n[3]+i4-1)*n[2]+i3-1)*n[1]+i2-1)*n[0]+i1);
    } //switch
    Error("_ind: Failed switch case");
    return 1;
    }

void ITensor::
directMultiply(const ITensor& other, 
               ProductProps& props, 
               IndexSet<Index>& new_index)
    {
    Counter u,  //uncontracted indices
            c;  //contracted indices

    int one = 1;
    array<int*,NMAX> li,
                     ri;
    li.assign(&one);
    ri.assign(&one);

    array<int,NMAX> nl,
                    nr;

    const int trn = this->is_.rn();
    const int orn = other.is_.rn();

    for(int j = 0; j < trn; ++j)
        {
        if(!props.contractedL[j+1])
            {
            ++u.rn; //(++u.r);
            u.n[u.rn] = is_[j].m();
            li[j] = &(u.i[u.rn]);
            new_index.addindexn(is_[j]);
            }
        else
            {
            for(int k = 0; k < orn; ++k)
                {
                if(is_[j] == other.is_[k])
                    {
                    ++c.rn; //(++c.r);
                    c.n[c.rn] = is_[j].m();
                    li[j] = &(c.i[c.rn]);
                    ri[k] = &(c.i[c.rn]);
                    break;
                    }
                }
            }
        nl[j] = is_[j].m();
        }

    for(int j = 0; j < orn; ++j)
        {
        if(!props.contractedR[j+1])
            {
            ++u.rn; //(++u.r);
            u.n[u.rn] = other.is_[j].m();
            ri[j] = &(u.i[u.rn]);
            new_index.addindexn(other.is_[j]);
            }
        nr[j] = other.is_[j].m();
        }

    Vector newdat(props.odimL*props.odimR);

    const Real* pL = this->p->v.Store()-1;
    const Real* pR = other.p->v.Store()-1;
    Real* pN = newdat.Store()-1;

    for(; u.notDone(); ++u)
        {
        Real& val = pN[u.ind];
        val = 0;
        for(c.reset(); c.notDone(); ++c)
            {
            val += pL[_ind(trn,nl,*li[0],*li[1],*li[2],*li[3],
                                  *li[4],*li[5],*li[6],*li[7])]
                 * pR[_ind(orn,nr,*ri[0],*ri[1],*ri[2],*ri[3],
                                  *ri[4],*ri[5],*ri[6],*ri[7])];
            }
        }

    if(!p.unique()) allocate();
    p->v = newdat;

    } //ITensor::directMultiply

#endif


ITensor& ITensor::
operator*=(const ITensor& other)
    {
    if(this == &other)
        {
        ITensor cp_oth(other);
        return operator*=(cp_oth);
        }

    if(this->isNull() || other.isNull())
        Error("Null ITensor in product");

    //Complex types are treated as just another index, of type ReIm
    //Multiplication is handled automatically with these simple tensor helpers
    if(hasindexn(Index::IndReIm()) && other.hasindexn(Index::IndReIm()) && 
	    !other.hasindexn(Index::IndReImP()) && !other.hasindex(Index::IndReImPP()) 
	    && !hasindex(Index::IndReImP()) && !hasindex(Index::IndReImPP()))
        {
        this->prime(ReIm);
        operator*=(ComplexProd() * primed(other,ReIm,2));
        return *this;
        }

    //These hold  regular new indices and the m==1 indices that appear in the result
    IndexSet<Index> new_index;

    array<const Index*,NMAX+1> new_index1_;
    int nr1_ = 0;

    //
    //Handle m==1 Indices
    //
    for(int k = is_.rn(); k < this->r(); ++k)
        {
        const Index& K = is_[k];
        if(!other.is_.hasindex1(K))
            new_index1_[++nr1_] = &K;
        }

    for(int j = other.is_.rn(); j < other.r(); ++j)
        {
        const Index& J = other.is_[j];
        if(!is_.hasindex1(J))
            new_index1_[++nr1_] = &J;
        }

    //
    //Special cases when one of the tensors
    //has only m==1 indices (effectively a scalar)
    //
    if(other.is_.rn() == 0)
        {
        scale_ *= other.scale_;
        scale_ *= other.p->v(1);
#ifdef DEBUG
        if((is_.rn()+nr1_) > NMAX) 
            {
            std::cout << "new r_ would be = " << is_.r() << "\n";
            std::cerr << "new r_ would be = " << is_.r() << "\n";
            Error("ITensor::operator*=: too many uncontracted indices in product (max is 8)");
            }
#endif
        for(int j = 1; j <= is_.rn(); ++j)
            new_index.addindexn(is_.index(j));
        //Keep current m!=1 indices, overwrite m==1 indices
        for(int j = 1; j <= nr1_; ++j) 
            new_index.addindex1( *(new_index1_[j]) );
        is_.swap(new_index);
        return *this;
        }
    else if(is_.rn() == 0)
        {
        scale_ *= other.scale_;
        scale_ *= p->v(1);
        p = other.p;
#ifdef DEBUG
        if((is_.rn()+nr1_) > NMAX) 
            {
            std::cout << "new r_ would be = " << is_.r() << "\n";
            std::cerr << "new r_ would be = " << is_.r() << "\n";
            Error("ITensor::operator*=: too many uncontracted indices in product (max is 8)");
            }
#endif
        for(int j = 1; j <= other.is_.rn(); ++j) 
            new_index.addindexn( other.is_.index(j) );
        for(int j = 1; j <= nr1_; ++j) 
            new_index.addindex1( *(new_index1_[j]) );
        is_.swap(new_index);
        return *this;
        }

    ProductProps props(*this,other);

#ifdef DEBUG
    if((is_.rn() + other.is_.rn() - 2*props.nsamen + nr1_) > NMAX) 
        {
        Print(*this);
        Print(other);
        Print(props.nsamen);
        cerr << "new m==1 indices\n";
        for(int j = 1; j <= nr1_; ++j) cerr << *(new_index1_.at(j)) << "\n";
        Error("ITensor::operator*=: too many uncontracted indices in product (max is 8)");
        }
#endif

    const
    bool do_matrix_multiply = (props.odimL*props.cdim*props.odimR) > 10000 ||
                              (is_.rn()+other.is_.rn()-2*props.nsamen) > 4 ||
                              is_.rn() > 4 ||
                              other.is_.rn() > 4;

    MatrixRefNoLink lref, rref;
    bool L_is_matrix,R_is_matrix;
    toMatrixProd(*this,other,props,lref,rref,
                 L_is_matrix,R_is_matrix,do_matrix_multiply);

    if(do_matrix_multiply || (L_is_matrix && R_is_matrix))
        {
        DO_IF_PS(++Prodstats::stats().c2;)

        //Do the matrix multiplication
        if(!p.unique()) allocate();

        p->v.ReDimension(rref.Nrows()*lref.Ncols());
        MatrixRef nref; 
        p->v.TreatAsMatrix(nref,rref.Nrows(),lref.Ncols());
        nref = rref*lref;

        //Handle m!=1 indices
        for(int j = 0; j < this->is_.rn(); ++j)
            if(!props.contractedL[j+1]) 
                new_index.addindexn( is_[j] );
        for(int j = 0; j < other.is_.rn(); ++j)
            if(!props.contractedR[j+1]) 
                new_index.addindexn( other.is_[j] );
        }
    else
        {
        directMultiply(other,props,new_index);
        }

    //Put in m==1 indices
    for(int j = 1; j <= nr1_; ++j) 
        new_index.addindex1( *(new_index1_.at(j)) );

    is_.swap(new_index);

    scale_ *= other.scale_;

    scaleOutNorm();

    return *this;
    } //ITensor::operator*=(ITensor)



ITensor& ITensor::
operator+=(const ITensor& other)
    {
    if(this == &other) 
        { 
        scale_ *= 2; 
        return *this; 
        }

    const bool complex_this = isComplex(*this);
    const bool complex_other = isComplex(other);
    if(!complex_this && complex_other)
        {
        return (*this = (*this * ITensor::Complex_1()) + other);
        }
    if(complex_this && !complex_other) 
        {
        return operator+=(other * ITensor::Complex_1());
        }

    if(fabs(is_.uniqueReal() - other.is_.uniqueReal()) > 1E-12)
        {
        cerr << format("this ur = %.10f, other.ur = %.10f\n")%is_.uniqueReal()%other.is_.uniqueReal();
        Print(*this);
        Print(other);
        Error("ITensor::operator+=: unique Reals don't match (different Index structure).");
        }

    if(this->scale_.sign() == 0)
        {
        *this = other;
        return *this;
        }

    if((other.scale_/scale_).isRealZero()) 
        { 
        return *this; 
        }

    if(this->p == other.p) 
        { 
        Real fac = 1+(other.scale_/scale_).real();
        scale_ *= fac;
        return *this;
        }

    solo();

    Vector& thisdat = p->v;
    const Vector& othrdat = other.p->v;

    Real scalefac = 1;
    if(scale_.magnitudeLessThan(other.scale_)) 
        {
        this->scaleTo(other.scale_); 
        }
    else
        {
        scalefac = (other.scale_/scale_).real();
        }

    bool same_ind_order = true;
    for(int j = 0; j < is_.rn(); ++j)
    if(is_[j] != other.is_[j])
        { 
        same_ind_order = false; 
        break; 
        }

    if(same_ind_order) 
        { 
        if(scalefac == 1)
            thisdat += othrdat; 
        else
            thisdat += scalefac*othrdat;
        return *this; 
        }

    Permutation P; 
    is_.getperm(other.is_,P);
    Counter c(other.is_);

    int *j[NMAX+1];
    for(int k = 1; k <= NMAX; ++k) j[P.dest(k)] = &(c.i[k]);
    //static int n[NMAX+1];
    int n[NMAX+1];
    for(int k = 1; k <= NMAX; ++k) 
        {
        n[P.dest(k)] = c.n[k];
        //n[k] = index_[k].m();
        }

//#ifdef STRONG_DEBUG
    //Real tot_this = thisdat.sumels();
    //Real tot_othr = othrdat.sumels();
//#endif

    if(scalefac == 1)
        {
        for(; c.notDone(); ++c)
            {
            thisdat((((((((*j[8]-1)*n[7]+*j[7]-1)*n[6]
            +*j[6]-1)*n[5]+*j[5]-1)*n[4]+*j[4]-1)*n[3]
            +*j[3]-1)*n[2]+*j[2]-1)*n[1]+*j[1])
            += othrdat(c.ind);
            }
        }
    else
        {
        for(; c.notDone(); ++c)
            {
            thisdat((((((((*j[8]-1)*n[7]+*j[7]-1)*n[6]
            +*j[6]-1)*n[5]+*j[5]-1)*n[4]+*j[4]-1)*n[3]
            +*j[3]-1)*n[2]+*j[2]-1)*n[1]+*j[1])
            += scalefac * othrdat(c.ind);
            }
        }


    /*
#ifdef STRONG_DEBUG
    Real new_tot = thisdat.sumels();
    Real compare = tot_this + scalefac*tot_othr;
    Real ref = Norm(thisdat);
    if(fabs(compare) > 1E-12 && fabs(new_tot-compare) > 1E-12 * ref)
	{
	Real di = new_tot - compare;
	cerr << format("new_tot = %f, compare = %f, dif = %f\n")%new_tot%compare%di;
	Error("Incorrect sum");
	}
#endif
    */

    return *this;
    } 

ITensor& ITensor::
operator-=(const ITensor& other)
    {
    if(this == &other) 
        { scale_ = 0; 
        return *this; 
        }
    scale_.negate();
    operator+=(other); 
    scale_.negate();
    return *this; 
    }

void ITensor::
fromMatrix11(const Index& i1, const Index& i2, const Matrix& M)
    {
    if(r() != 2) Error("fromMatrix11: incorrect rank");
    assert(hasindex(i1));
    assert(hasindex(i2));
    DO_IF_DEBUG(if(i1.m() != M.Nrows()) Error("fromMatrix11: wrong number of rows");)
    DO_IF_DEBUG(if(i2.m() != M.Ncols()) Error("fromMatrix11: wrong number of cols");)

    solo();
    scale_ = 1;

    MatrixRef dref; 
    if(i1 == is_[0])
        {
        p->v.TreatAsMatrix(dref,i2.m(),i1.m());
        dref = M.t();
        }
    else
        {
        p->v.TreatAsMatrix(dref,i1.m(),i2.m());
        dref = M;
        }
    }

void ITensor::
toMatrix11NoScale(const Index& i1, const Index& i2, Matrix& res) const
    {
    if(r() != 2) Error("toMatrix11: incorrect rank");
    assert(hasindex(i1));
    assert(hasindex(i2));
    res.ReDimension(i1.m(),i2.m());

    MatrixRef dref; 
    p->v.TreatAsMatrix(dref,is_[1].m(),is_[0].m());
    res = dref.t(i1==is_[0]); 
    }

void ITensor::
toMatrix11(const Index& i1, const Index& i2, Matrix& res) const
    { 
    toMatrix11NoScale(i1,i2,res); 
    res *= scale_.real(); 
    }

void ITensor::
toMatrix12NoScale(const Index& i1, const Index& i2, 
                  const Index& i3, Matrix& res) const
    {
    if(r() != 3) Error("toMatrix11: incorrect rank");
    assert(hasindex(i1));
    assert(hasindex(i2));
    assert(hasindex(i3));

    res.ReDimension(i1.m(),i2.m()*i3.m());

    const array<Index,NMAX> reshuf 
        = {{ i2, i3, i1,    Index::Null(), Index::Null(), 
             Index::Null(), Index::Null(), Index::Null() }};

    Permutation P; 
    is_.getperm(reshuf,P);

    Vector V;
    reshapeDat(P,V);
    res.TreatAsVector() = V;
    }

void ITensor::
toMatrix12(const Index& i1, const Index& i2, 
           const Index& i3, Matrix& res) const
    { 
    toMatrix12NoScale(i1,i2,i3,res); 
    res *= scale_.real(); 
    }

void ITensor::
fromMatrix12(const Index& i1, const Index& i2, 
             const Index& i3, const Matrix& M)
    {
    if(r() != 3) Error("fromMatrix12: incorrect rank");
    assert(hasindex(i1));
    assert(hasindex(i2));
    assert(hasindex(i3));

    ITensor Q(i3,i1,i2);
    Q.p->v = M.TreatAsVector();

    *this = Q;
    }

/*
// group i1,i2; i3,i4
void ITensor::toMatrix22(const Index& i1, const Index& i2, const Index& i3, const Index& i4,Matrix& res) const
{
    if(r() != 4) Error("toMatrix22: incorrect rank");
    assert(hasindex(i1));
    assert(hasindex(i2));
    assert(hasindex(i3));
    assert(hasindex(i4));
    int nrow = i1.m() * i2.m(), ncol = i3.m() * i4.m();
    if(nrow != res.Nrows()) Error("toMatrix22: wrong number of rows");
    if(ncol != res.Ncols()) Error("toMatrix22: wrong number of cols");
    res.ReDimension(nrow,ncol);
    const array<Index,NMAX+1> reshuf = {{ Index::Null(), i3, i4, i1, i2, Index::Null(), Index::Null(), Index::Null(), Index::Null() }};
    Permutation P; getperm(reshuf,P);
    Vector V; reshapeDat(P,V);
    res.TreatAsVector() = V;
    res *= scale_;
}

void ITensor::fromMatrix22(const Index& i1, const Index& i2, const Index& i3, const Index& i4,const Matrix& M)
{
    if(r() != 4) Error("fromMatrix22: incorrect rank");
    assert(hasindex(i1));
    assert(hasindex(i2));
    assert(hasindex(i3));
    assert(hasindex(i4));
    if(i3.m()*i4.m() != M.Ncols()) Error("fromMatrix22: wrong number of cols");
    if(i1.m()*i2.m() != M.Nrows()) Error("fromMatrix22: wrong number of rows");
    ITensor Q(i3,i4,i1,i2);
    Q.p->v = M.TreatAsVector();
    assignFrom(Q);
}



void ITensor::toMatrix21(const Index& i1, const Index& i2, const Index& i3, Matrix& res) const
{
    if(r() != 3) Error("toMatrix21: incorrect rank");
    assert(hasindex(i1));
    assert(hasindex(i2));
    res.ReDimension(i1.m()*i2.m(),i3.m());
    const array<Index,NMAX+1> reshuf = {{ Index::Null(), i3, i1, i2, Index::Null(), Index::Null(), Index::Null(), Index::Null(), Index::Null() }};
    Permutation P; getperm(reshuf,P);
    Vector V; reshapeDat(P,V);
    res.TreatAsVector() = V;
    res *= scale_;
}

void ITensor::toMatrix12(const Index& i1, const Index& i2, const Index& i3, Matrix& res) const
{
    if(r() != 3) Error("toMatrix12: incorrect rank");
    assert(hasindex(i1));
    assert(hasindex(i2));
    assert(hasindex(i3));
    res.ReDimension(i1.m(),i2.m()*i3.m());
    const array<Index,NMAX+1> reshuf = {{ Index::Null(), i2, i3, i1, Index::Null(), Index::Null(), Index::Null(), Index::Null(), Index::Null() }};
    Permutation P; getperm(reshuf,P);
    Vector V; reshapeDat(P,V);
    res.TreatAsVector() = V;
    res *= scale_;
}

void ITensor::fromMatrix21(const Index& i1, const Index& i2, const Index& i3, const Matrix& M)
{
    if(r() != 3) Error("fromMatrix21: incorrect rank");
    assert(hasindex(i1));
    assert(hasindex(i2));
    assert(hasindex(i3));
    if(i1.m()*i2.m() != M.Nrows()) Error("fromMatrix21: wrong number of rows");
    if(i3.m() != M.Ncols()) Error("fromMatrix21: wrong number of cols");
    ITensor Q(i3,i1,i2);
    Q.p->v = M.TreatAsVector();
    assignFrom(Q);
}

void ITensor::fromMatrix12(const Index& i1, const Index& i2, const Index& i3, const Matrix& M)
{
    if(r() != 3) Error("fromMatrix12: incorrect rank");
    assert(hasindex(i1) && hasindex(i2) && hasindex(i3));
    if(i1.m() != M.Nrows()) Error("fromMatrix12: wrong number of rows");
    if(i3.m()*i2.m() != M.Ncols()) Error("fromMatrix12: wrong number of cols");
    ITensor Q(i2,i3,i1);
    Q.p->v = M.TreatAsVector();
    assignFrom(Q);
}
*/

ostream& 
operator<<(ostream & s, const ITensor & t)
    {
    s << "ITensor r = " << t.r() << ": ";
    s << t.is_ << "\n";

    s << "  {log(scale)[incl in elems]=" << t.scale().logNum();


    if(t.isNull()) s << ", dat is null}\n";
    else 
        {
        s << ", L=" << t.vecSize();

        if(t.scale_.isFiniteReal())
            {
            Real nrm = t.norm();
            if(nrm >= 1E-2 && nrm < 1E5)
                s << format(", N=%.2f}\n") % nrm;
            else
                s << format(", N=%.1E}\n") % nrm;
            }
        else
            {
            s << ", N=too big} scale=" << t.scale() << "\n";
            }

        if(Global::printdat())
            {
            Real scale = 1.0;
            if(t.scale_.isFiniteReal()) scale = t.scale_.real();
            else s << "\n(omitting too large scale factor)" << endl;
            const Vector& v = t.p->v;
            Counter c(t.is_);
            for(; c.notDone(); ++c)
                {
                Real val = v(c.ind)*scale;
                if(fabs(val) > Global::printScale())
                    { s << "  " << c << (format(" %.10f\n") % val); }
                }
            }
        else 
            {
            s << "\n";
            }
        }
    return s;
    }



//
// ITDat
//

ITDat::
ITDat() 
    : 
    v(0)
    { }

ITDat::
ITDat(int size) 
    : 
    v(size)
    { 
    v = 0; 
    }

ITDat::
ITDat(const VectorRef& v_) 
    : 
    v(v_)
    { }

ITDat::
ITDat(Real r) 
    : 
    v(1)
    { 
    v = r; 
    }

ITDat::
ITDat(const ITDat& other) 
    : 
    v(other.v)
    { }

void ITDat:: 
read(std::istream& s) 
    { 
    int size = 0;
    s.read((char*) &size,sizeof(size));
    v.ReDimension(size);
    s.read((char*) v.Store(), sizeof(Real)*size);
    }


void ITDat::
write(std::ostream& s) const 
    { 
    const int size = v.Length();
    s.write((char*) &size, sizeof(size));
    s.write((char*) v.Store(), sizeof(Real)*size); 
    }

//
// commaInit
//

commaInit::
commaInit(ITensor& T,
          const Index& i1,
          const Index& i2,
          const Index& i3)
    : 
    T_(T),
    started_(false),
    c_(T.is_)
    { 
    if(T_.isNull()) 
        Error("Can't assign to null ITensor");

    boost::array<Index,NMAX> ii;
    ii.assign(Index::Null());

    if(i2 == Index::Null())
        {
        ii[0] = i1;
        }
    else
    if(i3 == Index::Null())
        {
        ii[0] = i2;
        ii[1] = i1;
        }
    else
        {
        ii[0] = i3;
        ii[1] = i2;
        ii[2] = i1;
        }
    try {
        T_.is_.getperm(ii,P_);
        }
    catch(const ITError& e)
        {
        Error("Not enough/wrong indices passed to commaInit");
        }

    T_.solo();
    T_.scaleTo(1);
    }

commaInit& commaInit::
operator<<(Real r)
    {
    started_ = true;
    return operator,(r);
    }


commaInit& commaInit::
operator,(Real r)
    {
    if(!started_)
        {
        Error("commaInit notation is T << #, #, #, ... ;");
        }
    if(c_.notDone()) 
        { T_.p->v(c_.ind) = r; ++c_; }
    else 
        { Error("Comma assignment list too long.\n"); }
    return *this;
    }

commaInit::
~commaInit()
    {
    T_.reshapeDat(P_);
    }

//
// Other methods defined in itensor.h
//

Real 
Dot(const ITensor& x, const ITensor& y)
    {
    ITensor res = x; 
    res *= y;
    if(res.r() != 0) 
        { 
        Print(x);
        Print(y);
        if(isComplex(x) || isComplex(y))
            {
            Error("Must use BraKet, not Dot, for complex ITensors");
            }
        Error("Bad Dot, product is not a scalar"); 
        }
    return res.toReal();
    }

void 
BraKet(const ITensor& x, const ITensor& y, Real& re, Real& im)
    {
    if(isComplex(x))
        {
        ITensor res = conj(x);
        res *= y;
        if(res.r() != 1) 
            {
            x.print("x");
            y.print("y");
            Error("Bad Dot, product not a complex scalar");
            }
        re = res(Index::IndReIm()(1));
        im = res(Index::IndReIm()(2));
        return;
        }
    else
    if(isComplex(y))
        {
        ITensor res = x;
        res *= y;
        if(res.r() != 1) 
            {
            x.print("x");
            y.print("y");
            Error("Bad Dot, product not a complex scalar");
            }
        re = res(Index::IndReIm()(1));
        im = res(Index::IndReIm()(2));
        return;
        }

    re = Dot(x,y);
    im = 0;
    }

