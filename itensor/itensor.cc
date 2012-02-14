#include "itensor.h"
using namespace std;
using boost::format;
using boost::array;
using boost::intrusive_ptr;

DatAllocator<ITDat> ITDat::allocator;

#ifdef DEBUG
#define ITENSOR_CHECK_NULL if(p == 0) Error("ITensor is null");
#else
#define ITENSOR_CHECK_NULL
#endif

Counter::
Counter() : rn_(0)
    {
    n.assign(1); n[0] = 0;
    reset(0);
    }

Counter::
Counter(const array<Index,NMAX+1>& ii,int rn,int r) 
    { init(ii,rn,r); }

Counter::
Counter(const IndexSet& is) { init(is); }

void Counter::
reset(int a)
    {
    i.assign(a);
    ind = 1;
    }

void Counter::
init(const array<Index,NMAX+1>& ii, int rn, int r)
    {
    rn_ = rn;
    r_ = r;
    n[0] = 0;
    for(int j = 1; j <= rn_; ++j) 
        { n[j] = ii[j].m(); }
    for(int j = rn_+1; j <= NMAX; ++j) 
        { n[j] = 1; }
    reset(1);
    }

void Counter::
init(const IndexSet& is)
    {
    rn_ = is.rn();
    r_ = is.r();
    n[0] = 0;
    for(int j = 1; j <= rn_; ++j) 
        { n[j] = is.index(j).m(); }
    for(int j = rn_+1; j <= NMAX; ++j) 
        { n[j] = 1; }
    reset(1);
    }

Counter& Counter::
operator++()
    {
    ++ind;
    ++i[1];
    if(i[1] > n[1])
    for(int j = 2; j <= rn_; ++j)
        {
        i[j-1] = 1;
        ++i[j];
        if(i[j] <= n[j]) break;
        }
    //set 'done' condition
    if(i[rn_] > n[rn_]) reset(0);
    return *this;
    }

bool Counter::
operator!=(const Counter& other) const
    {
    for(int j = 1; j <= NMAX; ++j)
        { if(i[j] != other.i[j]) return true; }
    return false;
    }

bool Counter::
operator==(const Counter& other) const
    { return !(*this != other); }

std::ostream&
operator<<(std::ostream& s, const Counter& c)
    {
    s << "("; 
    for(int i = 1; i < c.r_; ++i)
        {s << c.i[i] << " ";} 
    s << c.i[c.r_] << ")";
    return s;
    }

ostream& 
operator<<(ostream & s, const ITensor & t)
    {
    s << "log(scale)[incl in elems] = " << t.scale().logNum() 
      << ", r = " << t.r() << ": ";

    s << t.is_;

    if(t.isNull()) s << " (dat is null)\n";
    else 
        {
        if(t.scale_.isFiniteReal())
            {
            Real nrm = t.norm();
            if(nrm >= 1E-2)
                s << format(" (L=%d,N=%.2f)\n") % t.vecSize() % nrm;
            else
                s << format(" (L=%d,N=%.1E)\n") % t.vecSize() % nrm;
            }
        else
            {
            s << format(" (L=%d,N=too big)\n") 
                 % t.vecSize() << t.scale();
            }
        if(Globals::printdat())
            {
            Real scale = 1.0;
            if(t.scale_.isFiniteReal()) scale = t.scale_.real();
            else s << "\n(omitting too large scale factor)" << endl;
            const Vector& v = t.p->v;
            Counter c; t.initCounter(c);
            for(; c.notDone(); ++c)
                {
                Real val = v(c.ind)*scale;
                if(fabs(val) > Globals::printScale())
                    { s << c << " " << val << "\n"; }
                }
            }
        else 
            {
            s << "\n";
            }
        }
    return s;
    }

const Index& ITensor::
index(int j) const
    {
#ifdef DEBUG
    if(j > r()) 
        {
        Print(*this);
        Print(j);
        Error("j out of range (j > r_)");
        }
#endif
    return is_.index_[j];
    }

int ITensor::
m(int j) const
    {
#ifdef DEBUG
    if(j > r()) 
        {
        Print(*this);
        Print(j);
        Error("j out of range (j > r_)");
        }
#endif
    return is_.index_[j].m();
    }


//
// ITensor Constructors
//

ITensor::
ITensor()  
    : 
    p(0) 
    { }


ITensor::
ITensor(Real val) 
    { 
    allocate(1);
    p->v = val;
    is_.setUniqueReal();
    }

ITensor::
ITensor(const Index& i1) 
    :
    is_(i1)
	{ 
    allocate(i1.m());
    }

ITensor::
ITensor(const Index& i1, Real val) 
    :
    is_(i1)
	{ 
    allocate(i1.m());
    p->v = val; 
    }

ITensor::
ITensor(const Index& i1, const VectorRef& V) 
    : 
    p(new ITDat(V)),
    is_(i1)
	{ 
	if(i1.m() != V.Length()) 
	    Error("Mismatch of Index and Vector sizes.");
	}

ITensor::
ITensor(Index i1,Index i2) 
    :
    is_(i1,i2)
	{ 
    allocate(i1.m()*i2.m());
    }
    

ITensor::
ITensor(Index i1,Index i2,Real a) 
    :
    is_(i1,i2)
	{
    allocate(i1.m()*i2.m());
	if(is_.rn_ == 2) //then index order is i1, i2
	    {
	    const int nn = min(i1.m(),i2.m());
	    for(int i = 1; i <= nn; ++i) 
		p->v((i-1)*i1.m()+i) = a;
	    }
	else 
	    p->v(1) = a;
	}

ITensor::
ITensor(Index i1,Index i2,const MatrixRef& M) 
    :
    is_(i1,i2)
	{
    allocate(i1.m()*i2.m());
	if(i1.m() != M.Nrows() || i2.m() != M.Ncols()) 
	    Error("Mismatch of Index sizes and matrix.");
	MatrixRef dref; 
	p->v.TreatAsMatrix(dref,i2.m(),i1.m()); 
	dref = M.t();
	}

ITensor::
ITensor(Index i1, Index i2, Index i3,
            Index i4, Index i5, Index i6,
            Index i7, Index i8)
	{
	array<Index,NMAX> ii = {{ i1, i2, i3, i4, i5, i6, i7, i8 }};
	int size = 3;
	while(ii[size] != Index::Null()) ++size;
	int alloc_size; 
    is_ = IndexSet(ii,size,alloc_size);
	allocate(alloc_size);
	}

ITensor::
ITensor(const IndexVal& iv, Real fac) 
    :
    is_(iv.ind)
	{ 
    allocate(iv.ind.m());
	p->v(iv.i) = fac; 
	}

ITensor::
ITensor(const IndexVal& iv1, const IndexVal& iv2) 
    :
    is_(iv1.ind,iv2.ind)
	{ 
    allocate(iv1.ind.m()*iv2.ind.m());
	p->v((iv2.i-1)*iv1.ind.m()+iv1.i) = 1; 
	}

ITensor::
ITensor(const IndexVal& iv1, const IndexVal& iv2, 
        const IndexVal& iv3, const IndexVal& iv4, 
        const IndexVal& iv5, const IndexVal& iv6, 
        const IndexVal& iv7, const IndexVal& iv8)
	{
    //Construct ITensor
    array<Index,NMAX+1> ii = 
        {{ iv1.ind, iv2.ind, iv3.ind, iv4.ind, iv5.ind, 
           iv6.ind, iv7.ind, iv8.ind }};
    int size = 3; 
    while(size < NMAX && ii[size+1] != IndexVal::Null().ind) ++size;
    int alloc_size; 
    is_ = IndexSet(ii,size,alloc_size);
    allocate(alloc_size);

    //Assign specified element to 1
    array<int,NMAX+1> iv = 
        {{ iv1.i, iv2.i, iv3.i, iv4.i, iv5.i, iv6.i, iv7.i, iv8.i }};
    array<int,NMAX+1> ja; ja.assign(1);
    for(int k = 1; k <= rn(); ++k) //loop over indices of this ITensor
        for(int j = 0; j < size; ++j)  // loop over the given indices
    if(is_.index_[k] == ii[j]) 
        { ja[k] = iv[j]; break; }
    p->v(_ind(ja[1],ja[2],ja[3],ja[4],ja[5],ja[6],ja[7],ja[8])) = 1;
    }

ITensor::
ITensor(const std::vector<Index>& I) 
	{
    int alloc_size;
    is_ = IndexSet(I,I.size(),alloc_size);
	allocate(alloc_size);
	}

ITensor::
ITensor(const std::vector<Index>& I, const Vector& V) 
    : 
    p(new ITDat(V))
	{
    int alloc_size;
    is_ = IndexSet(I,I.size(),alloc_size);
	if(alloc_size != V.Length()) 
	    { Error("incompatible Index and Vector sizes"); }
	}


ITensor::
ITensor(const std::vector<Index>& I, const ITensor& other) 
    : 
    p(other.p), 
    scale_(other.scale_)
	{
    int alloc_size;
    is_ = IndexSet(I,I.size(),alloc_size);
	if(alloc_size != other.vecSize()) 
	    { Error("incompatible Index and ITensor sizes"); }
	}

ITensor::
ITensor(const std::vector<Index>& I, const ITensor& other, Permutation P) 
    : 
    p(0), 
    scale_(other.scale_)
    {
    int alloc_size;
    is_ = IndexSet(I,I.size(),alloc_size);
    if(alloc_size != other.vecSize()) 
        { Error("incompatible Index and ITensor sizes"); }
    if(P.is_trivial()) { p = other.p; }
    else               { allocate(); other.reshapeDat(P,p->v); }
    }

ITensor::
ITensor(ITmaker itm) 
	{
    is_ = IndexSet(Index::IndReIm());
    allocate(2);
    if(itm == makeComplex_1)  { p->v(1) = 1; }
    if(itm == makeComplex_i)  { p->v(2) = 1; }
    if(itm == makeConjTensor) { p->v(1) = 1; p->v(2) = -1; }
	}

void ITensor::
read(std::istream& s)
    { 
    bool isNull_;
    s.read((char*) &isNull_,sizeof(isNull_));
    if(isNull_) { *this = ITensor(); return; }

    is_.read(s);
    scale_.read(s);
    p = new ITDat(s);
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
val0() const 
	{ 
#ifdef DEBUG
    if(this->isNull())
        Error("ITensor is null");
#endif
    if(rn() != 0)
        {
        Print(*this);
        Error("ITensor is not a scalar");
        }

	try {
	    return p->v(1)*scale_.real(); 
	    }
	catch(TooBigForReal)
	    {
	    std::cout << "too big for real() in val0" << std::endl;
	    std::cerr << "too big for real() in val0" << std::endl;
	    std::cout << "p->v(1) is " << p->v(1) << std::endl;
	    std::cout << "scale is " << scale() << std::endl;
	    std::cout << "rethrowing" << std::endl;
	    throw;		// rethrow
	    }
	catch(TooSmallForReal)
	    {
	    std::cout << "warning: too small for real() in val0" << std::endl;
	    std::cerr << "warning: too small for real() in val0" << std::endl;
	    std::cout << "p->v(1) is " << p->v(1) << std::endl;
	    std::cout << "scale is " << scale() << std::endl;
	    return 0.0;
	    }
	return 0.0;
	}

Real ITensor::
val1(int i1) const
	{ 
#ifdef DEBUG
    if(this->isNull())
        Error("ITensor is null");
#endif
    if(rn() > 1)
        {
        Print(*this);
        Error("ITensor has rank > 1");
        }
    return p->v(i1)*scale_.real(); 
    }

Real& ITensor::
operator()()
	{ 
    if(rn() != 0)
        {
        std::cerr << format("# given = 0, rn_ = %d\n")%rn();
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
    if(rn() != 0)
        {
        std::cerr << format("# given = 0, rn_ = %d\n")%rn();
        Error("Not enough indices (requires all having m!=1)");
        }
    return scale_.real()*p->v(1);
    }

Real& ITensor::
operator()(const IndexVal& iv1)
	{
    if(rn() > 1) 
        {
        std::cerr << format("# given = 1, rn_ = %d\n")%rn();
        Error("Not enough m!=1 indices provided");
        }
    if(index(1) != iv1.ind)
        {
        Print(*this);
        Print(iv1);
        Error("Incorrect IndexVal argument to ITensor");
        }
    solo(); 
    scaleTo(1);
    return p->v(iv1.i);
	}

Real ITensor::
operator()(const IndexVal& iv1) const
	{
    ITENSOR_CHECK_NULL
    if(rn() > 1) 
        {
        std::cerr << format("# given = 1, rn() = %d\n")%rn();
        Error("Not enough m!=1 indices provided");
        }
    if(index(1) != iv1.ind)
        {
        Print(*this);
        Print(iv1);
        Error("Incorrect IndexVal argument to ITensor");
        }
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

void ITensor::
assignFrom(const ITensor& other)
    {
    if(this == &other) return;
    if(fabs(other.is_.uniqueReal() - is_.uniqueReal()) > 1E-12)
        {
        Print(*this); Print(other);
        Error("assignFrom: unique Real not the same"); 
        }
    Permutation P; 
    is_.getperm(other.is_,P);
    scale_ = other.scale_;
    if(p->count() != 1) { p = new ITDat(); }
#ifdef DO_ALT
    else { p->alt.clear(); }
#endif
    other.reshapeDat(P,p->v);
    }


void ITensor::
groupIndices(const array<Index,NMAX+1>& indices, int nind, 
             const Index& grouped, ITensor& res) const
    {
    array<bool,NMAX+1> isReplaced; isReplaced.assign(false);

    int tot_m = 1;
    int nn = 0; //number of m != 1 indices
    for(int j = 1; j <= nind; ++j) 
        {
        const Index& J = indices[j];
        if(J.m() != 1) ++nn;
        tot_m *= J.m();

        bool foundit = false;
        for(int k = 1; k <= r(); ++k) 
            { 
            if(index(k) == J) 
                {
                isReplaced[k] = true; 
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
    const int res_rn_ = rn() - nn + (nn == 0 ? 0 : 1);

    vector<Index> nindices; 
    nindices.reserve(r()-nind+1);
    Permutation P;
    int kk = 0, kr = 0;
    for(int j = 1; j <= rn(); ++j)
        {
        if(isReplaced[j])
            { 
            P.from_to(j,res_rn_+kr);
            kr += 1;
            }
        else 
            { 
            P.from_to(j,++kk);
            nindices.push_back(index(j)); 
            }
        }

    nindices.push_back(grouped);

    for(int j = rn()+1; j <= r(); ++j) 
        if(!isReplaced[j]) nindices.push_back(index(j));

    if(nn == 0) res = ITensor(nindices,*this);
    else        res = ITensor(nindices,*this,P); 
    }

void ITensor::
tieIndices(const array<Index,NMAX+1>& indices, int nind,
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
        const Index& K = is_.index_[k];
        for(int j = 1; j <= nind; ++j)
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
        for(int j = 1; j <= nind; ++j)
            cout << indices[j] << endl;
        Error("Couldn't find Index to tie");
        }

    //If tied indices have m==1, no work
    //to do; just replace indices
    if(tm == 1)
        {
        is_ = IndexSet(new_index_,new_r_,alloc_size,1);
        return;
        }

    int new_rn_ = rn()-nind+1;

    Counter nc(new_index_,new_rn_,new_r_);

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
    boost::intrusive_ptr<ITDat> np = new ITDat(alloc_size);
    Vector& resdat = np->v;

    const Vector& thisdat = p->v;
    for(; nc.notDone(); ++nc)
        {
        resdat(nc.ind) =
        thisdat(this->_ind(*ii[1],*ii[2],
                           *ii[3],*ii[4],
                           *ii[5],*ii[6],
                           *ii[7],*ii[8]));
        }

    is_.index_.swap(new_index_);
    is_.r_ = new_r_;
    is_.rn_ = new_rn_;
    is_.setUniqueReal();
    p.swap(np);

    } //ITensor::tieIndices

void ITensor::
tieIndices(const Index& i1, const Index& i2,
           const Index& tied)
    {
    array<Index,NMAX+1> inds =
        {{ Index::Null(), i1, i2, 
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
    array<Index,NMAX+1> inds =
        {{ Index::Null(), i1, i2, i3,
           Index::Null(), Index::Null(), 
           Index::Null(), Index::Null(), Index::Null() }};

    tieIndices(inds,3,tied);
    }

void ITensor::
tieIndices(const Index& i1, const Index& i2,
           const Index& i3, const Index& i4,
           const Index& tied)
    {
    array<Index,NMAX+1> inds =
        {{ Index::Null(), i1, i2, i3, i4,
           Index::Null(), Index::Null(), 
           Index::Null(), Index::Null() }};

    tieIndices(inds,4,tied);
    }

void ITensor::
trace(const array<Index,NMAX+1>& indices, int nind)
    {
    if(nind == 0) Error("No indices given");

    const int tm = indices[1].m();
    
    array<Index,NMAX+1> new_index_;

    //will count these up below
    int new_r_ = 0;
    int alloc_size = 1;

    array<bool,NMAX+1> traced;
    traced.assign(false);

    int nmatched = 0;
    for(int k = 1; k <= r(); ++k)
        {
        const Index& K = index(k);
        for(int j = 1; j <= nind; ++j)
        if(K == indices[j]) 
            { 
            if(indices[j].m() != tm)
                Error("Traced indices must have matching m's");
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
        for(int j = 1; j <= nind; ++j)
            cout << indices[j] << endl;
        Error("Couldn't find Index to trace");
        }

    //If traced indices have m==1, no work
    //to do; just replace indices
    if(tm == 1)
        {
        is_.r_ = new_r_;
        sortIndices(new_index_,is_.r_,is_.rn_,alloc_size,is_.index_,1);
        return;
        }

    int new_rn_ = rn()-nind;

    Counter nc(new_index_,new_rn_,new_r_);

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
    boost::intrusive_ptr<ITDat> np = new ITDat(alloc_size);
    Vector& resdat = np->v;

    const Vector& thisdat = p->v;
    for(; nc.notDone(); ++nc)
        {
        Real newval = 0;
        for(trace_ind = 1; trace_ind <= tm; ++trace_ind)
            {
            newval += 
            thisdat(this->_ind(*ii[1],*ii[2],
                               *ii[3],*ii[4],
                               *ii[5],*ii[6],
                               *ii[7],*ii[8]));
            }
        resdat(nc.ind) = newval;
        }

    is_.index_.swap(new_index_);
    is_.setUniqueReal();
    is_.r_ = new_r_;
    is_.rn_ = new_rn_;
    p.swap(np);

    } //ITensor::trace

void ITensor::
trace(const Index& i1, const Index& i2)
    {
    array<Index,NMAX+1> inds =
        {{ Index::Null(), i1, i2, 
           Index::Null(), Index::Null(), 
           Index::Null(), Index::Null(), 
           Index::Null(), Index::Null() }};

    trace(inds,2);
    }

void ITensor::
trace(const Index& i1, const Index& i2, const Index& i3)
    {
    array<Index,NMAX+1> inds =
        {{ Index::Null(), i1, i2, i3,
           Index::Null(), Index::Null(), 
           Index::Null(), Index::Null(), Index::Null() }};

    trace(inds,3);
    }

void ITensor::
trace(const Index& i1, const Index& i2,
      const Index& i3, const Index& i4)
    {
    array<Index,NMAX+1> inds =
        {{ Index::Null(), i1, i2, i3, i4,
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
    indices.reserve(r());
    int w = -1;
    for(int j = 1; j <= r(); ++j)
        {
        if(index(j) == small)
            {
            w = j;
            indices.push_back(big);
            }
        else 
            {
            indices.push_back(index(j));
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
    w = res.findindex(big);

    array<int,NMAX+1> inc;
    //Make sure all other inc's are zero
    inc.assign(0);
    inc.at(w) = start;

    Counter c; 
    initCounter(c);

    const Vector& thisdat = p->v;
    Vector& resdat = res.p->v;
    for(; c.notDone(); ++c)
        {
        resdat(res._ind(c.i[1]+inc[1],c.i[2]+inc[2],
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
    return (p == 0 ? 0 : p->v.Length()); 
    }

int ITensor::
maxSize() const 
    { 
    int ms = 1;
    for(int j = 1; j <= rn(); ++j)
        ms *= m(j);
    return ms;
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
    if(p->count() != 1) 
        { 
        boost::intrusive_ptr<ITDat> new_p(new ITDat(v)); 
        p.swap(new_p); 
        }
    else
        {
        p->v = v;
#ifdef DO_ALT
        p->alt.clear();
#endif
        }
	}

void ITensor::
reshape(const Permutation& P)
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

void ITensor::
Randomize() 
    { 
    solo(); 
    p->v.Randomize(); 
    }

void ITensor::
SplitReIm(ITensor& re, ITensor& im) const
	{
	re = *this; im = *this;
	if(!isComplex()) { im *= 0; return; }
	//re *= Index::IndReIm()(1); im *= Index::IndReIm()(2);

	re.mapindex(Index::IndReIm(),Index::IndReImP());
	im.mapindex(Index::IndReIm(),Index::IndReImP());
	re *= Index::IndReImP()(1);
	im *= Index::IndReImP()(2);
	}

Real ITensor::
sumels() const 
    { return p->v.sumels() * scale_.real(); }

Real ITensor::
norm() const 
    { return fabs(Norm(p->v) * scale_.real()); }

void ITensor::
pseudoInvertElems(Real cutoff)
    {
    solo();
    //Invert scale_
    scale_ = scale_.pow(-1);

    //Invert elems
    for(int j = 1; j <= p->v.Length(); ++j)
        {
        Real elem = p->v(j);
        p->v(j) = (fabs(elem) <= cutoff ? 0 : 1./elem);
        }
    }

void ITensor::
scaleOutNorm() const
	{
    Real f = Norm(p->v);
    if(fabs(f-1) < 1E-12) return;
    solo();
    if(f != 0) { p->v *= 1.0/f; scale_ *= f; }
	}

void ITensor::
scaleTo(LogNumber newscale) const
	{
    if(scale_ == newscale) return;
    solo();
    if(newscale.isRealZero()) 
        { p->v *= 0; }
    else 
        {
        scale_ /= newscale;
        if(scale_.isRealZero()) p->v *= 0;
        else p->v *= scale_.real();
        }
    scale_ = newscale;
	}

void ITensor::
print(std::string name,Printdat pdat) const 
    { 
    Globals::printdat() = (pdat==ShowData); 
    std::cerr << "\n" << name << " =\n" << *this << "\n"; 
    Globals::printdat() = false; 
    }

void ITensor::
initCounter(Counter& C) const 
    { 
    C.init(is_);
    }

void ITensor::
allocate(int dim) 
    { 
    p = new ITDat(dim); 
    }

void ITensor::
allocate() 
    { p = new ITDat(); }

#ifdef DO_ALT
void 
newAltDat(const Permutation& P) const
    { p->alt.push_back(PDat(P)); }

PDat& ITensor::
lastAlt() const 
    { return p->alt.back(); } 

#endif //DO_ALT

void ITensor::
solo() const
	{
    ITENSOR_CHECK_NULL
    if(p->count() != 1) 
        p = new ITDat(*p);
	}

int ITensor::
_ind(int i1, int i2, int i3, int i4, 
     int i5, int i6, int i7, int i8) const
    {
    ITENSOR_CHECK_NULL
    switch(rn())
    {
    case 0:
        return (1);
        break;
    case 1:
        return (i1);
        break;
    case 2:
        return ((i2-1)*m(1)+i1);
        break;
    case 3:
        return (((i3-1)*m(2)+i2-1)*m(1)+i1);
        break;
    case 4:
        return ((((i4-1)*m(3)+i3-1)*m(2)+i2-1)*m(1)+i1);
        break;
    case 5:
        return (((((i5-1)*m(4)+i4-1)*m(3)+i3-1)*m(2)+i2-1)
                        *m(1)+i1);
        break;
    case 6:
        return ((((((i6-1)*m(5)+i5-1)*m(4)+i4-1)*m(3)+i3-1)
                        *m(2)+i2-1)*m(1)+i1);
        break;
    case 7:
        return (((((((i7-1)*m(6)+i6-1)*m(5)+i5-1)*m(4)+i4-1)
                        *m(3)+i3-1)*m(2)+i2-1)*m(1)+i1);
        break;
    case 8:
        return ((((((((i8-1)*m(7)+i7-1)*m(6)+i6-1)*m(5)+i5-1)
                        *m(4)+i4-1)*m(3)+i3-1)*m(2)+i2-1)*m(1)+i1);
        break;
    } //switch(rn_)
    Error("ITensor::_ind: Failed switch case");
    return 1;
    }

int ITensor::
_ind2(const IndexVal& iv1, const IndexVal& iv2) const
    {
    if(rn() > 2) 
        {
        std::cerr << format("# given = 2, rn_ = %d\n")%rn();
        Error("Not enough m!=1 indices provided");
        }
    if(index(1) == iv1.ind && index(2) == iv2.ind)
        return ((iv2.i-1)*m(1)+iv1.i);
    else if(index(1) == iv2.ind && index(2) == iv1.ind)
        return ((iv1.i-1)*m(1)+iv2.i);
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
    array<const IndexVal*,NMAX+1> iv = 
        {{ 0, &iv1, &iv2, &iv3, &iv4, &iv5, &iv6, &iv7, &iv8 }};
    array<int,NMAX+1> ja; ja.assign(1);
    //Loop over the given IndexVals
    int j = 1, nn = 0;
    while(iv[j]->ind != Index::Null())
        {
        //Loop over indices of this ITensor
        bool matched = false;
        for(int k = 1; k <= r(); ++k)
            {
            if(index(k) == iv[j]->ind)
                {
                matched = true;
                if(k <= rn()) ++nn;
                ja[k] = iv[j]->i;
                break;
                }
            }
        if(!matched)
            {
            Print(*this);
            Print(*iv[j]);
            Error("Extra/incorrect IndexVal argument to ITensor");
            }
        ++j;
        }

    if(nn != rn())
        {
        Error("Too few m!=1 indices provided");
        }

    return _ind(ja[1],ja[2],ja[3],ja[4],ja[5],ja[6],ja[7],ja[8]);
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
    Counter c; 
    initCounter(c);
    array<int,NMAX+1> n;
    for(int j = 1; j <= c.rn_; ++j) n[ind[j]] = c.n[j];

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

    if(rn() == 2 && ind[1] == 2 && ind[2] == 1)
        {
        MatrixRef xref; 
        thisdat.TreatAsMatrix(xref,c.n[2],c.n[1]);
        rdat = Matrix(xref.t()).TreatAsVector();
        return; 
        }
    else if(rn() == 3)
        {
        DO_IF_PS(int idx = ((ind[1]-1)*3+ind[2]-1)*3+ind[3]; prodstats.perms_of_3[idx] += 1; )
        //Arranged loosely in order of frequency of occurrence
        Bif3(2,1,3) Loop3(i2,i1,i3)
        Bif3(2,3,1) Loop3(i2,i3,i1) //cyclic
        Bif3(3,1,2) Loop3(i3,i1,i2)
        //Bif3(1,3,2) Loop3(i1,i3,i2)
        //Bif3(3,2,1) Loop3(i3,i2,i1)
        }
    else if(rn() == 4)
        {
        DO_IF_PS(int idx = (((ind[1]-1)*4+ind[2]-1)*4+ind[3]-1)*4+ind[4]; prodstats.perms_of_4[idx] += 1; )
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
    else if(rn() == 5)
        {
        DO_IF_PS(int idx = ((((ind[1]-1)*5+ind[2]-1)*5+ind[3]-1)*5+ind[4]-1)*5+ind[5]; prodstats.perms_of_5[idx] += 1; )
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
    else if(rn() == 6)
        {
        DO_IF_PS(int idx = (((((ind[1]-1)*6+ind[2]-1)*6+ind[3]-1)*6+ind[4]-1)*6+ind[5]-1)*6+ind[6]; prodstats.perms_of_6[idx] += 1; )
        //Arranged loosely in order of frequency of occurrence
        Bif6(2,4,1,3,5,6) Loop6(i2,i4,i1,i3,i5,i6)
        Bif6(1,4,2,3,5,6) Loop6(i1,i4,i2,i3,i5,i6)
        Bif6(2,4,1,5,3,6) Loop6(i2,i4,i1,i5,i3,i6)
        Bif6(1,2,4,5,3,6) Loop6(i1,i2,i4,i5,i3,i6)
        Bif6(3,4,1,5,6,2) Loop6(i3,i4,i1,i5,i6,i2)
        }
    DO_IF_PS(prodstats.c4 += 1;)

    //The j's are pointers to the i's of xdat's Counter,
    //but reordered in a way appropriate for rdat
    array<int*,NMAX+1> j;
    for(int k = 1; k <= NMAX; ++k) { j[ind[k]] = &(c.i[k]); }

    //Catch-all loop that works for any tensor
    for(; c.notDone(); ++c)
        {
        rdat((((((((*j[8]-1)*n[7]+*j[7]-1)*n[6]+*j[6]-1)*n[5]+*j[5]-1)*n[4]+*j[4]-1)*n[3]+*j[3]-1)*n[2]+*j[2]-1)*n[1]+*j[1])
            = thisdat(c.ind);
        }
    } // ITensor::reshapeDat

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

    };

ProductProps::
ProductProps(const ITensor& L, const ITensor& R) 
    :
    nsamen(0), 
    cdim(1), 
    odimL(-1), 
    odimR(-1),
    lcstart(-1), 
    rcstart(-1)
    {
    for(int j = 1; j <= NMAX; ++j) 
        contractedL[j] = contractedR[j] = false;

    for(int j = 1; j <= L.rn(); ++j)
	for(int k = 1; k <= R.rn(); ++k)
	    if(L.index(j) == R.index(k))
		{
		if(lcstart == -1) lcstart = j;
		if(rcstart == -1) rcstart = k;

		++nsamen;
		pl.from_to(j,nsamen);
		pr.from_to(k,nsamen);

		contractedL[j] = contractedR[k] = true;

        cdim *= L.index(j).m();
		}

    odimL = L.p->v.Length()/cdim;
    odimR = R.p->v.Length()/cdim;
    }

//Converts ITensor dats into MatrixRef's that can be multiplied as rref*lref
//contractedL/R[j] == true if L/R.indexn(j) contracted
void 
toMatrixProd(const ITensor& L, const ITensor& R, ProductProps& props,
             MatrixRefNoLink& lref, MatrixRefNoLink& rref)
    {
    assert(L.p != 0);
    assert(R.p != 0);
    const Vector &Ldat = L.p->v, &Rdat = R.p->v;

    bool L_is_matrix = true, R_is_matrix = true;
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
        if(!(props.contractedL[1] || props.contractedL[L.rn()])) 
            {
            L_is_matrix = false; 
            }
        if(!(props.contractedR[1] || props.contractedR[R.rn()]))
            {
            R_is_matrix = false; 
            }
        }

    /*
    if(L_is_matrix)
        {
        cerr << "L is matrix, props.pl = " << props.pl << "\n";
        cerr << format("contractedL = %d %d %d %d %d\n")
                % props.contractedL[1]
                % props.contractedL[2]
                % props.contractedL[3]
                % props.contractedL[4]
                % props.contractedL[5];
        }
    else
        {
        cerr << "L is not matrix\n";
        }
    if(R_is_matrix)
        {
        cerr << "R is matrix, props.pl = " << props.pr << "\n";
        cerr << format("contractedR = %d %d %d %d %d\n")
                % props.contractedR[1]
                % props.contractedR[2]
                % props.contractedR[3]
                % props.contractedR[4]
                % props.contractedR[5];
        }
    else
        {
        cerr << "R is not matrix\n";
        }
        */

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
    else
        {
        bool done_with_L = false;
#ifdef DO_ALT
        //Not matrix, see if alternate dat is
        for(std::vector<PDat>::const_iterator it = L.p->alt.begin();
            it != L.p->alt.end();
            ++it)
            {
            const PDat& Alt = *it;
            bool front_matrix=true;
            for(int j = 1; j <= props.nsamen; ++j)
            if(!GET(props.contractedL,Alt.I.dest(j)))
                { front_matrix = false; break; }

            if(front_matrix) 
                { 
                Alt.v.TreatAsMatrix(lref,props.odimL,props.cdim); lref.ApplyTrans();
                done_with_L = true;
                L_is_matrix = true;
                break;
                }

            bool back_matrix=true;
            for(int j = L.rn(); j > (L.rn()-props.nsamen); --j)
            if(!GET(props.contractedL,Alt.I.dest(j)))
                { back_matrix = false; break; }

            if(back_matrix)
                {
                Alt.v.TreatAsMatrix(lref,props.cdim,props.odimL); 
                done_with_L = true;
                L_is_matrix = true;
                        break;
                }
            } //for int n
#endif
        //Finish making the permutation (stick non contracted inds on the back)
        if(!done_with_L)
            {
            int q = props.nsamen;
            for(int j = 1; j <= L.rn(); ++j)
            if(!props.contractedL[j]) props.pl.from_to(j,++q);
            if(L_is_matrix) Error("Calling reshapeDat although L is matrix.");
#ifdef DO_ALT
            L.newAltDat(props.pl);
            L.reshapeDat(props.pl,L.lastAlt().v);
            L.lastAlt().v.TreatAsMatrix(lref,props.odimL,props.cdim); lref.ApplyTrans();
#else
            Vector lv; L.reshapeDat(props.pl,lv);
            lv.TreatAsMatrix(lref,props.odimL,props.cdim); lref.ApplyTrans();
#endif
            done_with_L = true;
            }
        assert(done_with_L);
        }

    if(R_is_matrix) 
        {
        if(props.contractedR[1]) 
            { Rdat.TreatAsMatrix(rref,props.odimR,props.cdim); }
        else                    
            { Rdat.TreatAsMatrix(rref,props.cdim,props.odimR); rref.ApplyTrans(); }
        }
    else
        {
        bool done_with_R = false;
#ifdef DO_ALT
        //Not matrix, see if alternate dat is
        for(std::vector<PDat>::const_iterator it = R.p->alt.begin();
            it != R.p->alt.end();
            ++it)
            {
            const PDat& Alt = *it;
            bool front_matrix=true;
            for(int j = 1; j <= props.nsamen; ++j)
            if(!GET(props.contractedR,Alt.I.dest(j)))
                { front_matrix = false; break; }

            if(front_matrix) 
                { 
                Alt.v.TreatAsMatrix(rref,props.odimR,props.cdim); 
                done_with_R = true;
                R_is_matrix = true;
                break;
                }

            bool back_matrix=true;
            for(int j = R.rn(); j > (R.rn()-props.nsamen); --j)
            if(!GET(props.contractedR,Alt.I.dest(j)))
                { back_matrix = false; break; }

            if(back_matrix)
                {
                Alt.v.TreatAsMatrix(rref,props.cdim,props.odimR); rref.ApplyTrans();
                done_with_R = true;
                R_is_matrix = true;
                        break;
                }
            } //for int n
#endif
        //Finish making the permutation (stick non contracted inds on the back)
        if(!done_with_R)
            {
            int q = props.nsamen;
            for(int j = 1; j <= R.rn(); ++j)
            if(!props.contractedR[j]) props.pr.from_to(j,++q);
            if(R_is_matrix) Error("Calling reshape even though R is matrix.");
#ifdef DO_ALT
            R.newAltDat(props.pr);
            R.reshapeDat(props.pr,R.lastAlt().v);
            R.lastAlt().v.TreatAsMatrix(rref,props.odimR,props.cdim);
#else
            Vector rv; R.reshapeDat(props.pr,rv);
            rv.TreatAsMatrix(rref,props.odimR,props.cdim);
#endif
            done_with_R = true;
            }
        }

#ifdef COLLECT_PRODSTATS
    if(L.rn_ > R.rn_) 
        {
        ++prodstats.global[std::make_pair(L.rn_,R.rn_)];
        }
    else 
        {
        ++prodstats.global[std::make_pair(R.rn_,L.rn_)];
        }
    ++prodstats.total;
    if(L_is_matrix) ++prodstats.did_matrix;
    if(R_is_matrix) ++prodstats.did_matrix;
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

    //These hold the indices from other 
    //that will be added to this->index_
    int nr1_ = 0;
    static array<const Index*,NMAX+1> extra_index1_;

    //------------------------------------------------------------------
    //Handle m==1 Indices: set union
    for(int j = other.rn()+1; j <= other.r(); ++j)
        {
        const Index& J = other.is_.index_[j];
        bool this_has_index = false;
        for(int k = this->rn()+1; k <= this->r(); ++k)
            { 
            if(is_.index_[k] == J) 
                { 
                this_has_index = true; 
                break; 
                } 
            }
        if(!this_has_index) extra_index1_[++nr1_] = &J;
        }

    static array<Index,NMAX+1> new_index_;

    if(other.rn() == 0)
        {
        scale_ *= other.scale_;
        scale_ *= other.p->v(1);
        for(int j = 1; j <= nr1_; ++j) 
            { 
            is_.index_[is_.r_+j] = *(extra_index1_[j]); 
            }
        is_.r_ += nr1_;
        is_.setUniqueReal();
        return *this;
        }
    else if(rn() == 0)
        {
        scale_ *= other.scale_;
        scale_ *= p->v(1);
        p = other.p;
        is_.rn_ = other.is_.rn_;
        //Copy other's m!=1 indices
        for(int j = 1; j <= rn(); ++j) 
            {
            //cerr << "Copying [" << j << "] " << other.index_[j] << "\n";
            new_index_[j] = other.index(j);
            }
        //Move current m==1 indices past rn_
        for(int j = 1; j <= r(); ++j) 
            {
            //cerr << "Copying [" << rn_+j << "] " << index_[j] << "\n";
            new_index_[rn()+j] = index(j);
            }
        is_.r_ += is_.rn_;
        //Get the extra m==1 indices from other
        for(int j = 1; j <= nr1_; ++j) 
            {
            //cerr << "Copying [" << r_+j << "] " << *(extra_index1_[j]) << "\n";
            new_index_[r()+j] = *(extra_index1_[j]);
            }
        is_.r_ += nr1_;
        is_.index_.swap(new_index_);
        is_.setUniqueReal();
        return *this;
        }

    ProductProps props(*this,other);
    MatrixRefNoLink lref, rref;
    toMatrixProd(*this,other,props,lref,rref);

    if(p->count() != 1) { p = new ITDat(); }
#ifdef DO_ALT
    else { p->alt.clear(); }
#endif
    Vector& thisdat = p->v; 
    
    const int ni = lref.Ncols(), nj = lref.Nrows(), nk = rref.Nrows();
    Matrix L(lref), R(rref);
    thisdat.ReDimension(ni*nj*nk);
    
    for(int j = 1; j <= nj; ++j) 
    for(int k = 1; k <= nk; ++k) 
    for(int i = 1; i <= ni; ++i)
        { thisdat(((j-1)*nk+k-1)*ni+i) =  R(k,j) * L(j,i); }

    if((is_.r_ + other.is_.rn_ - props.nsamen + nr1_) > NMAX) 
        Error("ITensor::operator/=: too many indices in product.");

    //Handle m!=1 indices
    int nrn_ = 0;
    for(int j = 1; j <= is_.rn_; ++j)
        { if(!props.contractedL[j]) new_index_[++nrn_] = this->index(j); }
    for(int j = 1; j <= other.is_.rn_; ++j)
        { if(!props.contractedR[j]) new_index_[++nrn_] = other.index(j); }
    for(int j = 1; j <= is_.rn_; ++j)
        { if(props.contractedL[j])  new_index_[++nrn_] = this->index(j); }

    for(int j = rn()+1; j <= r(); ++j) new_index_[nrn_+j-is_.rn_] = index(j);
    is_.r_ = (r()-rn()) + nrn_;
    for(int j = 1; j <= nr1_; ++j) new_index_[r()+j] = *(extra_index1_[j]);
    is_.r_ += nr1_;

    is_.rn_ = nrn_;

    is_.index_.swap(new_index_);
    is_.setUniqueReal();
    
    scale_ *= other.scale_;

    scaleOutNorm();

    return *this;
    }


inline int 
ind4(int i4, int m3, int i3, int m2, int i2, int m1, int i1)
    {
    return (((i4-1)*m3+i3-1)*m2+i2-1)*m1+i1;
    }

//#define NEW_DIRECT_MULT

#ifndef NEW_DIRECT_MULT

void ITensor::
directMultiply(const ITensor& other, ProductProps& props, 
               int& new_rn_, array<Index,NMAX+1>& new_index_)
    {
    int am[NMAX+1], bm[NMAX+1], mcon[NMAX+1], mnew[NMAX+1];
    int *pa[NMAX+1], *pb[NMAX+1];
    int one = 1;
    for(int j = 1; j <= NMAX; ++j)
        {
        //Set *pa[j],*pb[j],mcon[j] and mnew[j]
        // to 1 unless set otherwise below
        pa[j] = pb[j] = &one, mcon[j] = mnew[j] = 1;
        am[j] = m(j);
        bm[j] = other.m(j);
        }
    int icon[NMAX+1], inew[NMAX+1];
    for(int j = 1; j <= this->rn(); ++j)
        if(!props.contractedL[j]) 
            {
            new_index_[++new_rn_] = index(j);
            mnew[new_rn_] = am[j];
            pa[j] = inew + new_rn_;
            }
        else
            {
            mcon[props.pl.dest(j)] = am[j];
            pa[j] = icon + props.pl.dest(j);
            }

    for(int j = 1; j <= other.rn(); ++j)
        if(!props.contractedR[j]) 
            {
            new_index_[++new_rn_] = other.index(j);
            mnew[new_rn_] = bm[j];
            pb[j] = inew + new_rn_;
            }
        else
            {
            mcon[props.pr.dest(j)] = bm[j];
            pb[j] = icon + props.pr.dest(j);
            }

    if(new_rn_ > 4) 
        {
        Print((*this));
        Print(other);
        cout << "new_rn_ is " << new_rn_ << endl;
        Error("new_rn_ too big for this part!");
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

    if(p->count() != 1) { p = new ITDat(); } 
    p->v = newdat;

    DO_IF_PS(++prodstats.c1;)

    } // directMultiply

#else

void ITensor::
directMultiply(const ITensor& other, ProductProps& props, 
               int& new_rn_, array<Index,NMAX+1>& new_index_)
    {
    //will count these up below
    int new_r_ = 0;
    int alloc_size = 1;

    array<bool,NMAX+1> traced;
    traced.assign(false);

    int nmatched = 0;
    for(int k = 1; k <= r(); ++k)
        {
        const Index& K = index(k);
        for(int j = 1; j <= nind; ++j)
        if(K == indices[j]) 
            { 
            if(indices[j].m() != tm)
                Error("Traced indices must have matching m's");
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
        for(int j = 1; j <= nind; ++j)
            cout << indices[j] << endl;
        Error("Couldn't find Index to trace");
        }

    //If traced indices have m==1, no work
    //to do; just replace indices
    if(tm == 1)
        {
        is_.r_ = new_r_;
        sortIndices(new_index_,is_.r_,is_.rn_,alloc_size,is_.index_,1);
        return;
        }

    int new_rn_ = is_.rn_-nind;

    Counter nc(new_index_,new_rn_,new_r_);

    //Set up ii pointers to link
    //elements of res to appropriate
    //elements of *this
    int trace_ind = 0;
    array<int*,NMAX+1> ii;
    int n = 1;
    for(int j = 1; j <= is_.r_; ++j)
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
    boost::intrusive_ptr<ITDat> np = new ITDat(alloc_size);
    Vector& resdat = np->v;

    const Vector& thisdat = p->v;
    for(; nc.notDone(); ++nc)
        {
        Real newval = 0;
        for(trace_ind = 1; trace_ind <= tm; ++trace_ind)
            {
            newval += 
            thisdat(this->_ind(*ii[1],*ii[2],
                               *ii[3],*ii[4],
                               *ii[5],*ii[6],
                               *ii[7],*ii[8]));
            }
        resdat(nc.ind) = newval;
        }

    is_.index_.swap(new_index_);
    is_.setUniqueReal();
    is_.r_ = new_r_;
    is_.rn_ = new_rn_;

    p.swap(np);

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

    //Complex types are treated as just another index, of type ReIm
    //Multiplication is handled automatically with these simple tensor helpers
    if(findindexn(Index::IndReIm()) && other.findindexn(Index::IndReIm()) && 
	    !other.findindexn(Index::IndReImP()) && !other.hasindex(Index::IndReImPP()) 
	    && !hasindex(Index::IndReImP()) && !hasindex(Index::IndReImPP()))
        {
        static ITensor primer(Index::IndReIm(),Index::IndReImP(),1.0);
        static ITensor primerP(Index::IndReIm(),Index::IndReImPP(),1.0);
        static ITensor prod(Index::IndReIm(),Index::IndReImP(),Index::IndReImPP());
        static bool first = true;
        if(first)
	    {
            IndexVal iv0(Index::IndReIm(),1), iv1(Index::IndReImP(),1), iv2(Index::IndReImPP(),1);
            iv0.i = 1; iv1.i = 1; iv2.i = 1; prod(iv0,iv1,iv2) = 1.0;
            iv0.i = 1; iv1.i = 2; iv2.i = 2; prod(iv0,iv1,iv2) = -1.0;
            iv0.i = 2; iv1.i = 2; iv2.i = 1; prod(iv0,iv1,iv2) = 1.0;
            iv0.i = 2; iv1.i = 1; iv2.i = 2; prod(iv0,iv1,iv2) = 1.0;
            first = false;
	    }
        operator*=(primer);
        operator*=(prod * (other * primerP));
        return *this;
        }

    //These hold  regular new indices and the m==1 indices that appear in the result
    static array<Index,NMAX+1> new_index_;
    static array<const Index*,NMAX+1> new_index1_;
    int nr1_ = 0;

    //
    //Handle m==1 Indices
    //
    for(int k = rn()+1; k <= this->r(); ++k)
        {
        const Index& K = is_.index_[k];
        for(int j = other.rn()+1; j <= other.r(); ++j)
	    if(other.index(j) == K) 
		goto skip_this;
        new_index1_[++nr1_] = &K;
        skip_this:;
        }
    for(int j = other.rn()+1; j <= other.r(); ++j)
        {
        const Index& J = other.index(j);
        for(int k = this->rn()+1; k <= this->r(); ++k)
	    if(index(k) == J) 
		goto skip_other;
        new_index1_[++nr1_] = &J;
        skip_other:;
        }

    //
    //Special cases when one of the tensors
    //has only m==1 indices (effectively a scalar)
    //
    if(other.rn() == 0)
        {
        scale_ *= other.scale_;
        scale_ *= other.p->v(1);
        is_.r_ = is_.rn_ + nr1_;
        if(is_.r_ > NMAX) 
            {
            std::cout << "new r_ would be = " << is_.r_ << "\n";
            std::cerr << "new r_ would be = " << is_.r_ << "\n";
            Error("ITensor::operator*=: too many uncontracted indices in product (max is 8)");
            }
        //Keep current m!=1 indices, overwrite m==1 indices
        for(int j = 1; j <= nr1_; ++j) 
            is_.index_[is_.rn_+j] = *(new_index1_[j]);
        is_.setUniqueReal();
        return *this;
        }
    else if(rn() == 0)
        {
        scale_ *= other.scale_;
        scale_ *= p->v(1);
        p = other.p;
        is_.rn_ = other.is_.rn_;
        is_.r_ = is_.rn_ + nr1_;
        if(is_.r_ > NMAX) 
            {
            std::cout << "new r_ would be = " << is_.r_ << "\n";
            std::cerr << "new r_ would be = " << is_.r_ << "\n";
            Error("ITensor::operator*=: too many uncontracted indices in product (max is 8)");
            }
        for(int j = 1; j <= is_.rn_; ++j) 
            new_index_[j] = other.is_.index_[j];
        for(int j = 1; j <= nr1_; ++j) 
            new_index_[is_.rn_+j] = *(new_index1_[j]);
        is_.index_.swap(new_index_);
        is_.setUniqueReal();
        return *this;
        }

    ProductProps props(*this,other);

    int new_rn_ = 0;

    /*
    bool do_matrix_multiply = (props.odimL*props.cdim*props.odimR) > 10000 
                              || (rn_+other.rn_-2*props.nsamen) > 4 
                              || rn_ > 4 
                              || other.rn_ > 4 ;

    if(do_matrix_multiply)
        {
    */

    DO_IF_PS(++prodstats.c2;)
    MatrixRefNoLink lref, rref;
    toMatrixProd(*this,other,props,lref,rref);

    //Do the matrix multiplication
    if(p->count() != 1) { p = new ITDat(); } 
#ifdef DO_ALT
    else { p->alt.clear(); }
#endif
    p->v.ReDimension(rref.Nrows()*lref.Ncols());
    MatrixRef nref; p->v.TreatAsMatrix(nref,rref.Nrows(),lref.Ncols());
    nref = rref*lref;

    //Fill in new_index_

    if((rn() + other.rn() - 2*props.nsamen + nr1_) > NMAX) 
        {
        Print(*this);
        Print(other);
        Print(props.nsamen);
        cerr << "new m==1 indices\n";
        for(int j = 1; j <= nr1_; ++j) cerr << *(new_index1_.at(j)) << "\n";
        Error("ITensor::operator*=: too many uncontracted indices in product (max is 8)");
        }

    //Handle m!=1 indices
    for(int j = 1; j <= this->rn(); ++j)
        { if(!props.contractedL[j]) new_index_[++new_rn_] = index(j); }
    for(int j = 1; j <= other.rn(); ++j)
        { if(!props.contractedR[j]) new_index_[++new_rn_] = other.index(j); }

    /*
    else
        {
        directMultiply(other,props,new_rn_,new_index_);
        }
    */

    is_.rn_ = new_rn_;

    //Put in m==1 indices
    is_.r_ = rn();
    for(int j = 1; j <= nr1_; ++j) 
        new_index_[++(is_.r_)] = *(new_index1_.at(j));

    is_.index_.swap(new_index_);
    is_.setUniqueReal();

    scale_ *= other.scale_;

    scaleOutNorm();

    return *this;
    } //ITensor::operator*=(ITensor)

void ITensor::
reshapeTo(const Permutation& P, ITensor& res) const
    {
    res.solo();

    res.is_ = is_;
    for(int k = 1; k <= r(); ++k) 
        res.is_.index_[P.dest(k)] = index(k);

    res.scale_ = scale_;

#ifdef DO_ALT //add a flag to solo to let it know not to even copy alt
    res.p->alt.clear();
#endif
    this->reshapeDat(P,res.p->v);
    }

ITensor& ITensor::
operator+=(const ITensor& other)
    {
    if(this == &other) { scale_ *= 2; return *this; }

    bool complex_this = isComplex();
    bool complex_other = other.isComplex();
    if(!complex_this && complex_other)
        {
        return (*this = (*this * ITensor::Complex_1()) + other);
        }
    if(complex_this && !complex_other) return operator+=(other * ITensor::Complex_1());

    if(fabs(is_.uniqueReal() - other.is_.uniqueReal()) > 1E-12)
        {
        cerr << format("this ur = %.10f, other.ur = %.10f\n")%is_.uniqueReal()%other.is_.uniqueReal();
        Print(*this);
        Print(other);
        Error("ITensor::operator+=: unique Reals don't match (different Index structure).");
        }

    if(this->scale_.isRealZero())
        {
        *this = other;
        return *this;
        }
    //if(other.scale_.isRealZero()) { return *this; }
    if((other.scale_/scale_).isRealZero()) { return *this; }

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
        /*
        LogNumber LNscalefac = other.scale_/scale_;
        if(LNscalefac.isRealZero())
        {
            //Other has no effect on this
            return *this;
        }
        */
        //scalefac = LNscalefac.real();
        scalefac = (other.scale_/scale_).real();
        }


#ifdef DO_ALT
    p->alt.clear();
#endif

    bool same_ind_order = true;
    for(int j = 1; j <= rn(); ++j)
    if(index(j) != other.index(j))
        { same_ind_order = false; break; }

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
    Counter c; other.initCounter(c);
    int *j[NMAX+1];
    for(int k = 1; k <= NMAX; ++k) j[P.dest(k)] = &(c.i[k]);
    static int n[NMAX+1];
    for(int k = 1; k <= NMAX; ++k) 
        {
        n[P.dest(k)] = c.n[k];
        //n[k] = index_[k].m();
        }

#ifdef STRONG_DEBUG
    //Real tot_this = thisdat.sumels();
    //Real tot_othr = othrdat.sumels();
#endif

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

    MatrixRef dref; p->v.TreatAsMatrix(dref,i2.m(),i1.m());
    if(rn() == 2)
    { dref = M.t(i1==index(1)); }
    else
    { dref = M.t(); }
    }

void ITensor::
toMatrix11NoScale(const Index& i1, const Index& i2, Matrix& res) const
    {
    if(r() != 2) Error("toMatrix11: incorrect rank");
    assert(hasindex(i1));
    assert(hasindex(i2));
    res.ReDimension(i1.m(),i2.m());

    MatrixRef dref; p->v.TreatAsMatrix(dref,i2.m(),i1.m());
    if(rn() == 2)
    { res = dref.t(i1==index(1)); }
    else
    { res = dref.t(); }
    }

void ITensor::
toMatrix11(const Index& i1, const Index& i2, Matrix& res) const
    { toMatrix11NoScale(i1,i2,res); res *= scale_.real(); }

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

void ITensor::
symmetricDiag11(const Index& i1, ITensor& D, ITensor& U, Index& mid) const
    {
    int mink,maxk;
    symmetricDiag11(i1,D,U,mid,mink,maxk);
    }

void ITensor::
symmetricDiag11(const Index& i1, ITensor& D, ITensor& U, Index& mid, int& mink, int& maxk) const
    {
    assert(hasindex(i1));
    assert(hasindex(primed(i1)));
    if(r() != 2) Error("symDiag11: rank must be 2");
    const int m = i1.m();

    MatrixRef ref; 
    p->v.TreatAsMatrix(ref,m,m);
    Matrix UU(m,m);
    Vector d(m);
    ref *= -1;
    EigenValues(ref,d,UU);
    ref *= -1;
    d *= -1;

    mid = Index((mid.isNull() ? "mid" : mid.rawname()),m,mid.type());
    U = ITensor(i1,mid,UU);
    D = ITensor(mid,d);
    D.scale_ = scale_;

    maxk = 1;
    mink = m;
    }

Real Dot(const ITensor& x, const ITensor& y, bool doconj)
    {
    if(x.isComplex())
        {
        ITensor res = (doconj ? conj(x) : x); res *= y;
        if(res.r() != 1) Error("Bad Dot 234234");
        return res.val0();
        }
    else if(y.isComplex())
        {
        ITensor res = x; res *= y;
        if(res.r() != 1) Error("Bad Dot 37298789");
        return res.val0();
        }

    ITensor res = x; res *= y;
    if(res.r() != 0) 
        { x.print("x"); y.print("y"); Error("bad Dot"); }
    return res.val0();
    }

void Dot(const ITensor& x, const ITensor& y, Real& re, Real& im, 
                bool doconj)
    {
    if(x.isComplex())
        {
        ITensor res = (doconj ? conj(x) : x); res *= y;
        if(res.r() != 1) error("Bad Dot 334234");
        re = res(Index::IndReIm()(1));
        im = res(Index::IndReIm()(2));
        return;
        }
    else if(y.isComplex())
        {
        ITensor res = x; res *= y;
        if(res.r() != 1) error("Bad Dot 47298789");
        re = res(Index::IndReIm()(1));
        im = res(Index::IndReIm()(2));
        return;
        }
    if(x.r() != y.r()) 
        {
        cerr << "x = " << x << "\n";
        cerr << "y = " << y << "\n";
        Error("bad Dot 122414");
        }
    ITensor res = x; res *= y;
    if(res.r() != 0) 
        {
        cerr << "x = " << x << "\n";
        cerr << "y = " << y << "\n";
        Error("bad Dot 20234");
        }
    re = res.val0();
    im = 0;
    }

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
