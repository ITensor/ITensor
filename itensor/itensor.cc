//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#include "itensor.h"
#include "lapack_wrap.h"
#include "simplematrix.h"

namespace itensor {

using std::ostream;
using std::cout;
using std::endl;
using std::sqrt;
using std::vector;

#ifdef DEBUG
#define ITENSOR_CHECK_NULL if(type_ == Null) Error("ITensor is null");
#else
#define ITENSOR_CHECK_NULL
#endif

void
vectormult(std::vector<Real>& v, Real fac)
    {
    dscal_wrapper(v.size(),fac,v.data());
    }

// Y += fac*X
void
vectordaxpy(std::vector<Real>& Y, 
            const std::vector<Real>& X,
            Real fac = 1)
    {
#ifdef DEBUG
    if(X.size() != Y.size()) Error("Mismatched sizes in vectordaxpy");
#endif
    daxpy_wrapper(Y.size(),fac,X.data(),1,Y.data(),1);
    }

//Both arguments and return value of _ind
//are zero-indexed
int
_ind(const IndexSet<Index>& is,
     int i1, int i2, int i3, int i4, 
     int i5, int i6, int i7, int i8);


int static
IT_TypeToInt(ITensor::Type t)
    {
    if(t == ITensor::Null)
        {
        return 0;
        }
    else
    if(t == ITensor::Dense)
        {
        return 1;
        }
    else
    if(t == ITensor::Diag)
        {
        return 2;
        }
    else
        {
        Error("Unrecognized value of ITensor::Type enum");
        }
    return 0;
    }

ITensor::Type static
IT_IntToType(int n)
    {
    if(n == 0)
        {
        return ITensor::Null;
        }
    else
    if(n == 1)
        {
        return ITensor::Dense;
        }
    else
    if(n == 2)
        {
        return ITensor::Diag;
        }
    else
        {
        Error("Unrecognized integer for conversion to ITensor::Type.");
        }
    return ITensor::Null;
    }

namespace detail {
void
reshape(const Permutation& P, 
        const IndexSet<Index>& is, 
        const vector<Real>& dat, 
        Real *res)
    {

    const auto& ind = P.store();

    //Make a counter for dat
    Counter c(is);
    array<int,NMAX+1> n;
    for(int j = 1; j <= c.rn; ++j) n[ind[j]] = c.n[j];

    //Special case loops
#define Loop6(q,z,w,k,y,s) {for(int i1 = 0; i1 < n[1]; ++i1) for(int i2 = 0; i2 < n[2]; ++i2)\
	for(int i3 = 0; i3 < n[3]; ++i3) for(int i4 = 0; i4 < n[4]; ++i4) for(int i5 = 0; i5 < n[5]; ++i5)\
    for(int i6 = 0; i6 < n[6]; ++i6)\
    res[ (((((i6)*n[5]+i5)*n[4]+i4)*n[3]+i3)*n[2]+i2)*n[1]+i1 ] =\
    dat[ (((((s)*c.n[5]+y)*c.n[4]+k)*c.n[3]+w)*c.n[2]+z)*c.n[1]+q ]; return; }

#define Loop5(q,z,w,k,y) {for(int i1 = 0; i1 < n[1]; ++i1) for(int i2 = 0; i2 < n[2]; ++i2)\
	for(int i3 = 0; i3 < n[3]; ++i3) for(int i4 = 0; i4 < n[4]; ++i4) for(int i5 = 0; i5 < n[5]; ++i5)\
    res[ ((((i5)*n[4]+i4)*n[3]+i3)*n[2]+i2)*n[1]+i1 ] = dat[ ((((y)*c.n[4]+k)*c.n[3]+w)*c.n[2]+z)*c.n[1]+q ]; return; }

#define Loop4(q,z,w,k) {for(int i1 = 0; i1 < n[1]; ++i1)  for(int i2 = 0; i2 < n[2]; ++i2)\
	for(int i3 = 0; i3 < n[3]; ++i3) for(int i4 = 0; i4 < n[4]; ++i4)\
	res[ (((i4)*n[3]+i3)*n[2]+i2)*n[1]+i1 ] = dat[ (((k)*c.n[3]+w)*c.n[2]+z)*c.n[1]+q ]; return; }

#define Loop3(q,z,w) {for(int i1 = 0; i1 < n[1]; ++i1)  for(int i2 = 0; i2 < n[2]; ++i2)\
	for(int i3 = 0; i3 < n[3]; ++i3) res[ ((i3)*n[2]+i2)*n[1]+i1 ] = dat[ ((w)*c.n[2]+z)*c.n[1]+q ]; return; }

#define Bif3(a,b,c) if(ind[1] == a && ind[2] == b && ind[3] == c)

#define Bif4(a,b,c,d) if(ind[1] == a && ind[2] == b && ind[3] == c && ind[4] == d)

#define Bif5(a,b,c,d,e) if(ind[1] == a && ind[2] == b && ind[3] == c && ind[4]==d && ind[5] == e)

#define Bif6(a,b,c,d,e,g) if(ind[1] == a && ind[2] == b && ind[3] == c && ind[4]==d && ind[5] == e && ind[6] == g)

    if(is.rn() == 2 && ind[1] == 2 && ind[2] == 1)
        {
        MatrixRef xref; 
        VectorRefNoLink vref(const_cast<Real*>(dat.data()),dat.size());
        vref.TreatAsMatrix(xref,c.n[2],c.n[1]);
        auto newvref = Matrix(xref.t()).TreatAsVector();
        std::copy(newvref.begin(),newvref.end(),res);
        return; 
        }
    else if(is.rn() == 3)
        {
        //DO_IF_PS(int idx = ((ind[1])*3+ind[2])*3+ind[3]; Prodstats::stats().perms_of_3[idx] += 1; )
        //Arranged loosely in order of frequency of occurrence
        Bif3(2,1,3) Loop3(i2,i1,i3)
        Bif3(2,3,1) Loop3(i2,i3,i1) //cyclic
        Bif3(3,1,2) Loop3(i3,i1,i2)
        //Bif3(1,3,2) Loop3(i1,i3,i2)
        //Bif3(3,2,1) Loop3(i3,i2,i1)
        }
    else if(is.rn() == 4)
        {
        //DO_IF_PS(int idx = (((ind[1])*4+ind[2])*4+ind[3])*4+ind[4]; Prodstats::stats().perms_of_4[idx] += 1; )
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
    else if(is.rn() == 5)
        {
        //DO_IF_PS(int idx = ((((ind[1])*5+ind[2])*5+ind[3])*5+ind[4])*5+ind[5]; Prodstats::stats().perms_of_5[idx] += 1; )
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
    else if(is.rn() == 6)
        {
        //DO_IF_PS(int idx = (((((ind[1])*6+ind[2])*6+ind[3])*6+ind[4])*6+ind[5])*6+ind[6]; Prodstats::stats().perms_of_6[idx] += 1; )
        //Arranged loosely in order of frequency of occurrence
        Bif6(2,4,1,3,5,6) Loop6(i2,i4,i1,i3,i5,i6)
        Bif6(1,4,2,3,5,6) Loop6(i1,i4,i2,i3,i5,i6)
        Bif6(2,4,1,5,3,6) Loop6(i2,i4,i1,i5,i3,i6)
        Bif6(1,2,4,5,3,6) Loop6(i1,i2,i4,i5,i3,i6)
        Bif6(3,4,1,5,6,2) Loop6(i3,i4,i1,i5,i6,i2)
        }
    //DO_IF_PS(Prodstats::stats().c4 += 1;)

    //The j's are pointers to the i's of xdat's Counter,
    //but reordered in a way appropriate for res
    array<int*,NMAX+1> j;
    for(int k = 1; k <= NMAX; ++k) 
        { 
        j[ind[k]] = &(c.i[k]); 
        }

    //Catch-all loops that work for any tensor
    switch(c.rn)
    {
    case 2:
        for(; c.notDone(); ++c)
            {
            res[(*j[2])*n[1]+*j[1]]
                = dat[c.ind];
            }
        return;
    case 3:
        for(; c.notDone(); ++c)
            {
            res[((*j[3])*n[2]+*j[2])*n[1]+*j[1]]
                = dat[c.ind];
            }
        return;
    case 4:
        for(; c.notDone(); ++c)
            {
            res[(((*j[4])*n[3]+*j[3])*n[2]+*j[2])*n[1]+*j[1]]
                = dat[c.ind];
            }
        return;
    case 5:
        for(; c.notDone(); ++c)
            {
            res[((((*j[5])*n[4]+*j[4])*n[3]+*j[3])*n[2]+*j[2])*n[1]+*j[1]]
                = dat[c.ind];
            }
        return;
    case 6:
        for(; c.notDone(); ++c)
            {
            res[(((((*j[6])*n[5]+*j[5])*n[4]+*j[4])*n[3]+*j[3])*n[2]+*j[2])*n[1]+*j[1]]
                = dat[c.ind];
            }
        return;
    case 7:
        for(; c.notDone(); ++c)
            {
            res[((((((*j[7])*n[6]+*j[6])*n[5]+*j[5])*n[4]+*j[4])*n[3]+*j[3])*n[2]+*j[2])*n[1]+*j[1]]
                = dat[c.ind];
            }
        return;
    default:
        for(; c.notDone(); ++c)
            {
            res[(((((((*j[8])*n[7]+*j[7])*n[6]+*j[6])*n[5]+*j[5])*n[4]+*j[4])*n[3]+*j[3])*n[2]+*j[2])*n[1]+*j[1]]
                = dat[c.ind];
            }
        return;
    } //switch(c.rn)

    } // reshape
};

void
reshape(const Permutation& P, 
        const IndexSet<Index>& is, 
        const vector<Real>& dat, 
        vector<Real>& res)
    {
    if(P.isTrivial())
        {
        res = dat;
        return;
        }
    res.resize(dat.size());
    detail::reshape(P,is,dat,res.data());
    }

void
reshape(const Permutation& P, 
        const IndexSet<Index>& is, 
        const vector<Real>& dat, 
        Vector& res)
    {
    if(P.isTrivial())
        {
        res = Vector(dat);
        return;
        }
    res.ReDimension(dat.size());
    detail::reshape(P,is,dat,res.Store());
    }

//
// ITensor
//

//
// ITensor Constructors
//

ITensor::
ITensor()  
    : 
    type_(Null)
    { }


ITensor::
ITensor(Real val) 
    :
    type_(Dense),
    scale_(1)
    { 
    allocate(1,val);
    }

ITensor::
ITensor(const Index& i1) 
    :
    type_(Dense),
    is_(i1),
    scale_(1)
	{ 
    allocate(i1.m());
    }

ITensor::
ITensor(const Index& i1, Real val) 
    :
    type_(Dense),
    is_(i1),
    scale_(1)
	{ 
    allocate(i1.m(),val);
    }

ITensor::
ITensor(const Index& i1, const VectorRef& V) 
    : 
    type_(Dense),
    r_(make_shared<ITDat>(V)),
    is_(i1),
    scale_(1)
	{ 
	if(i1.m() != V.Length()) 
	    Error("Mismatch of Index and Vector sizes.");
	}

ITensor::
ITensor(const Index& i1,const Index& i2) 
    :
    type_(Dense),
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
    type_ = Diag;
    const int dim = min(i1.m(),i2.m());
    allocate(dim,a);
	}

ITensor::
ITensor(const Index& i1,const Index& i2, const VectorRef& V) 
    :
    r_(make_shared<ITDat>(V)),
    is_(i1,i2),
    scale_(1)
	{
    type_ = Diag;
#ifdef DEBUG
    auto dim = min(i1.m(),i2.m());
    if(V.Length() != dim)
        Error("Diagonal vector must have length == min(i1.m(),i2.m())");
#endif
	}

ITensor::
ITensor(const Index& i1,const Index& i2,const MatrixRef& M) 
    :
    type_(Dense),
    is_(i1,i2),
    scale_(1)
	{
    allocate(i1.m()*i2.m());
	if(i1.m() != M.Nrows() || i2.m() != M.Ncols()) 
	    Error("Mismatch of Index sizes and matrix.");
	MatrixRef dref; 
    VectorRefNoLink vref(r_->data(),r_->size());
	vref.TreatAsMatrix(dref,i2.m(),i1.m()); 
	dref = M.t();
	}

ITensor::
ITensor(const Index& i1, const Index& i2, const Index& i3,
        const Index& i4, const Index& i5, const Index& i6,
        const Index& i7, const Index& i8)
    :
    type_(Dense),
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
ITensor(const IndexVal& iv)
    :
    type_(Dense),
    is_(iv.index),
    scale_(1)
	{ 
    allocate(iv.m());
#ifdef DEBUG
    if(iv.i < 1 || iv.i > int(r_->size())) 
        {
        Print(iv.i);
        Print(r_->size());
        Error("IndexVal out of range");
        }
#endif
	r_->v[iv.i-1] = 1; 
	}

ITensor::
ITensor(const IndexVal& iv1, const IndexVal& iv2) 
    :
    type_(Dense),
    is_(iv1.index,iv2.index),
    scale_(1)
	{ 
    allocate(iv1.m()*iv2.m());
    auto offset = (iv2.i-1)*iv1.m()+iv1.i-1;
#ifdef DEBUG
    if(offset < 0 || offset > r_->size()) Error("IndexVal out of range");
#endif
	r_->v[offset] = 1; 
	}

ITensor::
ITensor(const IndexVal& iv1, const IndexVal& iv2, 
        const IndexVal& iv3, const IndexVal& iv4, 
        const IndexVal& iv5, const IndexVal& iv6, 
        const IndexVal& iv7, const IndexVal& iv8)
    :
    type_(Dense),
    scale_(1)
	{
    //Construct ITensor
    array<Index,NMAX+1> ii = 
        {{ iv1.index, iv2.index, iv3.index, iv4.index, iv5.index, 
           iv6.index, iv7.index, iv8.index, Index::Null()}};
    int size = 3; 
    while(ii[size] != Index::Null()) ++size;
    int alloc_size = -1;
    is_ = IndexSet<Index>(ii,size,alloc_size,0);
    allocate(alloc_size);

    //Assign specified element to 1
    array<int,NMAX+1> iv = 
        {{ iv1.i, iv2.i, iv3.i, iv4.i, iv5.i, iv6.i, iv7.i, iv8.i, 1 }};
    array<int,NMAX> ja; 
    ja.fill(0);
    for(int k = 0; k < is_.rn(); ++k) //loop over indices of this ITensor
    for(int j = 0; j < size; ++j)      // loop over the given indices
        {
        if(is_[k] == ii[j]) 
            { ja[k] = iv[j]-1; break; }
        }

    r_->v[_ind(is_,ja[0],ja[1],ja[2],ja[3],ja[4],ja[5],ja[6],ja[7])] = 1;
    }

ITensor::
ITensor(const IndexSet<Index>& I) 
    :
    type_(Dense),
    is_(I),
    scale_(1)
	{
	allocate(is_.dim());
	}

ITensor::
ITensor(const IndexSet<Index>& I, const Vector& V) 
    : 
    type_(Dense),
    r_(make_shared<ITDat>(V)),
    is_(I),
    scale_(1)
	{
#ifdef DEBUG
    if(is_.dim() != V.Length())
	    Error("incompatible Index and Vector sizes");
#endif
	}


ITensor::
ITensor(const IndexSet<Index>& I, const ITensor& other) 
    : 
    type_(other.type_),
    r_(other.r_), 
    i_(other.i_), 
    is_(I),
    scale_(other.scale_)
	{
#ifdef DEBUG
	if(is_.dim() != other.is_.dim()) 
	    { Error("incompatible dimensions"); }
#endif
	}

ITensor::
ITensor(const IndexSet<Index>& I, const ITensor& other, const Permutation& P) 
    : 
    type_(other.type_),
    is_(I),
    scale_(other.scale_)
    {
#ifdef DEBUG
    if(is_.dim() != other.is_.dim()) 
        Error("incompatible dimensions");
#endif
    if(P.isTrivial()) 
        { 
        r_ = other.r_; 
        i_ = other.i_; 
        }
    else               
        { 
        allocate(); 
        reshape(P,other.is_,other.r_->v,r_->v); 

        if(other.i_)
            {
            allocateImag();
            reshape(P,other.is_,other.i_->v,i_->v); 
            }
        }
    }

ITensor& ITensor::
takeRealPart()
    {
    ITENSOR_CHECK_NULL
    i_.reset();
    return *this;
    }

ITensor& ITensor::
takeImagPart()
    {
    ITENSOR_CHECK_NULL
    if(!i_)
        {
        scale_ = LogNumber(0);
        return *this;
        }
    r_.swap(i_);
    i_.reset();
    return *this;
    }

Vector ITensor::
diag() const
    {
    if(this->isComplex()) 
        Error("diag() may only be called on real ITensors - try taking real or imaginary part first");
    Vector res;
    if(type_ == Diag)
        {
        res = Vector(r_->v);
        }
    else
    if(type_ == Dense)
        {
        res = Vector(minM(is_));
        for(int i = 0; i < res.Length(); ++i)
            {
            res(1+i) = r_->v[_ind(is_,i,i,i,i,i,i,i,i)];
            }
        }
    else
        {
        Error("diag: null ITensor");
        }

    if(scale_.isTooBigForReal())
        throw TooBigForReal("Scale too large for real in ITensor::diag()");
    res *= scale_.real0();

    return res;
    }

void ITensor::
read(std::istream& s)
    { 
    int tint = 0;
    s.read((char*) &tint,sizeof(tint));
    type_ = IT_IntToType(tint);

    if(type_ == Null) { *this = ITensor(); return; }

    is_.read(s);
    scale_.read(s);
    r_ = make_shared<ITDat>();
    r_->read(s);
    bool is_cplx = false;
    s.read((char*)&is_cplx,sizeof(is_cplx));
    if(is_cplx)
        {
        i_ = make_shared<ITDat>();
        i_->read(s);
        }
    }

void ITensor::
write(std::ostream& s) const 
    { 
    int tint = IT_TypeToInt(type_);
    s.write((char*) &tint,sizeof(tint));

    if(type_ == Null) return;

    is_.write(s);
    scale_.write(s);
    r_->write(s);
    bool is_cplx = isComplex();
    s.write((char*)&is_cplx,sizeof(is_cplx));
    if(is_cplx) i_->write(s);
    }


Real ITensor::
toReal() const 
	{ 
#ifdef DEBUG
    if(!this->valid())
        Error("ITensor is null");
#endif

    if(isComplex())
        Error("ITensor is complex");

    if(is_.rn() != 0)
        {
        Print(*this);
        Error("ITensor not a scalar");
        }

	try {
	    return r_->v.at(0)*scale_.real(); 
	    }
	catch(const TooBigForReal& e)
	    {
	    println("too big for real() in toReal");
	    println("r_->v[0] is ",r_->v.at(0));
	    println("scale is ",scale());
	    println("rethrowing");
	    throw e;
	    }
	catch(TooSmallForReal)
	    {
	    //cout << "warning: too small for real() in toReal" << endl;
	    //cout << "r_->v(1) is " << r_->v(1) << endl;
	    //cout << "scale is " << scale() << endl;
	    return 0.;
	    }
	return NAN; //shouldn't reach this line
	}

Complex ITensor::
toComplex() const
	{ 
    if(this->isComplex())
        {
        Real re, im;
        try {
            re = r_->v.at(0)*scale_.real(); 
            }
        catch(const TooBigForReal& e)
            {
            println("too big for real() in toReal");
            println("r_->v[0] is ",r_->v.at(0));
            println("scale is ",scale());
            println("rethrowing");
            throw e;
            }
        catch(TooSmallForReal)
            {
            re = 0.;
            }

        try {
            im = i_->v.at(0)*scale_.real(); 
            }
        catch(const TooBigForReal& e)
            {
            println("too big for real() in toReal");
            println("i_->v[0] is ",i_->v.at(0));
            println("scale is ",scale());
            println("rethrowing");
            throw e;
            }
        catch(TooSmallForReal)
            {
            im = 0.0;
            }
        return Complex(re,im);
        }
    return Complex(toReal(),0.);
    }

/*
Real& ITensor::
operator()()
	{ 
    if(is_.rn() != 0)
        {
        printfln("# given = 0, rn_ = %d\n",is_.rn());
        Error("Not enough indices (requires all having m!=1)");
        }
    solo(); 
    scaleTo(1);
    return r_->v(1);
    }

Real ITensor::
operator()() const
	{ 
    ITENSOR_CHECK_NULL
    if(is_.rn() != 0)
        {
        printfln("# given = 0, rn_ = %d\n",is_.rn());
        Error("Not enough indices (requires all having m!=1)");
        }
    return scale_.real()*r_->v(1);
    }
    */

Real& ITensor::
operator()(const IndexVal& iv1)
	{
#ifdef DEBUG
    if(is_.rn() > 1) 
        {
        printfln("# given = 1, rn_ = %d\n",is_.rn());
        Error("Not enough m!=1 indices provided");
        }
    if(is_[0] != iv1.index)
        {
        Print(*this);
        Print(iv1);
        Error("Incorrect IndexVal argument to ITensor");
        }
#endif
    solo(); 
    scaleTo(1);
    return r_->v.at(iv1.i-1);
	}

Real ITensor::
operator()(const IndexVal& iv1) const
	{
#ifdef DEBUG
    if(is_.rn() > 1) 
        {
        printfln("# given = 1, rn_ = %d\n",is_.rn());
        Error("Not enough m!=1 indices provided");
        }
    if(is_[0] != iv1.index)
        {
        Print(*this);
        Print(iv1);
        Error("Incorrect IndexVal argument to ITensor");
        }
#endif
    return scale_.real()*r_->v.at(iv1.i-1);
	}

Real& ITensor::
operator()(const IndexVal& iv1, const IndexVal& iv2) 
    {
    solo(); 
    scaleTo(1);

    if(type_ == Diag && iv1.i != iv2.i)
        {
        convertToDense();
        }

    if(type_ == Diag)
        {
        return r_->v.at(iv1.i-1);
        }

    return r_->v[_ind2(iv1,iv2)];
    }

Real ITensor::
operator()(const IndexVal& iv1, const IndexVal& iv2) const
    {
    ITENSOR_CHECK_NULL
    if(type_ == Diag)
        {
        if(iv1.i != iv2.i) return 0;
        return r_->v.at(iv1.i-1);
        }
    return scale_.real()*r_->v[_ind2(iv1,iv2)];
    }

Real& ITensor::
operator()(const IndexVal& iv1, const IndexVal& iv2, 
           const IndexVal& iv3, const IndexVal& iv4, 
           const IndexVal& iv5, const IndexVal& iv6,
           const IndexVal& iv7, const IndexVal& iv8)
    {
    solo(); 
    scaleTo(1);
    if(type_ == Diag)
        {
        const int di = _diag_ind8(iv1,iv2,iv3,iv4,iv5,iv6,iv7,iv8);
        if(di == -1)
            convertToDense();
        else
            return r_->v.at(di-1);
        }
    return r_->v[_ind8(iv1,iv2,iv3,iv4,iv5,iv6,iv7,iv8)];
    }

Real ITensor::
operator()(const IndexVal& iv1, const IndexVal& iv2, 
           const IndexVal& iv3, const IndexVal& iv4,
           const IndexVal& iv5,const IndexVal& iv6,
           const IndexVal& iv7,const IndexVal& iv8) const
    {
    ITENSOR_CHECK_NULL
    if(type_ == Diag)
        {
        const int di = _diag_ind8(iv1,iv2,iv3,iv4,iv5,iv6,iv7,iv8);
        if(di == -1) return 0;
        return scale_.real()*r_->v.at(di-1);
        }
    return scale_.real()*r_->v[_ind8(iv1,iv2,iv3,iv4,iv5,iv6,iv7,iv8)];
    }

//Process IndexVals for element access when ITensor
//has type_ == Diag
int ITensor::
_diag_ind8(const IndexVal& iv1, const IndexVal& iv2, 
           const IndexVal& iv3, const IndexVal& iv4,
           const IndexVal& iv5, const IndexVal& iv6,
           const IndexVal& iv7, const IndexVal& iv8) const
    {
    if(iv1 == IndexVal::Null())
        Error("Null IndexVal argument");

    array<const IndexVal*,NMAX> iv = 
        {{ &iv1, &iv2, &iv3, &iv4, &iv5, &iv6, &iv7, &iv8 }};

    //Loop over the given IndexVals
    for(int j = 1; j < is_.r(); ++j)
        {
        const IndexVal& J = *iv[j];
        if(J == IndexVal::Null()) break;
        if(J.i != iv1.i) return -1; //off-diagonal, signal with a -1
        //otherwise J.i == iv1.i, continue checking all non-Null IndexVals
        }
    return iv1.i;
    }


void ITensor::
groupIndices(const array<Index,NMAX+1>& indices, int nind, 
             const Index& grouped, ITensor& res) const
    {
    if(type_ == Diag)
        {
        Error("groupIndices not yet defined for type() == Diag");
        }

    array<int,NMAX+1> isReplaced; 
    isReplaced.fill(0);

    //Print(*this);

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
            if(is_.index(k) == J) 
                {
                isReplaced[k] = (J.m() == 1 ? -1 : nn);
                foundit = true; 
                break; 
                }
            }
        if(!foundit)
            {
            Print(*this);
            println("Couldn't find Index ",J," in ITensor.");
            println("indices:");
            for(int j = 1; j <= nind; ++j)
                println("  ",indices[j]);
            Error("bad request for Index to replace");
            }
        }
    if(tot_m != grouped.m()) Error("ITensor::groupAndReplace: \
                                    mismatched index sizes.");

    //Compute rn_ of res
    const int res_rn_ = is_.rn() - nn + (nn == 0 ? 0 : 1);

    IndexSet<Index> nindices; 
    Permutation P(NMAX+1);
    int nkept = 0; 
    for(int j = 1; j <= is_.rn(); ++j)
        {
        if(isReplaced[j] == 0)
            {
            P.setFromTo(j,++nkept);
            nindices.addindex(is_.index(j)); 
            }
        else
            {
            P.setFromTo(j,res_rn_+isReplaced[j]-1);
            }
        }

    nindices.addindex(grouped);

    for(int j = is_.rn()+1; j <= r(); ++j) 
        if(isReplaced[j] == 0) nindices.addindex(is_.index(j));

    if(nn == 0) 
        res = ITensor(nindices,*this);
    else        
        res = ITensor(nindices,*this,P); 
    }

void ITensor::
tieIndices(const array<Index,NMAX>& indices, int nind,
           const Index& tied)
    {
    if(type_ == Diag)
        {
        Error("tieIndices not yet defined for type() == Diag");
        }

    if(nind == 0) Error("No indices given");

    const int tm = tied.m();
    
    array<Index,NMAX+1> new_index_;
    new_index_[1] = tied;
    //will count these up below
    int new_r_ = 1;
    int alloc_size = tm;

    array<bool,NMAX+1> is_tied;
    is_tied.fill(false);

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
        println("indices = ");
        for(int j = 0; j < nind; ++j) println(indices[j]);
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
    array<const int*,NMAX+1> ii;
    int n = 2;
    for(int j = 1; j <= r(); ++j)
        {
        if(is_tied[j])
            ii[j] = &(nc.i[1]);
        else
            ii[j] = &(nc.i[n++]);
        }

    const int zero = 0;
    for(int j = r()+1; j <= NMAX; ++j)
        ii[j] = &zero;
    
    //Create the new dat
    auto np = make_shared<ITDat>(alloc_size,0);
    const auto& thisdat = r_->v;
    auto& newdat = np->v;
    for(; nc.notDone(); ++nc)
        {
        newdat[nc.ind] =
        thisdat[_ind(is_,*ii[1],*ii[2],
                         *ii[3],*ii[4],
                         *ii[5],*ii[6],
                         *ii[7],*ii[8])];
        }

    r_.swap(np);

    if(this->isComplex())
        {
        np = make_shared<ITDat>(alloc_size,0);
        const auto& thisidat = i_->v;
        auto& newdat = np->v;
        for(nc.reset(); nc.notDone(); ++nc)
            {
            newdat[nc.ind] =
            thisidat[_ind(is_,*ii[1],*ii[2],
                              *ii[3],*ii[4],
                              *ii[5],*ii[6],
                              *ii[7],*ii[8])];
            }
        i_.swap(np);
        }

    is_.swap(new_is_);

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

ITensor& ITensor::
trace(const array<Index,NMAX>& indices, int nind)
    {
    if(type_ == Diag)
        {
        Error("trace not yet defined for type() == Diag");
        }

    if(nind < 0)
        {
        nind = 0;
        while(indices[nind] != Index::Null()) ++nind;
        }

    if(nind == 0) Error("No indices given");

    const int tm = indices[0].m();
    
    array<Index,NMAX+1> new_index_;

    //will count these up below
    int new_r_ = 0;
    int alloc_size = 1;

    array<bool,NMAX+1> traced;
    traced.fill(false);

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
        println("indices = ");
        for(int j = 0; j < nind; ++j) println(indices[j]);
        Error("Couldn't find Index to trace");
        }

    IndexSet<Index> new_is_(new_index_,new_r_,alloc_size,1);

    //If traced indices have m==1, no work
    //to do; just replace indices
    if(tm == 1)
        {
        is_.swap(new_is_);
        return *this;
        }

    Counter nc(new_is_);

    //Set up ii pointers to link
    //elements of res to appropriate
    //elements of *this
    int trace_ind = 0;
    array<const int*,NMAX+1> ii;
    int n = 1;
    for(int j = 1; j <= r(); ++j)
        {
        if(traced[j])
            ii[j] = &(trace_ind);
        else
            ii[j] = &(nc.i[n++]);
        }

    const int zero = 0;
    for(int j = r()+1; j <= NMAX; ++j)
        ii[j] = &zero;
    
    //Create the new dat
    shared_ptr<ITDat> np = make_shared<ITDat>(alloc_size,0);
    auto& resdat = np->v;

    const auto& thisdat = r_->v;
    for(; nc.notDone(); ++nc)
        {
        Real newval = 0;
        for(trace_ind = 0; trace_ind < tm; ++trace_ind)
            {
            newval += 
            thisdat[_ind(is_,*ii[1],*ii[2],
                             *ii[3],*ii[4],
                             *ii[5],*ii[6],
                             *ii[7],*ii[8])];
            }
        resdat[nc.ind] = newval;
        }

    is_.swap(new_is_);
    r_.swap(np);

    if(this->isComplex())
        {
        const auto& thisidat = i_->v;
        for(nc.reset(); nc.notDone(); ++nc)
            {
            Real newval = 0;
            for(trace_ind = 0; trace_ind < tm; ++trace_ind)
                {
                newval += 
                thisidat[_ind(is_,*ii[1],*ii[2],
                                  *ii[3],*ii[4],
                                  *ii[5],*ii[6],
                                  *ii[7],*ii[8])];
                }
            resdat[nc.ind] = newval;
            }

        i_.swap(np);
        }

    return *this;
    } //ITensor::trace


ITensor& ITensor::
trace(const Index& i1, const Index& i2,
      const Index& i3, const Index& i4,
      const Index& i5, const Index& i6,
      const Index& i7, const Index& i8)
    {
    array<Index,NMAX> inds = {{ i1, i2, i3, i4,
                                i5, i6, i7, i8 }};
    trace(inds);
    return *this;
    }

void ITensor::
expandIndex(const Index& small, const Index& big, int start)
    {
    if(type_ == Diag)
        {
        Error("expandIndex not yet defined for type() == Diag");
        }

    if(this->isComplex())
        {
        ITensor r(*this),
                i(*this);
        r.takeRealPart();
        i.takeImagPart();
        r.expandIndex(small,big,start);
        i.expandIndex(small,big,start);
        i *= Complex_i;
        *this = r+i;
        return;
        }

#ifdef DEBUG
    if(small.m() > big.m())
        {
        Print(small);
        Print(big);
        Error("small Index must have smaller m() than big Index");
        }
    if(start >= big.m())
        {
        Print(start);
        Print(big.m());
        Error("start must be less than big.m()");
        }
#endif

    if(this->isComplex())
        {
        ITensor r(*this),
                i(*this);
        r.takeRealPart();
        i.takeImagPart();
        r.expandIndex(small,big,start);
        i.expandIndex(small,big,start);
        i *= Complex_i;
        *this = r+i;
        return;
        }

    IndexSet<Index> newinds; 
    bool found = false;
    for(int j = 1; j <= r(); ++j)
        {
        if(is_.index(j) == small)
            {
            newinds.addindex(big);
            found = true;
            }
        else 
            {
            newinds.addindex(is_.index(j));
            }
        }

    if(!found)
        {
        Print(*this);
        Print(small);
        Error("couldn't find index");
        }


    const int w = findindex(newinds,big);
    int inc = start;
    for(int n = 0; n < w; ++n)
        {
        inc *= newinds[n].m();
        }

    //Comparing nmax and omax determines whether
    //old dat fits into new dat sequentially, in which
    //case we can use std::copy
    const
	int nmax = 1+_ind(newinds,is_[0].m()-1,is_[1].m()-1, 
                              is_[2].m()-1,is_[3].m()-1, 
                              is_[4].m()-1,is_[5].m()-1, 
                              is_[6].m()-1,is_[7].m()-1);

    shared_ptr<ITDat> oldr(r_);
    allocate(newinds.dim());

    auto omax = oldr->v.size();
    const Real* const olddat = oldr->data();
    Real* const newdat = r_->data();

	if(nmax == omax)
	    {
#ifdef DEBUG
        if((inc+omax) > r_->size()) 
            {
            Print(inc);
            Print(omax);
            Print(inc+omax);
            Print(r_->size());
            Error("Mismatched sizes in expandIndex, copy case");
            }
#endif
        std::copy(olddat,olddat+omax,newdat+inc);
	    }
    else
        {
        Counter c(is_);
        for(; c.notDone(); ++c)
            {
#ifdef DEBUG
            auto ind = inc+_ind(newinds,c.i[1],c.i[2],
                                    c.i[3],c.i[4],
                                    c.i[5],c.i[6],
                                    c.i[7],c.i[8]);
            if(ind < 0 || ind >= int(r_->size())) Error("ind out of range in expandIndex");
            if(c.ind >= int(oldr->size())) Error("c.ind out of range in expandIndex");
#endif
            newdat[inc+_ind(newinds,c.i[1],c.i[2],
                                    c.i[3],c.i[4],
                                    c.i[5],c.i[6],
                                    c.i[7],c.i[8])]
            = olddat[c.ind];
            }
        }

    is_.swap(newinds);
    }


VectorRef ITensor::
assignToVec() const
    {
    ITENSOR_CHECK_NULL
    if(this->isComplex())
        Error("assignToVec defined only for real ITensor");
    Vector rv(r_->v);
    if(scale_.isRealZero()) 
        {
        rv *= 0;
        return rv;
        }
    rv *= scale_.real();
    return rv;
    }

void ITensor::
pseudoInvert(Real cutoff)
    {
    if(type_ != Diag)
        Error("pseudoInvert only defined for ITensor of type()==ITensor::Diag");
    if(this->isComplex())
        Error("pseudoInvert currently only defined for real ITensor");
    solo();
    scale_.pow(-1); //succeeds even if scale_ == 0
    for(int j = 1; j <= r_->size(); ++j)
        {
        if(r_->v.at(j-1) > cutoff)
            r_->v.at(j-1) = 1./r_->v.at(j-1);
        else
            r_->v.at(j-1) = 0;
        }
    }


void ITensor::
reshapeDat(const Permutation& P)
    {
    if(type_ == Diag)
        {
        Error("reshapeDat not yet defined for type() == Diag");
        }
    if(P.isTrivial()) return;
    solo();
    std::vector<Real> newdat;
    reshape(P,is_,r_->v,newdat);
    r_->v = std::move(newdat);
    if(i_)
        {
        reshape(P,is_,i_->v,newdat);
        i_->v = std::move(newdat);
        }
    }


void ITensor::
swap(ITensor& other)
    {
    const Type tmp = type_;
    type_ = other.type_;
    other.type_ = tmp;

    r_.swap(other.r_);
    i_.swap(other.i_);
    is_.swap(other.is_);
    scale_.swap(other.scale_);
    }

const Real* ITensor::
datStart() const
    {
    if(!r_) Error("ITensor is null");
    return r_->data();
    }

const Real* ITensor::
imagDatStart() const
    {
    if(!i_) Error("ITensor is real");
    return i_->data();
    }

void ITensor::
randomize(const Args& args) 
    { 
    solo(); 
    convertToDense();
    for(size_t j = 0; j < r_->v.size(); ++j)
        {
        r_->v[j] = Global::random();
        }
    if(i_ || args.getBool("Complex",false))
        {
        allocateImag(r_->v.size());
        for(size_t j = 0; j < i_->v.size(); ++j)
            {
            i_->v[j] = Global::random();
            }
        }
    }

void ITensor::
conj() 
    { 
    if(i_)
        {
        soloImag();
        vectormult(i_->v,-1);
        }
    }

Real
sumels(const ITensor& t)
    { 
    if(t.isComplex())
        Error("sumels defined only for real ITensor");
    return t.assignToVec().sumels();
    }

Real ITensor::
normNoScale() const 
    { 
    ITENSOR_CHECK_NULL
    if(!this->isComplex())
        {
        return Norm(VectorRefNoLink(r_->data(),r_->size()));
        }
    else
        {
        auto rref = VectorRefNoLink(r_->data(),r_->size());
        auto iref = VectorRefNoLink(i_->data(),i_->size());
        return sqrt(sqr(Norm(rref))+sqr(Norm(iref)));
        }
    }

Real ITensor::
norm() const 
    { 
    if(scale_.isTooBigForReal())
        {
        throw TooBigForReal("Scale too large for real in ITensor::norm()");
        }
    //If scale_ is too small to be converted to Real,
    //real0 method will return 0.0
    return fabs(scale_.real0())*normNoScale();
    }

LogNumber ITensor::
normLogNum() const 
    { 
    return LogNumber(log(normNoScale())+scale_.logNum(),+1);
    }


void ITensor::
scaleOutNorm()
    {
    Real f = normNoScale();
    //If norm already 1 return so
    //we don't have to call solo()
    if(fabs(f-1) < 1E-12) return;

    if(f != 0) 
        { 
        solo();
        vectormult(r_->v,1./f);
        if(i_) vectormult(i_->v,1./f);
        scale_ *= f; 
        }
    else //norm == zero
        {
        scale_ = LogNumber(1.);
        i_.reset();
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
    vectormult(r_->v,scale_.real0());
    if(i_) vectormult(i_->v,scale_.real0());
    scale_ = newscale;
    }


void ITensor::
allocate(size_t dim, Real val) 
    { 
    r_ = make_shared<ITDat>(dim,val);
    }

void ITensor::
allocate() 
    { 
    r_ = make_shared<ITDat>(); 
    }

void ITensor::
allocateImag(size_t dim,Real val) 
    { 
    i_ = make_shared<ITDat>(dim,val); 
    }

void ITensor::
allocateImag() 
    { 
    i_ = make_shared<ITDat>(); 
    }

void ITensor::
soloReal()
	{
    ITENSOR_CHECK_NULL
    if(!r_.unique())
        { 
        shared_ptr<ITDat> newr = make_shared<ITDat>();
        newr->v = r_->v;
        r_.swap(newr);
        }
    }

void ITensor::
soloImag()
    {
    if(!i_) return;

    if(!i_.unique())
        { 
        shared_ptr<ITDat> newi = make_shared<ITDat>();
        newi->v = i_->v;
        i_.swap(newi);
        }
	}

void ITensor::
solo()
    {
    soloReal();
    soloImag();
    }

void ITensor::
equalizeScales(ITensor& other)
    {
    if(scale_.sign() != 0)
        {
        other.scaleTo(scale_);
        }
    else //*this is equivalent to zero
        {
        soloReal();
        std::fill(r_->v.begin(),r_->v.end(),0.);
        scale_ = other.scale_;
        }
    }

ITensor& ITensor::
operator*=(Real fac)
    {
    if(fac == 0)
        {
        solo();
        std::fill(r_->v.begin(),r_->v.end(),0.);
        if(i_) std::fill(i_->v.begin(),i_->v.end(),0);
        return *this;
        }

    scale_ *= fac;
    return *this;
    }

ITensor& ITensor::
operator*=(Complex z)
    {
    if(z.real() == 0)
        {
        r_.swap(i_);
        if(!r_) allocate(i_->v.size());
        soloReal();
        vectormult(r_->v,-1);
        scale_ *= z.imag();
        return *this;
        }
    else
    if(z.imag() == 0)
        {
        operator*=(z.real());
        return *this;
        }
    else
    if(!this->isComplex())
        {
        allocateImag();
        i_->v = r_->v;
        if(fabs(z.real()) > fabs(z.imag()))
            {
            scale_ *= z.real();
            vectormult(i_->v,(z.imag()/z.real()));
            }
        else
            {
            soloReal();
            vectormult(r_->v,(z.real()/z.imag()));
            scale_ *= z.imag();
            }
        return *this;
        }

    //Else this is complex
    solo();
    VectorRefNoLink rref(r_->data(),r_->size());
    VectorRefNoLink iref(i_->data(),i_->size());
    //Vector newr = r_->v*z.real() - i_->v*z.imag();
    //Vector newi = r_->v*z.imag() + i_->v*z.real();
    Vector newr = rref*z.real()-iref*z.imag();
    Vector newi = rref*z.imag()+iref*z.real();
    r_->v.assign(newr.begin(),newr.end());
    i_->v.assign(newi.begin(),newi.end());

    return *this;
    }


int
_ind(const IndexSet<Index>& is,
     int i1, int i2, int i3, int i4, 
     int i5, int i6, int i7, int i8)
    {
    switch(is.rn())
    {
    case 0:
        return 0;
    case 1:
        return (i1);
    case 2:
        return ((i2)*is[0].m()+i1);
    case 3:
        return (((i3)*is[1].m()+i2)*is[0].m()+i1);
    case 4:
        return ((((i4)*is[2].m()+i3)*is[1].m()+i2)*is[0].m()+i1);
    case 5:
        return (((((i5)*is[3].m()+i4)*is[2].m()+i3)*is[1].m()+i2)
                        *is[0].m()+i1);
    case 6:
        return ((((((i6)*is[4].m()+i5)*is[3].m()+i4)*is[2].m()+i3)
                        *is[1].m()+i2)*is[0].m()+i1);
    case 7:
        return (((((((i7)*is[5].m()+i6)*is[4].m()+i5)*is[3].m()+i4)
                        *is[2].m()+i3)*is[1].m()+i2)*is[0].m()+i1);
    case 8:
        return ((((((((i8)*is[6].m()+i7)*is[5].m()+i6)*is[4].m()+i5)
                        *is[3].m()+i4)*is[2].m()+i3)*is[1].m()+i2)*is[0].m()+i1);
    } //switch
    Error("_ind: Failed switch case");
    return 0;
    }


int ITensor::
_ind2(const IndexVal& iv1, const IndexVal& iv2) const
    {
    if(is_.rn() > 2) 
        {
        printfln("# given = 2, rn_ = %d\n",is_.rn());
        Error("Not enough m!=1 indices provided");
        }
    if(is_[0] == iv1.index && is_[1] == iv2.index)
        return ((iv2.i-1)*is_[0].m()+iv1.i-1);
    else if(is_[0] == iv2.index && is_[1] == iv1.index)
        return ((iv1.i-1)*is_[0].m()+iv2.i-1);
    else
        {
        Print(*this);
        Print(iv1);
        Print(iv2);
        Error("Incorrect IndexVal argument to ITensor");
        return 0;
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
    array<int,NMAX> ja; 
    ja.fill(0);
    //Loop over the given IndexVals
    int nn = 0;
    for(int j = 0; j < is_.r(); ++j)
        {
        const IndexVal& J = *iv[j];
        if(J == IndexVal::Null()) break;
        if(J.m() != 1) ++nn;
        if(J.i == 1) continue;
        //Loop over indices of this ITensor
        for(int k = 0; k < is_.r(); ++k)
            {
            if(is_[k] == J.index)
                {
                ja[k] = J.i-1;
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
    Permutation pl,
                pr;

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
    rcstart(100),
    pl(NMAX+1),
    pr(NMAX+1)
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
		pl.setFromTo(j,nsamen);
		pr.setFromTo(k,nsamen);

		contractedL[j] = contractedR[k] = true;

        cdim *= L.is_.index(j).m();

        //matchL.setFromTo(k,j-lcstart+1);
		}
    //Finish making pl
    int q = nsamen;
    for(int j = 1; j <= L.is_.rn(); ++j)
        if(!contractedL[j]) pl.setFromTo(j,++q);
    //Finish making pr and matchL
    q = nsamen;
    for(int j = 1; j <= R.is_.rn(); ++j)
        if(!contractedR[j]) 
            {
            ++q;
            pr.setFromTo(j,q);
            //matchL.setFromTo(j,q);
            }

    odimL = L.r_->size()/cdim;
    odimR = R.r_->size()/cdim;
    }

//Converts ITensor dats into MatrixRef's that can be multiplied as rref*lref
//contractedL/R[j] == true if L/R.indexn(j) contracted
void 
toMatrixProd(const ITensor& L, const ITensor& R, 
             std::vector<Real>& newLdat,
             std::vector<Real>& newRdat,
             ProductProps& props,
             SimpleMatrixRef& lref, SimpleMatrixRef& rref, 
             bool& L_is_matrix, bool& R_is_matrix, bool doReshape)
    {
#ifdef DEBUG
    if(L.type() == ITensor::Diag)
        Error("toMatrixProd not implemented for ITensor of type Diag (L)");
    if(R.type() == ITensor::Diag)
        Error("toMatrixProd not implemented for ITensor of type Diag (R)");
    if(!L) Error("L null in toMatrixProd");
    if(!R) Error("R null in toMatrixProd");
#endif
    const auto &Ldat = L.r_->v;
    const auto &Rdat = R.r_->v;

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
            lref = SimpleMatrixRef(Ldat.data(),props.odimL,props.cdim);
#ifdef DEBUG
            if(not lref.readOnly()) Error("lref should be readOnly");
#endif
            lref.ApplyTrans(); 
            }
        else 
            { 
            lref = SimpleMatrixRef(Ldat.data(),props.cdim,props.odimL);
#ifdef DEBUG
            if(not lref.readOnly()) Error("lref should be readOnly");
#endif
            }
        }
    else //L not matrix, need to reshape to make lref
        {
        reshape(props.pl,L.is_,L.r_->v,newLdat);
        lref = SimpleMatrixRef(newLdat.data(),props.odimL,props.cdim);
        lref.ApplyTrans(); 
#ifdef DEBUG
            if(lref.readOnly()) Error("lref should not be readOnly");
#endif
        }

    if(R_is_matrix) 
        {
        if(props.contractedR[1]) 
            { 
            rref = SimpleMatrixRef(Rdat.data(),props.odimR,props.cdim);
#ifdef DEBUG
            if(not rref.readOnly()) Error("rref should be readOnly");
#endif
            }
        else                    
            { 
            rref = SimpleMatrixRef(Rdat.data(),props.cdim,props.odimR);
            rref.ApplyTrans(); 
#ifdef DEBUG
            if(not rref.readOnly()) Error("rref should be readOnly");
#endif
            }
        }
    else //R not matrix, need to reshape to make rref
        {
        reshape(props.pr,R.is_,R.r_->v,newRdat);
        rref = SimpleMatrixRef(newRdat.data(),props.odimR,props.cdim);
#ifdef DEBUG
            if(rref.readOnly()) Error("rref should not be readOnly");
#endif
        }

#ifdef COLLECT_PRODSTATS
    /*
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
    */
#endif
    }


//Non-contracting product: Cikj = Aij Bkj (no sum over j)
ITensor& ITensor::
operator/=(const ITensor& other)
    {
    if(type_ == Diag)
        Error("Non-contracting product not yet implemented for type Diag (this)");
    if(other.type_ == Diag)
        Error("Non-contracting product not yet implemented for type Diag (other)");

    if(this == &other)
        {
        ITensor cp_oth(other);
        return operator/=(cp_oth);
        }

    if(scale_.isZero() || other.scale_.isZero())
        {
        scale_ = LogNumber(0);
        return *this;
        }

    if(this->isComplex())
        {
        if(other.isComplex())
            {
            //Both complex
            ITensor rt(*this),
                    it(*this),
                    ro(other),
                    io(other);
            rt.takeRealPart();
            it.takeImagPart();
            ro.takeRealPart();
            io.takeImagPart();

            *this = rt / ro;
            *this -= it / io;

            ITensor ir = rt / io;
            ir += it / ro;

            equalizeScales(ir);

            i_.swap(ir.r_);

            return *this;
            }
        else
            {
            //This complex, other real
            ITensor rr(*this);
            rr.takeRealPart();
            rr /= other;
            ITensor ir(*this);
            ir.takeImagPart();
            ir /= other;
            rr.equalizeScales(ir);
            r_.swap(rr.r_);
            i_.swap(ir.r_);
            scale_ = rr.scale_;
            is_.swap(rr.is_);
            return *this;
            }
        }
    else
    if(other.isComplex())
        {
        //This real, other complex
        ITensor ri = operator/(*this,imagPart(other));
        operator/=(realPart(other));
        equalizeScales(ri);
        i_.swap(ri.r_);
        return *this;
        }

    //------------------------------------------------------------------
    //Handle m==1 Indices: set union

    if(other.is_.rn() == 0)
        {
        scale_ *= other.scale_;
        scale_ *= other.r_->v.at(0);
        for(int j = other.is_.rn()+1; j <= other.r(); ++j)
            {
            const Index& J = other.is_.index(j);
            if(!hasindex(is_,J))
                is_.addindex(J);
            }
        return *this;
        }
    else if(is_.rn() == 0)
        {
        scale_ *= other.scale_;
        scale_ *= r_->v.at(0);
        r_ = other.r_;
        IndexSet<Index> new_is(other.is_);
        for(int j = 1; j <= r(); ++j) 
            {
            const Index& J = is_.index(j);
            if(!hasindex(new_is,J))
                new_is.addindex(J);
            }
        is_.swap(new_is);
        return *this;
        }

    ProductProps props(*this,other);
    SimpleMatrixRef lref, 
                    rref;
    bool L_is_matrix = false,
         R_is_matrix = false;
    vector<Real> newLdat,newRdat;
    toMatrixProd(*this,other,newLdat,newRdat,props,lref,rref,L_is_matrix,R_is_matrix);

    const int ni = lref.Ncols(), 
              nj = lref.Nrows(), 
              nk = rref.Nrows();

    auto nsize = ni*nj*nk;
    auto np = make_shared<ITDat>(nsize,0);
    auto &thisdat = np->v; 
    for(int j = 1; j <= nj; ++j) 
    for(int k = 1; k <= nk; ++k) 
    for(int i = 1; i <= ni; ++i)
        { 
        thisdat[((j-1)*nk+k-1)*ni+i-1] =  rref(k,j) * lref(j,i); 
        }

    r_.swap(np);

    IndexSet<Index> new_index;

    //Handle m!=1 indices

    for(int j = 0; j < is_.rn(); ++j)
        if(!props.contractedL[j+1]) 
            new_index.addindex(this->is_[j]);

    for(int j = 0; j < other.is_.rn(); ++j)
        if(!props.contractedR[j+1]) 
            new_index.addindex(other.is_[j]);

    for(int j = 0; j < is_.rn(); ++j)
        if(props.contractedL[j+1])  
            new_index.addindex(this->is_[j]);

    //Handle m==1 indices

    for(int j = is_.rn(); j < r(); ++j) 
        new_index.addindex(is_[j]);

    for(int j = other.is_.rn()+1; j <= other.r(); ++j)
        {
        const Index& J = other.is_.index(j);
        if(!hasindex(is_,J)) 
            new_index.addindex(J);
        }

    is_.swap(new_index);
    
    scale_ *= other.scale_;

    scaleOutNorm();

    return *this;
    }


void static
directMultiply(const ITensor& L,
               const ITensor& R, 
               ProductProps& props, 
               vector<Real>& newdat,
               IndexSet<Index>& new_index)
    {
    Counter u,  //uncontracted indices
            c;  //contracted indices

    const int zero = 0;

    const int* li[NMAX];
    const int* ri[NMAX];
    for(int n = 0; n < NMAX; ++n)
        {
        li[n] = &zero;
        ri[n] = &zero;
        }

    array<int,NMAX> nl,
                    nr;
    std::fill(nl.begin(),nl.end(),0);
    std::fill(nr.begin(),nr.end(),0);


    const IndexSet<Index>& Lis = L.indices();
    const IndexSet<Index>& Ris = R.indices();

    const int trn = Lis.rn();
    const int orn = Ris.rn();

    for(int j = 0; j < trn; ++j)
        {
        if(!props.contractedL[j+1])
            {
            ++u.rn; //(++u.r);
            u.n[u.rn] = Lis[j].m();
            li[j] = &(u.i[u.rn]);
            new_index.addindex(Lis[j]);
            }
        else
            {
            for(int k = 0; k < orn; ++k)
                {
                if(Lis[j] == Ris[k])
                    {
                    ++c.rn; //(++c.r);
                    c.n[c.rn] = Lis[j].m();
                    li[j] = &(c.i[c.rn]);
                    ri[k] = &(c.i[c.rn]);
                    break;
                    }
                }
            }
        nl[j] = Lis[j].m();
        }

    for(int j = 0; j < orn; ++j)
        {
        if(!props.contractedR[j+1])
            {
            ++u.rn; //(++u.r);
            u.n[u.rn] = Ris[j].m();
            ri[j] = &(u.i[u.rn]);
            new_index.addindex(Ris[j]);
            }
        nr[j] = Ris[j].m();
        }

    newdat.resize(props.odimL*props.odimR);

    const Real* pL = L.datStart();
    const Real* pR = R.datStart();
    Real* pN = newdat.data();

    for(; u.notDone(); ++u)
        {
        Real& val = pN[u.ind];
        val = 0;
        for(c.reset(); c.notDone(); ++c)
            {
            val += pL[((((((((*li[7])*nl[6]+*li[6])*nl[5]+*li[5])*nl[4]+*li[4])
                      *nl[3]+*li[3])*nl[2]+*li[2])*nl[1]+*li[1])*nl[0]+*li[0])]
                 * pR[((((((((*ri[7])*nr[6]+*ri[6])*nr[5]+*ri[5])*nr[4]+*ri[4])
                      *nr[3]+*ri[3])*nr[2]+*ri[2])*nr[1]+*ri[1])*nr[0]+*ri[0])];
            }
        }

    } //directMultiply


void
contractDiagDense(const ITensor& S, const ITensor& T, ITensor& res)
    {
#ifdef DEBUG
    if(!(S.type_ == ITensor::Diag && T.type_ == ITensor::Dense))
        Error("contractDiagDense assumes first argument Diag, second Dense");
#endif

    res.type_ = ITensor::Dense;

    //This is set to true if some of the indices
    //of res come from S.
    //If false, there is an extra loop in the sum
    //tracing over the elements of S.
    bool res_has_Sind = false; 

    //The ti pointers connect
    //the indices of T to either
    //the Counter created below 
    //or the diagonal index of S
    //
    //The ri pointer does the same
    //but for res
    const int zero = 0;
    array<const int*,NMAX+1> ti,
                             ri; 

    for(int n = 0; n <= NMAX; ++n)
        {
        ti[n] = &zero;
        ri[n] = &zero;
        }

    //Index that will loop over 
    //the diagonal elems of S
    int diag_ind = 0;
    const int dsize = S.r_->size();

    //Create a Counter that only loops
    //over the free Indices of T
    Counter tc;

    res.is_.clear();
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
    array<int,NMAX+1> tcon,
                      scon;
    tcon.fill(0);
    scon.fill(0);
    int ncon = 0; //number contracted

    //Analyze contracted Indices
    for(int i = 1; i <= S.r(); ++i)
    for(int j = 1; j <= T.r(); ++j)
        if(S.is_.index(i) == T.is_.index(j))
            {
            scon[i] = j;
            tcon[j] = i;

            ++ncon;
            }

    //Put uncontracted m != 1 Indices
    //of S into res
    for(int i = 1; i <= S.is_.rn(); ++i)
        if(scon[i] == 0)
            {
            res.is_.addindex(S.is_[i-1]);
            alloc_size *= S.is_[i-1].m();
            res_has_Sind = true;

            //Link ri pointer to diagonal of S
            ri[res.is_.r()] = &(diag_ind);
            }

    //Put uncontracted m != 1 Indices
    //of T into res
    for(int i = 1; i <= T.is_.rn(); ++i)
        if(tcon[i] == 0)
            {
            res.is_.addindex(T.is_[i-1]);
            alloc_size *= T.is_[i-1].m();

            //Init appropriate elements
            //of Counter tc
            tc.n[++tc.rn] = T.is_[i-1].m();
            ++tc.r;
            //Link up ti pointer
            ti[i] = &(tc.i[tc.rn]);

            //Link ri pointer to free index of T
            ri[res.is_.r()] = &(tc.i[tc.rn]);
            }
        else
            {
            //If contracted, will
            //be summed with diag of S
            ti[i] = &(diag_ind);
            }

    //Put uncontracted m == 1 Indices
    //of S into res
    for(int i = S.is_.rn()+1; i <= S.r(); ++i)
        if(scon[i] == 0)
            {
            res.is_.addindex(S.is_[i-1]);
            }

    //Put uncontracted m == 1 Indices
    //of T into res
    for(int i = T.is_.rn()+1; i <= T.r(); ++i)
        if(tcon[i] == 0)
            {
            res.is_.addindex(T.is_.index(i));
            }

#ifdef DEBUG
    if(res.is_.r() != (S.r()+T.r() - 2*ncon))
        {
        Print(res.is_);
        printfln("res.is_.r() = %d != (S.r()+T.r()-2*ncon) = %d",
                 res.is_.r(),(S.r()+T.r()-2*ncon));
        Error("Incorrect rank");
        }
#endif

    res.scale_ = S.scale_ * T.scale_;

    //If S has dimension 1
    //it is just a scalar.
    //res may have different m==1 
    //Indices than T, though.
    if(S.is_.rn() == 0)
        {
        res.r_ = T.r_;
        res.i_ = T.i_;
        res *= S.r_->v.at(0);
        return;
        }

    //Allocate a new dat for res if necessary
    if(!res || !res.r_.unique())
        { 
        res.r_ = make_shared<ITDat>(alloc_size,0); 
        }
    else
        {
        res.r_->v.assign(alloc_size,0);
        }

    //Finish initting Counter tc
    for(int k = tc.rn+1; k <= NMAX; ++k)
        {
        tc.n[k] = 1;
        }


    const auto &Tdat = T.r_->v;
    auto &resdat = res.r_->v;

    if(res_has_Sind)
        {
        for(tc.reset(); tc.notDone(); ++tc)
        for(diag_ind = 0; diag_ind < dsize; ++diag_ind)
            {
            resdat[_ind(res.is_,*ri[1],*ri[2],
                                *ri[3],*ri[4],
                                *ri[5],*ri[6],
                                *ri[7],*ri[8])]
             = S.r_->v[diag_ind] 
               * Tdat[_ind(T.is_,*ti[1],*ti[2],
                                 *ti[3],*ti[4],
                                 *ti[5],*ti[6],
                                 *ti[7],*ti[8])];
            }
        }
    else
        {
        for(tc.reset(); tc.notDone(); ++tc)
            {
            Real val = 0;
            for(diag_ind = 0; diag_ind < dsize; ++diag_ind)
                {
                val +=
                S.r_->v[diag_ind] 
                * Tdat[_ind(T.is_,*ti[1],*ti[2],
                                  *ti[3],*ti[4],
                                  *ti[5],*ti[6],
                                  *ti[7],*ti[8])];
                }
            resdat[_ind(res.is_,*ri[1],*ri[2],
                                *ri[3],*ri[4],
                                *ri[5],*ri[6],
                                *ri[7],*ri[8])]
            = val;
            }
        }

    } // contractDiagDense

void
contractDiagDiag(const ITensor& A, const ITensor& B, ITensor& res)
    {
#ifdef DEBUG
    if(!(A.type_ == ITensor::Diag && B.type_ == ITensor::Diag))
        Error("contractDiagDense assumes both arguments Diag");
#endif

    res.is_ = A.is_*B.is_;

    const bool has_common_inds = (res.r() < (A.r()+B.r()));

    res.scale_ = A.scale_*B.scale_;

    if(has_common_inds)
        {
        res.type_ = ITensor::Diag;
        res.allocate();
        auto& rdat = res.r_->v;
        const auto& Adat = A.r_->v;
        const auto& Bdat = B.r_->v;
        rdat = Adat;
        for(int j = 0; j < rdat.size(); ++j)
            {
            rdat[j] *= Bdat[j];
            }
        }
    else //no indices in common
        {
        //TODO: can this be made more efficient, taking advantage of sparsity?
        res = A;
        res.convertToDense();
        res *= B;
        }

    } // contractDiagDiag


ITensor& ITensor::
operator*=(const ITensor& other)
    {
    if(!this->valid() || !other.valid())
        Error("Null ITensor in product");

    if(this == &other)
        {
        ITensor cp_oth(other);
        return operator*=(cp_oth);
        }

    /* Following code suspicious: doesn't modify indices...
    if(scale_.isZero() || other.scale_.isZero())
        {
        scale_ = 0;
        return *this;
        }
        */

    if(this->isComplex())
        {
        if(other.isComplex())
            {
            //Both complex
            ITensor rt(*this),
                    it(*this),
                    ro(other),
                    io(other);
            rt.takeRealPart();
            it.takeImagPart();
            ro.takeRealPart();
            io.takeImagPart();

            *this = rt * ro;
            *this -= it * io;

            ITensor ir = rt * io;
            ir += it * ro;

            equalizeScales(ir);

            i_.swap(ir.r_);

            return *this;
            }
        else
            {
            //This complex, other real
            ITensor ir(imagPart(*this));
            ir *= other;
            takeRealPart();
            operator*=(other);
            equalizeScales(ir);
            i_.swap(ir.r_);
            return *this;
            }
        }
    else
    if(other.isComplex())
        {
        //This real, other complex
        ITensor oi(imagPart(other));
        ITensor cp_this(*this);
        cp_this *= oi;
        operator*=(realPart(other));
        equalizeScales(cp_this);
        i_.swap(cp_this.r_);
        return *this;
        }

    //Handle Diag/Dense cases requiring conversion
    if(type_==Diag && other.type_==Dense)
        {
        ITensor res;
        contractDiagDense(*this,other,res);
        this->swap(res);
        return *this;
        }
    else
    if(type_==Dense && other.type_==Diag)
        {
        ITensor res;
        contractDiagDense(other,*this,res);
        this->swap(res);
        return *this;
        }
    else
    if(type_==Diag && other.type_==Diag)
        {
        ITensor res;
        contractDiagDiag(*this,other,res);
        this->swap(res);
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
        if(!hasindex(other,K))
            new_index1_[++nr1_] = &K;
        }

    for(int j = other.is_.rn(); j < other.r(); ++j)
        {
        const Index& J = other.is_[j];
        if(!hasindex(*this,J))
            new_index1_[++nr1_] = &J;
        }

    //
    //Special cases when one of the tensors
    //has only m==1 indices (effectively a scalar)
    //
    if(other.is_.rn() == 0)
        {
        scale_ *= other.scale_;
        scale_ *= other.r_->v.at(0);
#ifdef DEBUG
        if((is_.rn()+nr1_) > NMAX) 
            {
            println("new r_ would be = ",is_.r());
            Error("ITensor::operator*=: too many uncontracted indices in product (max is 8)");
            }
#endif
        for(int j = 1; j <= is_.rn(); ++j)
            new_index.addindex(is_.index(j));
        //Keep current m!=1 indices, overwrite m==1 indices
        for(int j = 1; j <= nr1_; ++j) 
            new_index.addindex( *(new_index1_[j]) );
        is_.swap(new_index);
        return *this;
        }
    else if(is_.rn() == 0)
        {
        scale_ *= other.scale_;
        scale_ *= r_->v.at(0);
        r_ = other.r_;
#ifdef DEBUG
        if((is_.rn()+nr1_) > NMAX) 
            {
            println("new r_ would be = ",is_.r());
            Error("ITensor::operator*=: too many uncontracted indices in product (max is 8)");
            }
#endif
        for(int j = 1; j <= other.is_.rn(); ++j) 
            new_index.addindex( other.is_.index(j) );
        for(int j = 1; j <= nr1_; ++j) 
            new_index.addindex( *(new_index1_[j]) );
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
        println("new m==1 indices");
        for(int j = 1; j <= nr1_; ++j) println(*(new_index1_.at(j)));
        Error("ITensor::operator*=: too many uncontracted indices in product (max is 8)");
        }
#endif

    //Using a long int here because int was overflowing
    long int complexity = props.odimL;
    complexity *= props.cdim;
    complexity *= props.odimR;
    
    const
    bool do_matrix_multiply = (complexity > 1000);

    SimpleMatrixRef lref, 
                    rref;
    bool L_is_matrix = false,
         R_is_matrix = false;
    vector<Real> Lrs_store,Rrs_store; //storage in case reshape is needed
    toMatrixProd(*this,other,Lrs_store,Rrs_store,props,lref,rref,
                 L_is_matrix,R_is_matrix,do_matrix_multiply);

    if(do_matrix_multiply || (L_is_matrix && R_is_matrix))
        {
        //DO_IF_PS(++Prodstats::stats().c2;)

        //Do the matrix multiplication
        auto nsize = rref.Nrows()*lref.Ncols();
        auto np = make_shared<ITDat>(nsize,0);

        SimpleMatrixRef nref(np->data(),rref.Nrows(),lref.Ncols());
        mult_add(rref,lref,nref,0);

        r_.swap(np);
        
        //Handle m!=1 indices
        for(int j = 0; j < this->is_.rn(); ++j)
            if(!props.contractedL[j+1]) 
                new_index.addindex( is_[j] );
        for(int j = 0; j < other.is_.rn(); ++j)
            if(!props.contractedR[j+1]) 
                new_index.addindex( other.is_[j] );
        }
    else
        {
        auto np = make_shared<ITDat>();
        directMultiply(*this,other,props,np->v,new_index);
        r_.swap(np);
        }


    //Put in m==1 indices
    for(int j = 1; j <= nr1_; ++j) 
        new_index.addindex( *(new_index1_.at(j)) );

    is_.swap(new_index);

    scale_ *= other.scale_;

    scaleOutNorm();

    return *this;
    } //ITensor::operator*=(ITensor)


bool static
checkSameIndOrder(const IndexSet<Index> is1,
                  const IndexSet<Index> is2)
    {
    for(int j = 0; j < is1.rn(); ++j)
    if(is1[j] != is2[j])
        { 
        return false;
        }
    return true;
    }


ITensor& ITensor::
operator+=(const ITensor& other)
    {
    if(!this->valid())
        {
        Error("Calling += on null ITensor");
        return *this;
        }
    if(!other.valid())
        {
        Error("Right-hand side of ITensor += is null");
        return *this;
        }

    if(this == &other) 
        { 
        scale_ *= 2; 
        return *this; 
        }

    if(this->scale_.isZero())
        {
        *this = other;
        return *this;
        }

    //Handle Diag/Dense cases requiring conversion
    if(type_==Dense && other.type_==Diag)
        {
        ITensor cp_o(other);
        cp_o += *this;
        swap(cp_o);
        return *this;
        }
    else
    if(type_==Diag && other.type_==Dense)
        {
        convertToDense();
        }

    const bool bothDiag = (type_==Diag && other.type_==Diag);

    const
    bool same_ind_order = (bothDiag || checkSameIndOrder(is_,other.is_));

    const bool complex_this = this->isComplex();
    const bool complex_other = other.isComplex();
    if(!complex_this && complex_other)
        {
        operator+=(realPart(other));
        if(same_ind_order)
            {
            i_ = other.i_;
            }
        else
            {
            Permutation P(NMAX+1); 
            getperm(is_,other.is_,P);
            allocateImag();
            reshape(P,other.is_,other.i_->v,i_->v);
            }
        if(scale_ != other.scale_)
            {
            soloImag();
            LogNumber nscale = other.scale_/scale_;
            vectormult(i_->v,nscale.real0());
            }
        return *this;
        }
    else
    if(complex_this && !complex_other) 
        {
        ITensor rr(*this);
        rr.takeRealPart();
        rr += other;
        rr.scaleTo(scale_);
        r_.swap(rr.r_);
        return *this;
        }
    else
    if(complex_this && complex_other)
        {
        ITensor r = realPart(*this) + realPart(other);
        ITensor i = imagPart(*this) + imagPart(other);
        i.scaleTo(r.scale_);
        scale_ = r.scale_;
        r_.swap(r.r_);
        i_.swap(i.r_);
        return *this;
        }


    //if(is_ != other.is_)
    //    {
    //    printfln("this ur = %.10f, other.ur = %.10f\n",is_.uniqueReal(),other.is_.uniqueReal());
    //    Print(*this);
    //    Print(other);
    //    Print(this->is_);
    //    printfln("this indexset uniqueReal = %.20f",is_.uniqueReal());
    //    Print(other.is_);
    //    printfln("other indexset uniqueReal = %.20f",other.is_.uniqueReal());
    //    Error("ITensor::operator+=: different Index structure");
    //    }



    Real scalefac = 1;
    if(scale_.magnitudeLessThan(other.scale_)) 
        {
        this->scaleTo(other.scale_); 
        }
    else
        {
        scalefac = (other.scale_/scale_).real();
        }

    solo();

    auto &thisdat = r_->v;
    const auto &othrdat = other.r_->v;

    if(same_ind_order) 
        { 
        vectordaxpy(thisdat,othrdat,scalefac);
        return *this; 
        }

    Permutation P(NMAX+1); 
    getperm(is_,other.is_,P);
    Counter c(other.is_);

    const int* j[NMAX+1];
    for(int k = 1; k <= NMAX; ++k) 
        {
        j[P.dest(k)] = &(c.i[k]);
        }

    int n[NMAX+1];
    for(int k = 1; k <= NMAX; ++k) 
        {
        n[P.dest(k)] = c.n[k];
        }

    if(scalefac == 1)
        {
        for(; c.notDone(); ++c)
            {
            thisdat[(((((((*j[8])*n[7]+*j[7])*n[6]
            +*j[6])*n[5]+*j[5])*n[4]+*j[4])*n[3]
            +*j[3])*n[2]+*j[2])*n[1]+*j[1]]
            += othrdat[c.ind];
            }
        }
    else
        {
        for(; c.notDone(); ++c)
            {
            thisdat[(((((((*j[8])*n[7]+*j[7])*n[6]
            +*j[6])*n[5]+*j[5])*n[4]+*j[4])*n[3]
            +*j[3])*n[2]+*j[2])*n[1]+*j[1]]
            += scalefac * othrdat[c.ind];
            }
        }

    return *this;
    } 

ITensor& ITensor::
operator-=(const ITensor& other)
    {
    if(this == &other) 
        { 
        scale_ = LogNumber(0); 
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
    if(type_ == Diag)
        Error("fromMatrix not implemented for ITensor type Diag");

    DO_IF_DEBUG(if(i1.m() != M.Nrows()) Error("fromMatrix11: wrong number of rows");)
    DO_IF_DEBUG(if(i2.m() != M.Ncols()) Error("fromMatrix11: wrong number of cols");)

    solo();
    scale_ = LogNumber(1);
    is_ = IndexSet<Index>(i1,i2);

    MatrixRef dref; 
    VectorRefNoLink vref(r_->data(),r_->size());
    if(i1 == is_[0])
        {
        vref.TreatAsMatrix(dref,i2.m(),i1.m());
        dref = M.t();
        }
    else
        {
        vref.TreatAsMatrix(dref,i1.m(),i2.m());
        dref = M;
        }
    i_.reset();
    }

void ITensor::
toMatrix11NoScale(const Index& i1, const Index& i2, Matrix& res) const
    {
    if(type_ == Diag)
        Error("toMatrix not implemented for ITensor type Diag");

    if(r() != 2) 
        {
        Print(i1);
        Print(i2);
        for(auto i : this->indices()) printfln("%s: %s",i.id(),i);
        Print(*this);
        Error("toMatrix11: incorrect rank");
        }

    if(this->isComplex())
        {
        Error("toMatrix11 defined only for real ITensor");
        }
#ifdef DEBUG
    if(!hasindex(*this,i1))
        {
        Print(i1);
        Error("ITensor does not have row Index provided.");
        }
    if(!hasindex(*this,i2))
        {
        Print(i2);
        Error("ITensor does not have column Index provided.");
        }
#endif
    res.ReDimension(i1.m(),i2.m());

    MatrixRef dref; 
    VectorRefNoLink vref(const_cast<Real*>(r_->data()),r_->size());
    vref.TreatAsMatrix(dref,is_[1].m(),is_[0].m());
    res = dref.t(i1==is_[0]); 
    }

void ITensor::
toMatrix11(const Index& i1, const Index& i2, Matrix& res) const
    { 
    if(type_ == Diag)
        Error("toMatrix not implemented for ITensor type Diag");

    toMatrix11NoScale(i1,i2,res); 
    res *= scale_.real(); 
    }

void ITensor::
toMatrix12NoScale(const Index& i1, const Index& i2, 
                  const Index& i3, Matrix& res) const
    {
    if(type_ == Diag)
        Error("toMatrix not implemented for ITensor type Diag");

    if(r() != 3) 
        {
        Print(i1);
        Print(i2);
        Print(i3);
        for(auto i : this->indices()) printfln("%s: %s",i.id(),i);
        Print(*this);
        Error("toMatrix12: incorrect rank");
        }
    if(this->isComplex())
        Error("toMatrix12 defined only for real ITensor");
    assert(hasindex(*this,i1));
    assert(hasindex(*this,i2));
    assert(hasindex(*this,i3));

    res.ReDimension(i1.m(),i2.m()*i3.m());

    const array<Index,NMAX> reshuf 
        = {{ i2, i3, i1,    Index::Null(), Index::Null(), 
             Index::Null(), Index::Null(), Index::Null() }};

    Permutation P(NMAX+1); 
    getperm(is_,reshuf,P);

    Vector V;
    reshape(P,is_,r_->v,V);
    res.TreatAsVector() = V;
    }

void ITensor::
toMatrix12(const Index& i1, const Index& i2, 
           const Index& i3, Matrix& res) const
    { 
    if(type_ == Diag)
        Error("toMatrix not implemented for ITensor type Diag");

    toMatrix12NoScale(i1,i2,i3,res); 
    res *= scale_.real(); 
    }

void ITensor::
fromMatrix12(const Index& i1, const Index& i2, 
             const Index& i3, const Matrix& M)
    {
    if(type_ == Diag)
        Error("fromMatrix not implemented for ITensor type Diag");

    ITensor Q(i3,i1,i2);
    auto mref = M.TreatAsVector();
    Q.r_->v.assign(mref.begin(),mref.end());
    *this = Q;
    }

// group i1,i2; i3,i4
void ITensor::
toMatrix22(const Index& i1, const Index& i2, const Index& i3, const Index& i4,Matrix& res) const
    {
    if(type_ == Diag)
        Error("toMatrix not implemented for ITensor type Diag");

    if(r() != 4) Error("toMatrix22: incorrect rank");
    if(this->isComplex())
        Error("toMatrix22 defined only for real ITensor");
    assert(hasindex(*this,i1));
    assert(hasindex(*this,i2));
    assert(hasindex(*this,i3));
    assert(hasindex(*this,i4));
    const 
    int nrow = i1.m() * i2.m(), 
        ncol = i3.m() * i4.m();
    res.ReDimension(nrow,ncol);
    const array<Index,NMAX> reshuf = {{ i3, i4, i1, i2, Index::Null(), Index::Null(), Index::Null(), Index::Null() }};
    Permutation P(NMAX+1); 
    getperm(is_,reshuf,P);
    Vector V; 
    reshape(P,is_,r_->v,V);
    res.TreatAsVector() = V;
    res *= scale_.real0();
    }


void ITensor::
fromMatrix22(const Index& i1, const Index& i2, const Index& i3, const Index& i4, const Matrix& M)
    {
    if(type_ == Diag)
        Error("fromMatrix not implemented for ITensor type Diag");

    if(i3.m()*i4.m() != M.Ncols()) Error("fromMatrix22: wrong number of cols");
    if(i1.m()*i2.m() != M.Nrows()) Error("fromMatrix22: wrong number of rows");
    ITensor Q(i3,i4,i1,i2);
    auto mref = M.TreatAsVector();
    Q.r_->v.assign(mref.begin(),mref.end());
    *this = Q;
    }



/*
void ITensor::toMatrix21(const Index& i1, const Index& i2, const Index& i3, Matrix& res) const
{
    if(r() != 3) Error("toMatrix21: incorrect rank");
    assert(hasindex(*this,i1));
    assert(hasindex(*this,i2));
    res.ReDimension(i1.m()*i2.m(),i3.m());
    const array<Index,NMAX+1> reshuf = {{ Index::Null(), i3, i1, i2, Index::Null(), Index::Null(), Index::Null(), Index::Null(), Index::Null() }};
    Permutation P; 
    getperm(is_,reshuf,P);
    Vector V; reshapeDat(P,V);
    res.TreatAsVector() = V;
    res *= scale_;
}

void ITensor::toMatrix12(const Index& i1, const Index& i2, const Index& i3, Matrix& res) const
{
    if(r() != 3) Error("toMatrix12: incorrect rank");
    assert(hasindex(*this,i1));
    assert(hasindex(*this,i2));
    assert(hasindex(*this,i3));
    res.ReDimension(i1.m(),i2.m()*i3.m());
    const array<Index,NMAX+1> reshuf = {{ Index::Null(), i2, i3, i1, Index::Null(), Index::Null(), Index::Null(), Index::Null(), Index::Null() }};
    Permutation P; 
    getperm(is_,reshuf,P);
    Vector V; reshapeDat(P,V);
    res.TreatAsVector() = V;
    res *= scale_;
}

void ITensor::fromMatrix21(const Index& i1, const Index& i2, const Index& i3, const Matrix& M)
{
    if(r() != 3) Error("fromMatrix21: incorrect rank");
    assert(hasindex(*this,i1));
    assert(hasindex(*this,i2));
    assert(hasindex(*this,i3));
    if(i1.m()*i2.m() != M.Nrows()) Error("fromMatrix21: wrong number of rows");
    if(i3.m() != M.Ncols()) Error("fromMatrix21: wrong number of cols");
    ITensor Q(i3,i1,i2);
    Q.p->v = M.TreatAsVector();
    assignFrom(Q);
}

void ITensor::fromMatrix12(const Index& i1, const Index& i2, const Index& i3, const Matrix& M)
{
    if(r() != 3) Error("fromMatrix12: incorrect rank");
    assert(hasindex(*this,i1) && hasindex(*this,i2) && hasindex(*this,i3));
    if(i1.m() != M.Nrows()) Error("fromMatrix12: wrong number of rows");
    if(i3.m()*i2.m() != M.Ncols()) Error("fromMatrix12: wrong number of cols");
    ITensor Q(i2,i3,i1);
    Q.p->v = M.TreatAsVector();
    assignFrom(Q);
}
*/

void ITensor::
convertToDense()
    {
    ITENSOR_CHECK_NULL
    if(type_ == Diag)
        {
        solo();
        const int dim = is_.dim(); //dense dimension
        shared_ptr<ITDat> oldr = r_;
        allocate(dim);
        const int ds = oldr->size();
        for(int j = 0; j < ds; ++j)
            {
            r_->v[_ind(is_,j,j,j,j,j,j,j,j)] = oldr->v[j];
            }

        if(this->isComplex())
            {
            shared_ptr<ITDat> oldi = i_;
            allocateImag(dim);
            for(int j = 0; j < ds; ++j)
                {
                i_->v[_ind(is_,j,j,j,j,j,j,j,j)] = oldi->v[j];
                }
            }
        }
    type_ = Dense;
    }

ostream& 
operator<<(ostream & s, const ITensor& t)
    {
    s << "ITensor r = " << t.r() << ": ";
    s << t.indices() << "\n";

    s << "  {log(scale)[incl in elems]=" << t.scale().logNum();

    const
    bool iscplx = t.isComplex();
    const 
    bool isdiag = (t.type() == ITensor::Diag);

    if(!t) s << ", dat is null}\n";
    else 
        {
        s << ", L=" << t.indices().dim();

        if(t.scale().isFiniteReal())
            {
            Real nrm = t.norm();
            if(nrm >= 1E-2 && nrm < 1E5)
                {
                s << format(",N=%.2f",nrm);
                }
            else
                {
                s << format(",N=%.1E",nrm);
                }
            }
        else
            {
            s << ",N=too big,scale=" << t.scale();
            }

        s << format("%s%s}\n",(isdiag ? ",D" : ""), (iscplx ? ",C" : ""));

        const bool ff_set = (std::ios::floatfield & s.flags()) != 0;

        if(ff_set || Global::printdat())
            {
            Real scale = 1.0;
            if(t.scale().isFiniteReal()) scale = t.scale().real();
            else s << "  (omitting too large scale factor)" << endl;

            if(t.r() == 0)
                {
                const Real rval = t.r_->v.at(0)*scale;
                if(!iscplx)
                    {
                    s << format("  %.10f\n",rval);
                    }
                else
                    {
                    const Real ival = t.i_->v.at(0)*scale;
                    const char sgn = (ival > 0 ? '+' : '-');
                    s << format("  %.10f%s%.10fi\n",rval,sgn,fabs(ival));
                    }
                return s;
                }

            if(t.type() == ITensor::Diag)
                {
                const int ds = t.indices().front().m();
                for(int j = 1; j <= ds; ++j)
                    {
                    const Real rval = t.r_->v.at(j-1)*scale;
                    if(!iscplx)
                        {
                        if(fabs(rval) > Global::printScale())
                            {
                            s << "  (" << j;
                            for(int n = 2; n <= t.r(); ++n)
                                s << "," << j;
                            s << format(") %.10f\n",rval);
                            }
                        }
                    else
                        {
                        const Real ival = t.i_->v.at(j-1)*scale;
                        if(sqrt(sqr(rval)+sqr(ival)) > Global::printScale())
                            {
                            const char sgn = (ival > 0 ? '+' : '-');
                            s << "  (" << j;
                            for(int n = 2; n <= t.r(); ++n)
                                s << "," << j;
                            s << format(") %.10f%s%.10fi\n",rval,sgn,fabs(ival));
                            }
                        }
                    }
                }
            else
            if(t.type() == ITensor::Dense)
                {
                if(!iscplx)
                    {
                    const Real* pv = t.datStart();
                    Counter c(t.indices());
                    for(; c.notDone(); ++c)
                        {
                        Real val = pv[c.ind]*scale;
                        if(fabs(val) > Global::printScale())
                            { s << "  " << c << format(" %.10f\n",val); }
                        }
                    }
                else //t is complex
                    {
                    const Real* pr = t.datStart();
                    const Real* pi = t.imagDatStart();
                    Counter c(t.indices());
                    for(; c.notDone(); ++c)
                        {
                        Real rval = pr[c.ind]*scale;
                        Real ival = pi[c.ind]*scale;
                        const char sgn = (ival > 0 ? '+' : '-');
                        if(sqrt(sqr(rval)+sqr(ival)) > Global::printScale())
                            { 
                            s << "  " << c << format(" %.10f%s%.10fi\n",rval,sgn,fabs(ival)); 
                            }
                        }
                    }
                }
            }
        }
    return s;
    }



//
// ITDat
//

ITDat::
ITDat() 
    { }

ITDat::
ITDat(size_t size, 
      Real val) 
    : 
    v(size,val)
    { 
    }

ITDat::
ITDat(const VectorRef& v_) 
    : 
    v(v_.begin(),v_.end())
    { }

ITDat::
ITDat(const ITDat& other) 
    : 
    v(other.v)
    { }

void ITDat:: 
read(std::istream& s) 
    { 
    size_t size = 0;
    s.read((char*) &size,sizeof(size));
    v.resize(size);
    s.read((char*) v.data(), sizeof(Real)*size);
    }


void ITDat::
write(std::ostream& s) const 
    { 
    size_t size = v.size();
    s.write((char*) &size, sizeof(size));
    s.write((char*) v.data(), sizeof(Real)*size); 
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
    c_(T.is_),
    P_(NMAX+1)
    { 
    if(!T_) 
        Error("Can't assign to null ITensor");

    if(T.type() == ITensor::Diag)
        Error("commaInit not yet implemented for ITensor type Diag");


    array<Index,NMAX> ii;
    ii.fill(Index::Null());

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
        getperm(T.is_,ii,P_);
        }
    catch(const ITError& e)
        {
        Error("Not enough and/or wrong indices passed to commaInit");
        }

    T_.solo();
    T_.scaleTo(1);
    }

commaInit& commaInit::
operator=(Real r)
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
        { T_.r_->v[c_.ind] = r; ++c_; }
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
        if(x.isComplex() || y.isComplex())
            {
            throw ITError("Must use BraKet, not Dot, for complex ITensors");
            }
        throw ITError("Bad call to Dot, product is not a scalar"); 
        }
    return res.toReal();
    }

Complex 
BraKet(const ITensor& x, const ITensor& y)
    {
    if(x.isComplex())
        {
        ITensor res = dag(x);
        res *= y;
        if(res.r() != 0) 
            {
            throw ITError("Bad call to BraKet, product not a scalar");
            }
        return res.toComplex();
        }
    else
    if(y.isComplex())
        {
        ITensor res = x;
        res *= y;
        if(res.r() != 0) 
            {
            throw ITError("Bad call to BraKet, product not a scalar");
            }
        return res.toComplex();
        }

    return Complex(Dot(x,y),0.);
    }

bool
isZero(const ITensor& T, const Args& args)
    {
    if(T.scale().isZero())
        return true;
    //done with all fast checks
    if(args.getBool("Fast",false)) return false;
    if(T.normNoScale() == 0) return false;
    return false;
    }

}; //namespace itensor
