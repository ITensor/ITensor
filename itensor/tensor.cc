#include "tensor.h"
#include <iomanip>
#include "cputime.h"

Real ran1();

IndexVal Index::operator()(int i) const 
{ return IndexVal(*this,i); }

ITensor IndexVal::operator*(const IndexVal& oth) const { ITensor t(*this); return (t *= oth); }

ITensor IndexVal::operator*(Real fac) const { return ITensor(*this,fac); }

inline ITensor operator*(Real fac, const IndexVal& iv) { return ITensor(iv,fac); }

ostream& operator<<(ostream & s, const ITensor & t)
{
    s << "logfac(incl in elems) = " << t.logfac() << ", r = " << t.r() << ": ";
    int i = 0;
    for(i = 1; i <= t.r_n(); ++i) { s << t.indexn(i) << (i != t.r_n() ? ", " : "; "); }
    foreach(const Index& I, t.index1()) { s << I << (i != t.r() ? ", " : ""); ++i; }
    if(t.is_null()) s << " (dat is null)\n";
    else 
    {
        s << format(" (L=%d,N=%.2f)\n") % t.vec_size() % t.norm();
        if(printdat)
        {
            Counter c(t);
            for(; c != Counter::done; ++c)
            {
                assert(c.ind > 0);
                assert(c.ind <= t.Length());
                if(fabs(t.dat()(c.ind)) > 1E-10)
                { s << c << " " << t.dat()(c.ind)*exp(t.logfac()) << "\n"; }
            }
        }
        else s << "\n";
    }

    return s;
}

void ITensor::ReshapeDat(const Permutation& p, Vector& rdat) const
{
    assert(this->p != 0);
    const Vector& thisdat = this->p->v;

    int firstdif = -1;
    for(int j = 1; j <= rn; ++j)
	if(p.ind[j] != j)
    {
	    firstdif = j;
	    break;
    }

    //Permutation is trivial
    if(firstdif == -1)
	{
        DO_IF_PS(++prodstats.c2;)
        rdat = thisdat;
        return;
	}

    DO_IF_PS(++prodstats.c1;)

    rdat.ReDimension(thisdat.Length());

    //Make a counter for thisdat
    Counter c(*this);
    array<int,NMAX+1> n;
    for(int j = 1; j <= rn; ++j) n[p.ind[j]] = c.n[j];

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

#define Bif3(a,b,c) if(p.ind[1] == a && p.ind[2] == b && p.ind[3] == c)

#define Bif4(a,b,c,d) if(p.ind[1] == a && p.ind[2] == b && p.ind[3] == c && p.ind[4] == d)

#define Bif5(a,b,c,d,e) if(p.ind[1] == a && p.ind[2] == b && p.ind[3] == c && p.ind[4]==d && p.ind[5] == e)

#define Bif6(a,b,c,d,e,g) if(p.ind[1] == a && p.ind[2] == b && p.ind[3] == c && p.ind[4]==d && p.ind[5] == e && p.ind[6] == g)

    if(rn == 2 && p.ind[1] == 2 && p.ind[2] == 1)
	{
        MatrixRef xref; thisdat.TreatAsMatrix(xref,c.n[2],c.n[1]);
        rdat = Matrix(xref.t()).TreatAsVector();
        return; 
	}
    else if(rn == 3)
	{
        DO_IF_PS(int idx = ((p.ind[1]-1)*3+p.ind[2]-1)*3+p.ind[3]; prodstats.perms_of_3[idx] += 1; )
        //Arranged loosely in order of frequency of occurrence
        Bif3(2,1,3) Loop3(i2,i1,i3)
        Bif3(2,3,1) Loop3(i2,i3,i1) //cyclic
        Bif3(3,1,2) Loop3(i3,i1,i2)
        //Bif3(1,3,2) Loop3(i1,i3,i2)
        //Bif3(3,2,1) Loop3(i3,i2,i1)
	}
    else if(rn == 4)
	{
        DO_IF_PS(int idx = (((p.ind[1]-1)*4+p.ind[2]-1)*4+p.ind[3]-1)*4+p.ind[4]; prodstats.perms_of_4[idx] += 1; )
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
    else if(rn == 5)
	{
        DO_IF_PS(int idx = ((((p.ind[1]-1)*5+p.ind[2]-1)*5+p.ind[3]-1)*5+p.ind[4]-1)*5+p.ind[5]; prodstats.perms_of_5[idx] += 1; )
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
    else if(rn == 6)
	{
        DO_IF_PS(int idx = (((((p.ind[1]-1)*6+p.ind[2]-1)*6+p.ind[3]-1)*6+p.ind[4]-1)*6+p.ind[5]-1)*6+p.ind[6]; prodstats.perms_of_6[idx] += 1; )
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
    for(int k = 1; k <= NMAX; ++k) { j[p.ind[k]] = &(c.i[k]); }

    //Catch-all loop that works for any tensor
    for( ; c != Counter::done ; ++c)
    {
        rdat((((((((*j[8]-1)*n[7]+*j[7]-1)*n[6]+*j[6]-1)*n[5]+*j[5]-1)*n[4]+*j[4]-1)*n[3]+*j[3]-1)*n[2]+*j[2]-1)*n[1]+*j[1])
            = thisdat(c.ind);
    }
}

//Converts ITensor dats into MatrixRef's that can be multiplied as rref*lref
//contractedL/R[j] == true if L/R.indexn(j) contracted
void toMatrixProd(const ITensor& L, const ITensor& R, array<bool,NMAX+1>& contractedL, array<bool,NMAX+1>& contractedR,
                           MatrixRefNoLink& lref, MatrixRefNoLink& rref)
{
    assert(L.p != 0);
    assert(R.p != 0);
    const Vector &Ldat = L.dat(), &Rdat = R.dat();

    for(int j = 1; j <= NMAX; ++j) { contractedL[j] = contractedR[j] = false; }

    int nsamen = 0, //How many m!=1 Indices this and other share
          cdim = 1;   //dimension of contracted inds

    Permutation pl, pr;

    int q = 0; 
    int lcstart = -1, rcstart = -1;
    for(int j = 1; j <= L.rn; ++j)
    for(int k = 1; k <= R.rn; ++k)
    if(L._indexn[j] == R._indexn[k])
    {
        if(lcstart < 0) lcstart = j;
        if(rcstart < 0) rcstart = k;

        ++nsamen;

        contractedL[j] = contractedR[k] = true;

        ++q;
        pl.ind[j] = q;
        pr.ind[k] = q;

        cdim *= L._indexn[j].m();
    }
    const int odimL = Ldat.Length()/cdim;
    const int odimR = Rdat.Length()/cdim;

    bool L_is_matrix = true, R_is_matrix = true;
    if(nsamen != 0)
    {
        //Check that contracted inds are contiguous
        for(int i = 0; i < nsamen; ++i) 
        {
            if(!contractedL[lcstart+i]) L_is_matrix = false;
            if(!contractedR[rcstart+i]) R_is_matrix = false;
        }
        //Check that contracted inds are all at beginning or end of _indexn
        if(!(contractedL[1] || contractedL[L.rn])) L_is_matrix = false; 
        if(!(contractedR[1] || contractedR[R.rn])) R_is_matrix = false;
    }

    if(L_is_matrix)  
    {
        if(contractedL[1]) 
        { Ldat.TreatAsMatrix(lref,odimL,cdim); lref.ApplyTrans(); }
        else { Ldat.TreatAsMatrix(lref,cdim,odimL); }
    }
    else
    {
        bool done_with_L = false;
#ifdef DO_ALT
        //Not matrix, see if alternate dat is
        foreach(const PDat& Alt, L.p->alt)
        {
            bool front_matrix=true;
            for(int j = 1; j <= nsamen; ++j)
            if(!GET(contractedL,Alt.I.ind[j]))
            { front_matrix = false; break; }

            if(front_matrix) 
            { 
                Alt.v.TreatAsMatrix(lref,odimL,cdim); lref.ApplyTrans();
                done_with_L = true;
                L_is_matrix = true;
                //DO_IF_PS(++prodstats.c3;)
                break;
            }

            bool back_matrix=true;
            for(int j = L.rn; j > (L.rn-nsamen); --j)
            if(!GET(contractedL,Alt.I.ind[j]))
            { back_matrix = false; break; }

            if(back_matrix)
            {
                Alt.v.TreatAsMatrix(lref,cdim,odimL); 
                done_with_L = true;
                L_is_matrix = true;
                DO_IF_PS(++prodstats.c3;)
                break;
            }
        } //for int n
#endif
        //Finish making the permutation (stick non contracted inds on the back)
        if(!done_with_L)
        {
            q = nsamen;
            for(int j = 1; j <= L.rn; ++j)
            if(!contractedL[j]) pl.ind[j] = ++q;
            if(L_is_matrix) Error("Calling ReshapeDat although L is matrix.");
#ifdef DO_ALT
            L.newAltDat(pl);
            L.ReshapeDat(pl,L.lastAlt().v);
            L.lastAlt().v.TreatAsMatrix(lref,odimL,cdim); lref.ApplyTrans();
#else
            Vector lv; L.ReshapeDat(pl,lv);
            lv.TreatAsMatrix(lref,odimL,cdim); lref.ApplyTrans();
#endif
            done_with_L = true;
        }
        assert(done_with_L);
    }

    if(R_is_matrix) 
    {
        if(contractedR[1]) { Rdat.TreatAsMatrix(rref,odimR,cdim); }
        else                    
        { Rdat.TreatAsMatrix(rref,cdim,odimR); rref.ApplyTrans(); }
    }
    else
    {
        bool done_with_R = false;
#ifdef DO_ALT
        //Not matrix, see if alternate dat is
        foreach(const PDat& Alt, R.p->alt)
        {
            bool front_matrix=true;
            for(int j = 1; j <= nsamen; ++j)
            if(!GET(contractedR,Alt.I.ind[j]))
            { front_matrix = false; break; }

            if(front_matrix) 
            { 
                Alt.v.TreatAsMatrix(rref,odimR,cdim); 
                done_with_R = true;
                R_is_matrix = true;
                //DO_IF_PS(++prodstats.c3;)
                break;
            }

            bool back_matrix=true;
            for(int j = R.rn; j > (R.rn-nsamen); --j)
            if(!GET(contractedR,Alt.I.ind[j]))
            { back_matrix = false; break; }

            if(back_matrix)
            {
                Alt.v.TreatAsMatrix(rref,cdim,odimR); rref.ApplyTrans();
                done_with_R = true;
                R_is_matrix = true;
                DO_IF_PS(++prodstats.c3;)
                break;
            }
        } //for int n
#endif
        //Finish making the permutation (stick non contracted inds on the back)
        if(!done_with_R)
        {
            q = nsamen;
            for(int j = 1; j <= R.rn; ++j)
            if(!contractedR[j]) pr.ind[j] = ++q;
            if(R_is_matrix) Error("Calling reshape even though R is matrix.");
#ifdef DO_ALT
            R.newAltDat(pr);
            R.ReshapeDat(pr,R.lastAlt().v);
            R.lastAlt().v.TreatAsMatrix(rref,odimR,cdim);
#else
            Vector rv; R.ReshapeDat(pr,rv);
            rv.TreatAsMatrix(rref,odimR,cdim);
#endif
            done_with_R = true;
        }
    }

#ifdef COLLECT_PRODSTATS
    if(L.rn > R.rn) ++prodstats.global[make_pair(L.rn,R.rn)];
    else ++prodstats.global[make_pair(R.rn,L.rn)];
    ++prodstats.total;
    if(L_is_matrix) ++prodstats.did_matrix;
    if(R_is_matrix) ++prodstats.did_matrix;
#endif

}

//Non-contracting product: Cikj = Aij Bkj (no sum over j)
ITensor& ITensor::operator/=(const ITensor& other)
{
    //------------------------------------------------------------------
    //Handle m==1 Indices
    if(!other._index1.empty()) {
    if(_index1.empty()) _index1 = other._index1;
    else
    {
        //Compute set union
        sort(_index1.begin(),_index1.end());
        sort(other._index1.begin(),other._index1.end());
        list<Index> symdiff(_index1.size()+other._index1.size());
        list<Index>::iterator it =
        set_union(_index1.begin(),_index1.end(),
        other._index1.begin(),other._index1.end(),symdiff.begin());
        _index1.assign(symdiff.begin(),it);
    } }

    if(other.rn == 0)
    {
        _logfac += other._logfac; _neg = (_neg^other._neg);
        set_unique_Real();
        return operator*=(other.p->v(1));
    }

    array<bool,NMAX+1> contractedL, contractedR; MatrixRefNoLink lref, rref;
    toMatrixProd(*this,other,contractedL,contractedR,lref,rref);

    if(p->count() != 1) { p = new Internal::ITDat(0); }
#ifdef DO_ALT
    else { p->alt.clear(); }
#endif
    Vector& thisdat = p->v; 
    
    const int ni = lref.Ncols(), nj = lref.Nrows(), nk = rref.Nrows();
    Matrix L(lref), R(rref);
    thisdat.ReDimension(ni*nj*nk);
    
    //cerr << format("L is %d x %d\n")%L.Nrows()%L.Ncols();

    for(int j = 1; j <= nj; ++j) for(int k = 1; k <= nk; ++k) for(int i = 1; i <= ni; ++i)
    { thisdat(((j-1)*nk+k-1)*ni+i) =  R(k,j) * L(j,i); }

    array<Index,NMAX+1> _new_indexn;
    int new_rn = 0;

    int nsamen = 0;
    for(int j = 1; j <= rn; ++j)
    if(!contractedL[j]) { _new_indexn[++new_rn] = _indexn[j]; }
    else { ++nsamen; } 

    if((rn + other.rn - nsamen) > NMAX) Error("ITensor::operator/=: rn too big in product.");

    for(int j = 1; j <= other.rn; ++j)
    if(!contractedR[j]) { _new_indexn[++new_rn] = other._indexn[j]; }

    for(int j = 1; j <= rn; ++j)
    if(contractedL[j]) { _new_indexn[++new_rn] = _indexn[j]; }

    _indexn = _new_indexn;
    rn = new_rn;
    
    _logfac += other._logfac; _neg = (_neg^other._neg);
    set_unique_Real();

    return *this;
}


ITensor& ITensor::operator*=(const ITensor& other)
{
    /* This code is buggy!
    if(ur == other.ur) //Perform a trace
    {
        Real res = p->v*p->v;
        if(p->count() != 1) { p = new ITDat(0); }
        p->v.ReDimension(1);
        p->v(1) = res;
        rn = 0; _index1.clear();
        _logfac += other._logfac; _neg = (_neg^other._neg);
        set_unique_Real();
        return *this;
    }
    */

    //Complex types are treated as just another index, of type ReIm
    //Multiplication is handled automatically with these simple tensor helpers
    if(findindexn(IndReIm) && other.findindexn(IndReIm) && !other.findindexn(IndReImP)
	   && !other.findindexn(IndReImPP))
	{
        static ITensor primer(IndReIm,IndReImP,1.0);
        static ITensor primerP(IndReIm,IndReImPP,1.0);
        static ITensor prod(IndReIm,IndReImP,IndReImPP);
        static bool first = true;
        if(first)
        {
            IndexVal iv0(IndReIm,1), iv1(IndReImP,1), iv2(IndReImPP,1);
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

    //Handle m==1 Indices
    if(!other._index1.empty()) {
    if(_index1.empty()) _index1 = other._index1;
    else
    {
        //Compute symmetric difference
        sort(_index1.begin(),_index1.end());
        sort(other._index1.begin(),other._index1.end());
        list<Index> symdiff(_index1.size()+other._index1.size());
        list<Index>::iterator it =
        set_symmetric_difference(_index1.begin(),_index1.end(),
        other._index1.begin(),other._index1.end(),symdiff.begin());
        _index1.assign(symdiff.begin(),it);
    } }

    if(other.rn == 0)
    {
        _logfac += other._logfac; _neg = (_neg^other._neg);
        set_unique_Real();
        return operator*=(other.p->v(1));
    }

    array<bool,NMAX+1> contractedL, contractedR; MatrixRefNoLink lref, rref;
    toMatrixProd(*this,other,contractedL,contractedR,lref,rref);

    //Do the matrix multiplication
    if(p->count() != 1) { p = new Internal::ITDat(0); } 
#ifdef DO_ALT
    else { p->alt.clear(); }
#endif
    p->v.ReDimension(rref.Nrows()*lref.Ncols());
    MatrixRef nref; p->v.TreatAsMatrix(nref,rref.Nrows(),lref.Ncols());
    nref = rref*lref;

    //Create new _indexn
    int new_rn = 0;
    for(int j = 1; j <= this->rn; ++j)
    if(!contractedL[j])  _indexn[++new_rn] = _indexn[j];
    for(int j = 1; j <= other.rn; ++j)
    if(!contractedR[j]) _indexn[++new_rn] = other._indexn[j];

    DO_IF_DEBUG(if(new_rn > NMAX) Error("rn too large in product");)

    rn = new_rn;

    _logfac += other._logfac; _neg = (_neg^other._neg);
    donormlog();
    set_unique_Real();

    return *this;

} //ITensor::operator*=(ITensor)

void ITensor::Reshape(const Permutation& p, ITensor& res) const
{
    res.rn = rn;
    res._logfac = _logfac; res._neg = _neg;
    for(int k = 1; k <= rn; ++k) res._indexn[p.ind[k]] = _indexn[k];
    res._index1.assign(_index1.begin(),_index1.end());
    //res.set_unique_Real();
#ifdef DO_ALT
    res.p->alt.clear();
#endif
    this->ReshapeDat(p,res.ncdat());
}

void ITensor::getperm(const ITensor& other, Permutation& P)
{
    if(other.rn != rn)
	{
        cerr << format("this rn = %d, other rn = %d\n")%rn%other.rn;
        this->print("this");
        other.print("other");
        Error("getperm: rn not the same");
	}
    for(int j = 1; j <= rn; ++j)
	{
        bool got_one = false;
        for(int k = 1; k <= rn; ++k)
        if(other._indexn[j] == _indexn[k])
        {
            P.ind[j] = k;
            got_one = true;
            break;
        }
        if(!got_one)
        {
            cerr << "j = " << j << "\n";
            cerr << "this = " << *this << "\n";
            cerr << "other = " << other << "\n";
            Error("getpermBtoA: no matching index");
        }
	}
}

void ITensor::Assign(const ITensor& other)
{
    if(this == &other) return;
    Permutation P; getperm(other,P);
    _logfac = other._logfac; _neg = other._neg;
    //_index1.assign(other._index1.begin(),other._index1.end());
    if(p->count() != 1) { p = new Internal::ITDat(0); }
#ifdef DO_ALT
    else { p->alt.clear(); }
#endif
    other.ReshapeDat(P,p->v);
}

// group i1,i2; i3,i4
void ITensor::toMatrix22(const Index& i1, const Index& i2, const Index& i3, const Index& i4,Matrix& res, Real& lfac) const	// doesn't put in logfac
{
    if(r() != 4) Error("toMatrix22: incorrect rank");
    assert(hasindex(i1));
    assert(hasindex(i2));
    assert(hasindex(i3));
    assert(hasindex(i4));
    int nrow = i1.m() * i2.m(), ncol = i3.m() * i4.m();
    if(nrow != res.Nrows()) Error("toMatrix22: wrong number of rows");
    if(ncol != res.Ncols()) Error("toMatrix22: wrong number of cols");
    ITensor Q(i3,i4,i1,i2);
    Q.Assign(*this);
    res.ReDimension(nrow,ncol);
    res.TreatAsVector() = Q.dat();
    lfac = _logfac;
}
void ITensor::toMatrix22(const Index& i1, const Index& i2, const Index& i3, const Index& i4,Matrix& res) const	// puts in logfac
{ Real lfac; toMatrix22(i1,i2,i3,i4,res,lfac); res *= exp(lfac); }

void ITensor::fromMatrix22(const Index& i1, const Index& i2, const Index& i3, const Index& i4,const Matrix& res)	// doesn't put in logfac
{
    if(r() != 4) Error("fromMatrix22: incorrect rank");
    assert(hasindex(i1));
    assert(hasindex(i2));
    assert(hasindex(i3));
    assert(hasindex(i4));
    if(i3.m()*i4.m() != res.Ncols()) Error("fromMatrix22: wrong number of cols");
    if(i1.m()*i2.m() != res.Nrows()) Error("fromMatrix22: wrong number of rows");
    ITensor Q(i3,i4,i1,i2);
    Q.ncdat() = res.TreatAsVector();
    Assign(Q);
}

void ITensor::toMatrix11(const Index& i1, const Index& i2, Matrix& res, Real& lfac) const	// doesn't put in logfac
{
    if(r() != 2) Error("toMatrix11: incorrect rank");
    assert(hasindex(i1));
    assert(hasindex(i2));
    res.ReDimension(i1.m(),i2.m());
    if(i1 == indexn(2))
    { res.TreatAsVector() = dat(); } 
    else
	{
        ITensor Q(i2,i1);
        Q.Assign(*this);
        res.TreatAsVector() = Q.dat();
	}
    lfac = _logfac;
}
void ITensor::toMatrix11(const Index& i1, const Index& i2, Matrix& res) const	// puts in logfac
{ Real lfac; toMatrix11(i1,i2,res,lfac); res *= exp(lfac); }

void ITensor::fromMatrix11(const Index& i1, const Index& i2, const Matrix& res)	// doesn't put in logfac
{
    if(r() != 2) Error("fromMatrix11: incorrect rank");
    assert(hasindex(i1));
    assert(hasindex(i2));
    if(i1.m() != res.Nrows()) Error("fromMatrix11: wrong number of rows");
    if(i2.m() != res.Ncols()) Error("fromMatrix11: wrong number of cols");
    if(i1 == indexn(2))
    { ncdat() = res.TreatAsVector(); }
    else
	{
        ITensor Q(i2,i1);
        Q.ncdat() = res.TreatAsVector();
        Assign(Q);
	}
}

void ITensor::toMatrix21(const Index& i1, const Index& i2, const Index& i3, Matrix& res, Real& lfac) const // doesn't put in logfac
{
    if(r() != 3) Error("toMatrix21: incorrect rank");
    assert(hasindex(i1));
    assert(hasindex(i2));
    ITensor Q(i3,i1,i2);
    Q.Assign(*this);
    res.ReDimension(i1.m()*i2.m(),i3.m());
    res.TreatAsVector() = Q.dat();
    lfac = _logfac;
}
void ITensor::toMatrix21(const Index& i1, const Index& i2, const Index& i3, Matrix& res) const	// puts in logfac
{ Real lfac; toMatrix21(i1,i2,i3,res,lfac); res *= exp(lfac); }

void ITensor::toMatrix12(const Index& i1, const Index& i2, const Index& i3, Matrix& res, Real& lfac) const	// doesn't put in logfac
{
    if(r() != 3) Error("toMatrix12: incorrect rank");
    assert(hasindex(i1));
    assert(hasindex(i2));
    assert(hasindex(i3));
    ITensor Q(i2,i3,i1);
    Q.Assign(*this);
    res.ReDimension(i1.m(),i2.m()*i3.m());
    res.TreatAsVector() = Q.dat();
    lfac = _logfac;
}
void ITensor::toMatrix12(const Index& i1, const Index& i2, const Index& i3, Matrix& res) const	// puts in logfac
{ Real lfac; toMatrix12(i1,i2,i3,res,lfac); res *= exp(lfac); }

void ITensor::fromMatrix21(const Index& i1, const Index& i2, const Index& i3, const Matrix& res)
{
    if(r() != 3) Error("fromMatrix21: incorrect rank");
    assert(hasindex(i1));
    assert(hasindex(i2));
    assert(hasindex(i3));
    if(i1.m()*i2.m() != res.Nrows()) Error("fromMatrix21: wrong number of rows");
    if(i3.m() != res.Ncols()) Error("fromMatrix21: wrong number of cols");
    ITensor Q(i3,i1,i2);
    Q.ncdat() = res.TreatAsVector();
    Assign(Q);
}

void ITensor::fromMatrix12(const Index& i1, const Index& i2, const Index& i3, const Matrix& res)
{
    if(r() != 3) Error("fromMatrix12: incorrect rank");
    assert(hasindex(i1) && hasindex(i2) && hasindex(i3));
    if(i1.m() != res.Nrows()) Error("fromMatrix12: wrong number of rows");
    if(i3.m()*i2.m() != res.Ncols()) Error("fromMatrix12: wrong number of cols");
    ITensor Q(i2,i3,i1);
    Q.ncdat() = res.TreatAsVector();
    Assign(Q);
}

Index index_in_common(const ITensor& A, const ITensor& B, IndexType t)
{
    foreach(const Index& I, A.indexn())
	if(I.type() == t) if(B.hasindexn(I)) return I;

    foreach(const Index& I, A.index1())
	if(I.type() == t) if(B.hasindex1(I)) return I;

    //cerr << "\n"; A.print(false,"A"); B.print(false,"B");
    //Error("index_in_common: no common Index found");
    return Index();
}

ITensor& ITensor::operator+=(const ITensor& other)
{
    Real dlogfac = other._logfac - _logfac;
    if(dlogfac < -200.0) return *this; // no effect from other

    solo_dosign();
    int sign = (other._neg ? -1 : 1);

    Vector& thisdat = p->v;
    const Vector& othrdat = other.p->v;

    if(this == &other)
    {
        _logfac += log(2); //multiply by 2, without touching p->v
        return *this;
    }

    bool complex_this = is_complex();
    bool complex_other = other.is_complex();
    if(!complex_this && complex_other)
    {
        return (*this = (*this * Complex_1) + other);
    }
    if(complex_this && !complex_other) return operator+=(other * Complex_1);

#ifdef DO_ALT
    p->alt.clear();
#endif

    if(rn != other.rn) 
    {
        cerr << "*this = " << *this << "\n";
        cerr << "other = " << other << "\n";
        Error("ITensor::operator+=: mismatched number of Indices.");
    }

    bool same_ind_order = true;
    for(int j = 1; j <= rn; j++)
    if(_indexn[j] != other._indexn[j])
    {
        same_ind_order = false;
        break;
    }
    if(same_ind_order)
    {
        if(dlogfac < 0.0)	
        {
            thisdat += sign * exp(dlogfac) * othrdat;
            return *this;
        }
        if(dlogfac > 200.0)
        {
            thisdat = othrdat;
            this->_neg = other._neg;
        }
        else
        {
            thisdat *= exp(-dlogfac);
            thisdat += sign * othrdat;
        }
        _logfac = other._logfac;
        return *this;
    }

    Permutation P; getperm(other,P);
    int *j[NMAX+1];
    Counter c(other);
    for(int m = 1; m <= NMAX; ++m) j[P.ind[m]] = &(c.i[m]);
    assert(other.p != 0);
    if(dlogfac < 0.0)	
	{
        Real f = sign * exp(dlogfac);
        for( ; c != Counter::done ; ++c)
            thisdat((((((((*j[8]-1)*this->m(7)+*j[7]-1)*this->m(6)+*j[6]-1)*this->m(5)
            +*j[5]-1)*this->m(4)+*j[4]-1)*this->m(3)+*j[3]-1)*this->m(2)+*j[2]-1)*this->m(1)+*j[1])
            += f * othrdat(c.ind);
        return *this;
	}
    if(dlogfac > 200.0)
    {
        for( ; c != Counter::done ; ++c)
            thisdat((((((((*j[8]-1)*this->m(7)+*j[7]-1)*this->m(6)+*j[6]-1)*this->m(5)
            +*j[5]-1)*this->m(4)+*j[4]-1)*this->m(3)+*j[3]-1)*this->m(2)+*j[2]-1)*this->m(1)+*j[1])
            = othrdat(c.ind);
        this->_neg = other._neg;
    }
    else
	{
        thisdat *= exp(-dlogfac);
        if(other._neg)
        {
        for( ; c != Counter::done ; ++c)
            thisdat((((((((*j[8]-1)*this->m(7)+*j[7]-1)*this->m(6)+*j[6]-1)*this->m(5)
            +*j[5]-1)*this->m(4)+*j[4]-1)*this->m(3)+*j[3]-1)*this->m(2)+*j[2]-1)*this->m(1)+*j[1])
            -= othrdat(c.ind);
        }
        else
        {
        for( ; c != Counter::done ; ++c)
            thisdat((((((((*j[8]-1)*this->m(7)+*j[7]-1)*this->m(6)+*j[6]-1)*this->m(5)
            +*j[5]-1)*this->m(4)+*j[4]-1)*this->m(3)+*j[3]-1)*this->m(2)+*j[2]-1)*this->m(1)+*j[1])
            += othrdat(c.ind);
        }
	}
    _logfac = other._logfac;
    return *this;
} 

void ITensor::SplitReIm(ITensor& re, ITensor& im) const
{
    ITensor cop(*this);
    if(!is_complex())
	{
        re = cop;
        im = cop; im *= 0;
        return;
	}
    cop._logfac = 0;
    ITensor repart(IndReIm), impart(IndReIm);  
    repart.p->v(1) = 1;
    impart.p->v(2) = 1;
    re = cop; re /= repart;
    im = cop; im /= impart;
    assert(fabs(re._logfac) < 1E-10);
    assert(fabs(im._logfac) < 1E-10);
    re._logfac = _logfac; re._neg = _neg;
    im._logfac = _logfac; im._neg = _neg;
}

void ITensor::conj()
{
    if(!is_complex()) return;
    operator/=(ConjTensor);
}

ITensor operator*(const ITensor& t, const Combiner& c)
{
    int j;
    ITensor res;

    if((j = t.findindex1(c.right())) != 0)
    {
        res = t;
        res.removeindex1(j);
        //All of c's left indices must be m==1, so add them all
        foreach(const Index& I, c.left1()) res.addindex1(I);
        return res;
    }
    else if((j = t.findindexn(c.right())) != 0)
	{
        vector<Index> nindices; nindices.reserve(t.r_n()+c.rln()-1);
        for(int i = 1; i < j; ++i)
            nindices.push_back(t.indexn(i));
        foreach(const Index& I, c.leftn())
            nindices.push_back(I);
        for(int i = j+1; i <= t.r_n(); ++i)
            nindices.push_back(t.indexn(i));
        foreach(const Index& I, c.left1()) nindices.push_back(I);
        foreach(const Index& I, t.index1()) nindices.push_back(I);
        const bool do_allocate = false;
        res = ITensor(nindices,do_allocate);
        res.ncdat() = t.dat();
        res.setlogfac(t.logfac());
        return res;
	}

    vector<Index> nindices; nindices.reserve(t.r_n()-c.rln()+1);
    Permutation P;
    for(int i = 1; i <= c.rln(); ++i)
	if((j = t.findindexn(c.leftn(i))) == 0)
    {
	    cerr << "t = " << t << "\n";
	    cerr << "c = " << c << "\n";
        cerr << "Couldn't find 'left' Index " << c.leftn(i) << " in ITensor t." << endl;
	    Error("operator*(ITensor,Combiner): bad Combiner ITensor product");
    }
	else
    {
	    P.ind[j] = t.r_n() - c.rln() + i;
    }

    int k = 1;
    for(int i = 1; i <= t.r_n(); ++i)
	if(c.findindexn(t.indexn(i)) == 0) 
    {
        P.ind[i] = k++;
        nindices.push_back(t.indexn(i));
    }

    nindices.push_back(c.right());

    vector<Index> res_index1 = t.index1();
    foreach(const Index& L, c.left1())
    {
        vector<Index>::iterator it = find(res_index1.begin(),res_index1.end(),L);
        if(it == res_index1.end())
        {
            cout << "t = " << t << "\n"; cout << "c = " << c << "\n";
            cout << "Couldn't find 'left' Index " << L << " in ITensor t." << endl;
            Error("operator*(ITensor,Combiner): bad Combiner ITensor product");
        }
        res_index1.erase(it);
    }

    const bool do_allocate = false;
    res = ITensor(nindices,do_allocate);
    t.ReshapeDat(P,res.ncdat());
    res.addindex1(res_index1);
    res.setlogfac(t.logfac()); res *= (t.neg() ? -1 : 1);

    return res;
}

void Combiner::toITensor(ITensor& res) 
{
    if(right().m() > 16) 
    { cerr << "\n\n" << "WARNING: too large of an m in IQCombiner::toIQTensor(). May be inefficient!\n\n"; }

    ITensor Delta(right(),right().primed(),1);

    //Use the delta tensor to convert this IQCombiner into an IQTensor
    res = (*this) * Delta;
    res.noprimeind(right().primed());
}

