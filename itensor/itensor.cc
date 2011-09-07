#include "itensor.h"

ostream& operator<<(ostream & s, const ITensor & t)
{
    s << "log(scale)[incl in elems] = " << t.scale().logNum() 
      << ", r = " << t.r() << ": ";

    int i = 1; 
    for(; i < t.r(); ++i) { s << t.index(i) << ", "; } 
    if(t.r() != 0) { s << t.index(i); }

    if(t.is_null()) s << " (dat is null)\n";
    else 
    {
        s << format(" (L=%d,N=%.2f)\n") % t.vec_size() % t.norm();
        if(printdat)
        {
            const Real scale = t.scale_;
            const Vector& v = t.p->v;
            Counter c; t.initCounter(c);
            for(; c != Counter::done; ++c)
            {
                Real val = v(c.ind)*scale;
                if(fabs(val) > 1E-10)
                { s << c << " " << val << "\n"; }
            }
        }
        else { s << "\n"; }
    }
    return s;
}

void ITensor::groupIndices(const array<Index,NMAX+1>& indices, int nind, 
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
        for(int k = 1; k <= r_; ++k) 
        { 
            if(index_[k] == J) 
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
    const int res_rn_ = (nn == 0 ? rn_-nn : rn_-nn+1);

    vector<Index> nindices; nindices.reserve(r_-nind+1);
    Permutation P;
    int kk = 0, kr = 0;
    for(int j = 1; j <= rn_; ++j)
    {
        if(isReplaced[j])
        { 
            P.from_to(j,res_rn_+kr);
            kr += 1;
        }
        else 
        { 
            P.from_to(j,++kk);
            nindices.push_back(index_[j]); 
        }

    }

    nindices.push_back(grouped);

    for(int j = rn_+1; j <= r_; ++j) 
        if(!isReplaced[j]) nindices.push_back(index_[j]);

    if(nn == 0) res = ITensor(nindices,*this);
    else        res = ITensor(nindices,*this,P); 
}

void ITensor::reshapeDat(const Permutation& P, Vector& rdat) const
{
    assert(p != 0);
    const Vector& thisdat = p->v;

    if(P.is_trivial())
	{
        DO_IF_PS(++prodstats.c2;)
        rdat = thisdat;
        return;
	}

    DO_IF_PS(++prodstats.c1;)

    rdat.ReDimension(thisdat.Length());
    rdat = 0;

    const Permutation::int9& ind = P.ind();

    //Make a counter for thisdat
    Counter c; initCounter(c);
    array<int,NMAX+1> n;
    for(int j = 1; j <= rn_; ++j) n[ind[j]] = c.n[j];

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

    if(rn_ == 2 && ind[1] == 2 && ind[2] == 1)
	{
        MatrixRef xref; thisdat.TreatAsMatrix(xref,c.n[2],c.n[1]);
        rdat = Matrix(xref.t()).TreatAsVector();
        return; 
	}
    else if(rn_ == 3)
	{
        DO_IF_PS(int idx = ((ind[1]-1)*3+ind[2]-1)*3+ind[3]; prodstats.perms_of_3[idx] += 1; )
        //Arranged loosely in order of frequency of occurrence
        Bif3(2,1,3) Loop3(i2,i1,i3)
        Bif3(2,3,1) Loop3(i2,i3,i1) //cyclic
        Bif3(3,1,2) Loop3(i3,i1,i2)
        //Bif3(1,3,2) Loop3(i1,i3,i2)
        //Bif3(3,2,1) Loop3(i3,i2,i1)
	}
    else if(rn_ == 4)
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
    else if(rn_ == 5)
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
    else if(rn_ == 6)
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
    for( ; c != Counter::done ; ++c)
    {
        rdat((((((((*j[8]-1)*n[7]+*j[7]-1)*n[6]+*j[6]-1)*n[5]+*j[5]-1)*n[4]+*j[4]-1)*n[3]+*j[3]-1)*n[2]+*j[2]-1)*n[1]+*j[1])
            = thisdat(c.ind);
    }
}

//Converts ITensor dats into MatrixRef's that can be multiplied as rref*lref
//contractedL/R[j] == true if L/R.indexn(j) contracted
void toMatrixProd(const ITensor& L, const ITensor& R, int& nsamen, int& cdim, array<bool,NMAX+1>& contractedL, array<bool,NMAX+1>& contractedR,
                           MatrixRefNoLink& lref, MatrixRefNoLink& rref)
{
    assert(L.p != 0);
    assert(R.p != 0);
    const Vector &Ldat = L.p->v, &Rdat = R.p->v;

    for(int j = 1; j <= NMAX; ++j) { contractedL[j] = contractedR[j] = false; }

    nsamen = 0; //How many m!=1 Indices this and other share
    cdim = 1;   //dimension of contracted inds

    Permutation pl, pr;

    int q = 0; 
    int lcstart = -1, rcstart = -1;
    for(int j = 1; j <= L.rn_; ++j)
    for(int k = 1; k <= R.rn_; ++k)
    if(L.index_[j] == R.index_[k])
    {
        if(lcstart == -1) lcstart = j;
        if(rcstart == -1) rcstart = k;

        ++nsamen;

        contractedL[j] = contractedR[k] = true;

        ++q;
        pl.from_to(j,q);
        pr.from_to(k,q);

        cdim *= L.index_[j].m();
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
        if(!(contractedL[1] || contractedL[L.rn_])) L_is_matrix = false; 
        if(!(contractedR[1] || contractedR[R.rn_])) R_is_matrix = false;
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
            if(!GET(contractedL,Alt.I.dest(j)))
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
            for(int j = L.rn_; j > (L.rn_-nsamen); --j)
            if(!GET(contractedL,Alt.I.dest(j)))
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
            for(int j = 1; j <= L.rn_; ++j)
            if(!contractedL[j]) pl.from_to(j,++q);
            if(L_is_matrix) Error("Calling reshapeDat although L is matrix.");
#ifdef DO_ALT
            L.newAltDat(pl);
            L.reshapeDat(pl,L.lastAlt().v);
            L.lastAlt().v.TreatAsMatrix(lref,odimL,cdim); lref.ApplyTrans();
#else
            Vector lv; L.reshapeDat(pl,lv);
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
            if(!GET(contractedR,Alt.I.dest(j)))
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
            for(int j = R.rn_; j > (R.rn_-nsamen); --j)
            if(!GET(contractedR,Alt.I.dest(j)))
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
            for(int j = 1; j <= R.rn_; ++j)
            if(!contractedR[j]) pr.from_to(j,++q);
            if(R_is_matrix) Error("Calling reshape even though R is matrix.");
#ifdef DO_ALT
            R.newAltDat(pr);
            R.reshapeDat(pr,R.lastAlt().v);
            R.lastAlt().v.TreatAsMatrix(rref,odimR,cdim);
#else
            Vector rv; R.reshapeDat(pr,rv);
            rv.TreatAsMatrix(rref,odimR,cdim);
#endif
            done_with_R = true;
        }
    }

#ifdef COLLECT_PRODSTATS
    if(L.rn_ > R.rn_) ++prodstats.global[make_pair(L.rn_,R.rn_)];
    else ++prodstats.global[make_pair(R.rn_,L.rn_)];
    ++prodstats.total;
    if(L_is_matrix) ++prodstats.did_matrix;
    if(R_is_matrix) ++prodstats.did_matrix;
#endif

}

//Non-contracting product: Cikj = Aij Bkj (no sum over j)
ITensor& ITensor::operator/=(const ITensor& other)
{
    int nrn_ = 0, nr1_ = 0;
    //These hold the indices from other 
    //that will be added to this->index_
    static array<const Index*,NMAX+1> extra_indexn_;
    static array<const Index*,NMAX+1> extra_index1_;

    //------------------------------------------------------------------
    //Handle m==1 Indices: set union
    for(int j = other.rn_+1; j <= other.r_; ++j)
    {
        const Index& J = other.index_[j];
        bool this_has_index = false;
        for(int k = this->rn_+1; k <= this->r_; ++k)
        { if(index_[k] == J) { this_has_index = true; break; } }

        if(!this_has_index) extra_index1_[++nr1_] = &J;
    }

    if(other.rn_ == 0)
    {
        scale_ *= other.scale_;
        scale_ *= other.p->v(1);
        assert(r_+nr1_ <= NMAX);
        for(int j = 0; j < nr1_; ++j) index_[r_+1+j] = *(extra_index1_[j]);
        r_ += nr1_;
        set_unique_Real();
        return *this;
    }
    else if(rn_ == 0)
    {
        scale_ *= other.scale_;
        scale_ *= p->v(1);
        p = other.p;
        rn_ = other.rn_;
        //Move current m==1 indices past rn_
        for(int j = 1; j <= r_; ++j)  index_[rn_+j] = index_[j];
        //Copy other's m!=1 indices
        for(int j = 1; j <= rn_; ++j) index_[j] = other.index_[j];
        r_ += rn_;
        //Get the extra m==1 indices from other
        for(int j = 0; j < nr1_; ++j) index_[r_+1+j] = *(extra_index1_[j]);
        r_ += nr1_;
        set_unique_Real();
        return *this;
    }

    int nsamen, cdim;
    array<bool,NMAX+1> contractedL, contractedR; MatrixRefNoLink lref, rref;
    toMatrixProd(*this,other,nsamen,cdim,contractedL,contractedR,lref,rref);

    if(p->count() != 1) { p = new ITDat(); }
#ifdef DO_ALT
    else { p->alt.clear(); }
#endif
    Vector& thisdat = p->v; 
    
    const int ni = lref.Ncols(), nj = lref.Nrows(), nk = rref.Nrows();
    Matrix L(lref), R(rref);
    thisdat.ReDimension(ni*nj*nk);
    
    for(int j = 1; j <= nj; ++j) for(int k = 1; k <= nk; ++k) for(int i = 1; i <= ni; ++i)
    { thisdat(((j-1)*nk+k-1)*ni+i) =  R(k,j) * L(j,i); }

    //Handle m!=1 indices
    for(int j = 1; j <= other.rn_; ++j)
    { if(!contractedR[j]) extra_indexn_[++nrn_] = &(other.index_[j]); }

    if((r_ + nrn_ + nr1_) > NMAX) Error("ITensor::operator/=: too many indices in product.");

    //Move current m==1 indices to the right
    for(int j = rn_+1; j <= r_; ++j)  index_[nrn_+j] = index_[j];
    //Fill in new m!=1 indices
    for(int j = 0; j < nrn_; ++j) index_[rn_+1+j] = *(extra_indexn_[j]);
    rn_ += nrn_;
    r_ += nrn_;
    //Fill in new m==1 indices
    for(int j = 0; j < nr1_; ++j) index_[r_+1+j] = *(extra_index1_[j]);
    r_ += nr1_;
    
    scale_ *= other.scale_;

    doNormLog();

    set_unique_Real();

    return *this;
}


ITensor& ITensor::operator*=(const ITensor& other)
{
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

    //This holds the m==1 indices that appear in the result
    int nr1_ = 0;
    static array<const Index*,NMAX+1> new_index1_;

    //Handle m==1 Indices

#define USE_GOTO //profile to see which is faster, if any difference

#ifdef USE_GOTO
    for(int k = rn_+1; k <= this->r_; ++k)
    {
        const Index& K = index_[k];
        for(int j = other.rn_+1; j <= other.r_; ++j)
        { if(other.index_[j] == K) { goto skip_this; } }
        new_index1_[++nr1_] = &K;
        skip_this:;
    }
    for(int j = other.rn_+1; j <= other.r_; ++j)
    {
        const Index& J = other.index_[j];
        for(int k = this->rn_+1; k <= this->r_; ++k)
        { if(index_[k] == J) { goto skip_other; } }
        new_index1_[++nr1_] = &J;
        skip_other:;
    }
#else
    for(int k = rn_+1; k <= this->r_; ++k)
    {
        const Index& K = index_[k];
        bool other_has_index = false;
        for(int j = other.rn_+1; j <= other.r_; ++j)
        { if(other.index_[j] == K) { other_has_index = true; break; } }

        if(!other_has_index) new_index1_[++nr1_] = &K;
    }
    for(int j = other.rn_+1; j <= other.r_; ++j)
    {
        const Index& J = other.index_[j];
        bool this_has_index = false;
        for(int k = this->rn_+1; k <= this->r_; ++k)
        { if(index_[k] == J) { this_has_index = true; break; } }

        if(!this_has_index) new_index1_[++nr1_] = &J;
    }
#endif


    if(other.rn_ == 0)
    {
        scale_ *= other.scale_;
        scale_ *= other.p->v(1);
        r_ = rn_ + nr1_;
        assert(r_ <= NMAX);
        //Keep current m!=1 indices, overwrite m==1 indices
        for(int j = 1; j <= nr1_; ++j) index_[rn_+j] = *(new_index1_[j]);
        set_unique_Real();
        return *this;
    }
    else if(rn_ == 0)
    {
        scale_ *= other.scale_;
        scale_ *= p->v(1);
        p = other.p;
        rn_ = other.rn_;
        r_ = rn_ + nr1_;
        assert(r_ <= NMAX);
        //Put in new m==1 Indices (backwards since using pointers)
        for(int j = nr1_; j > 0; --j) index_[rn_+j] = *(new_index1_[j]);
        //Fill in other's m!=1 Indices
        for(int j = 1; j <= rn_; ++j) index_[j] = other.index_[j];
        set_unique_Real();
        return *this;
    }

    int nsamen,cdim;
    array<bool,NMAX+1> contractedL, contractedR; MatrixRefNoLink lref, rref;
    toMatrixProd(*this,other,nsamen,cdim,contractedL,contractedR,lref,rref);

    //Do the matrix multiplication
    if(p->count() != 1) { p = new ITDat(); } 
#ifdef DO_ALT
    else { p->alt.clear(); }
#endif
    p->v.ReDimension(rref.Nrows()*lref.Ncols());
    MatrixRef nref; p->v.TreatAsMatrix(nref,rref.Nrows(),lref.Ncols());
    nref = rref*lref;

    //Create new index_
    static array<Index,NMAX+1> new_index_;

#ifndef NDEBUG
    if((rn_ + other.rn_ - 2*nsamen + nr1_) > NMAX) 
    {
        Print(*this);
        Print(other);
        Print(nsamen);
        cerr << "new m==1 indices\n";
        for(int j = 1; j <= nr1_; ++j) cerr << *(new_index1_.at(j)) << "\n";
        Error("ITensor::operator*=: too many uncontracted indices in product (max is 8)");
    }
#endif

    //Handle m!=1 indices
    int new_rn_ = 0;
    for(int j = 1; j <= this->rn_; ++j)
    { if(!contractedL[j]) new_index_[++new_rn_] = index_[j]; }
    for(int j = 1; j <= other.rn_; ++j)
    { if(!contractedR[j]) new_index_[++new_rn_] = other.index_[j]; }
    rn_ = new_rn_;

    //Put in m==1 indices
    r_ = rn_;
    for(int j = 1; j <= nr1_; ++j) new_index_[++r_] = *(new_index1_.at(j));

    index_.swap(new_index_);

    scale_ *= other.scale_;
    doNormLog();

    set_unique_Real();

    return *this;

} //ITensor::operator*=(ITensor)

void ITensor::reshapeTo(const Permutation& P, ITensor& res) const
{
    res.rn_ = rn_;
    res.r_ = r_;
    res.scale_ = scale_;
    for(int k = 1; k <= r_; ++k) res.index_[P.dest(k)] = index_[k];
    res.ur = ur;
    res.solo();
#ifdef DO_ALT //add a flag to solo to let it know not to even copy alt
    res.p->alt.clear();
#endif
    this->reshapeDat(P,res.p->v);
}

ITensor& ITensor::operator+=(const ITensor& other)
{
    if(this == &other) { scale_ *= 2; return *this; }

    bool complex_this = is_complex();
    bool complex_other = other.is_complex();
    if(!complex_this && complex_other)
    {
        return (*this = (*this * Complex_1) + other);
    }
    if(complex_this && !complex_other) return operator+=(other * Complex_1);

    if(fabs(ur - other.ur) > 1E-12)
    {
        cerr << format("this ur = %.10f, other.ur = %.10f\n")%ur%other.ur;
        Print(*this);
        Print(other);
        Error("ITensor::operator+=: unique Reals don't match (different Index structure).");
    }

    intrusive_ptr<ITDat> curr_p = p;
    const Vector* othrdat = 0;

    if(scale_.magnitudeLessThan(other.scale_)) 
    { 
        this->setScale(other.scale_); 
        solo();
        othrdat = &(other.p->v);
    }
    else
    { 
        //Would be simpler to call other.setScale(this->scale_)
        //but the following prevents other from having to call solo()
        if(p->count() != 1) {  p = new ITDat(other.p->v); }
        else                { p->v = other.p->v; }
        p->v *= (other.scale_/this->scale_);
        othrdat = &(curr_p->v);
    }

    Vector& thisdat = p->v;

#ifdef DO_ALT
    p->alt.clear();
#endif

    bool same_ind_order = true;
    for(int j = 1; j <= rn_; ++j)
    if(index_[j] != other.index_[j])
    { same_ind_order = false; break; }

    if(same_ind_order) { thisdat += *othrdat; return *this; }

    Permutation P; getperm(other.index_,P);
    int *j[NMAX+1];
    Counter c; other.initCounter(c);
    for(int m = 1; m <= NMAX; ++m) j[P.dest(m)] = &(c.i[m]);

    for(; c != Counter::done ; ++c)
    {
        thisdat((((((((*j[8]-1)*c.n[7]+*j[7]-1)*c.n[6]+*j[6]-1)*c.n[5]
        +*j[5]-1)*c.n[4]+*j[4]-1)*c.n[3]+*j[3]-1)*c.n[2]+*j[2]-1)*c.n[1]+*j[1])
        += (*othrdat)(c.ind);
    }

    return *this;
} 

void ITensor::fromMatrix11(const Index& i1, const Index& i2, const Matrix& M)
{
    if(r() != 2) Error("fromMatrix11: incorrect rank");
    assert(hasindex(i1));
    assert(hasindex(i2));
    DO_IF_DEBUG(if(i1.m() != M.Nrows()) Error("fromMatrix11: wrong number of rows");)
    DO_IF_DEBUG(if(i2.m() != M.Ncols()) Error("fromMatrix11: wrong number of cols");)

    solo();
    scale_ = 1;

    MatrixRef dref; p->v.TreatAsMatrix(dref,i2.m(),i1.m());
    if(rn_ == 2)
    { dref = M.t(i1==index(1)); }
    else
    { dref = M.t(); }
}

void ITensor::toMatrix11NoScale(const Index& i1, const Index& i2, Matrix& res) const
{
    if(r() != 2) Error("toMatrix11: incorrect rank");
    assert(hasindex(i1));
    assert(hasindex(i2));
    res.ReDimension(i1.m(),i2.m());

    MatrixRef dref; p->v.TreatAsMatrix(dref,i2.m(),i1.m());
    if(rn_ == 2)
    { res = dref.t(i1==index(1)); }
    else
    { res = dref.t(); }
}
void ITensor::toMatrix11(const Index& i1, const Index& i2, Matrix& res) const
{ toMatrix11NoScale(i1,i2,res); res *= scale_; }

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
    const array<Index,NMAX+1> reshuf = {{ IndNull, i3, i4, i1, i2, IndNull, IndNull, IndNull, IndNull }};
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
    const array<Index,NMAX+1> reshuf = {{ IndNull, i3, i1, i2, IndNull, IndNull, IndNull, IndNull, IndNull }};
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
    const array<Index,NMAX+1> reshuf = {{ IndNull, i2, i3, i1, IndNull, IndNull, IndNull, IndNull, IndNull }};
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


