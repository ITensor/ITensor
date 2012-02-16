#include "itsparse.h"
using namespace std;
using boost::format;
using boost::array;

ITSparse::
ITSparse()
    :
    scale_(0)
    { }

ITSparse::
ITSparse(const Index& i1)
    :
    is_(i1),
    scale_(0)
    { 
    }

ITSparse::
ITSparse(const Index& i1, Real d)
    :
    is_(i1),
    scale_(d)
    { 
    }

ITSparse::
ITSparse(const Index& i1, const Vector& diag)
    :
    diag_(diag),
    is_(i1)
    { 
#ifdef DEBUG
    if(diag_.Length() != i1.m())
        {
        Print(i1);
        Print(diag_.Length());
        Error("Mismatched Index and Vector size");
        }
#endif
    }

ITSparse::
ITSparse(const Index& i1, const Index& i2, Real d)
    :
    is_(i1,i2),
    scale_(d)
    { 
    }

ITSparse::
ITSparse(const Index& i1, const Index& i2, const Vector& diag)
    :
    diag_(diag),
    is_(i1,i2)
    { 
#ifdef DEBUG
    if(diag_.Length() != is_.minM())
        {
        Print(is_);
        Print(is_.minM());
        Print(diag_.Length());
        Error("Vector size must be same as smallest m (> 1)");
        }
#endif
    }

ITSparse::
ITSparse(const Index& i1, const Index& i2, 
              const Index& i3, Real d)
    :
    is_(i1,i2,i3),
    scale_(d)
    {
    }

ITSparse::
ITSparse(const Index& i1, const Index& i2, 
              const Index& i3, const Vector& diag)
    :
    diag_(diag),
    is_(i1,i2,i3)
    {
#ifdef DEBUG
    if(diag_.Length() != is_.minM())
        {
        Print(is_);
        Print(is_.minM());
        Print(diag_.Length());
        Error("Vector size must be same as smallest m (> 1)");
        }
#endif
    }

ITSparse::
ITSparse(const Index& i1, const Index& i2, 
              const Index& i3, const Index& i4, Real d)
    :
    is_(i1,i2,i3,i4),
    scale_(d)
    {
    }

int ITSparse::
diagSize() const
    {
    if(diagAllSame())
        return is_.minM();
    else
        return diag_.Length();
    }

Vector ITSparse::
diag() const
    {
    if(diagAllSame())
        {
        Vector res(diagSize());
        res = scale_.real();
        return res;
        }
    else
        {
        return diag_*scale_.real();
        }
    }

ITSparse& ITSparse::
operator+=(const ITSparse& other)
    {
    if(this == &other)
        {
        operator*=(2);
        return *this;
        }

    if(fabs(is_.uniqueReal() - other.is_.uniqueReal()) > 1E-12)
        {
        cerr << format("this ur = %.10f, other.ur = %.10f\n")%is_.uniqueReal()%other.is_.uniqueReal();
        Print(*this);
        Print(other);
        Error("ITSparse::operator+=: unique Reals don't match (different Index structure).");
        }

    const bool this_allsame = this->diagAllSame();
    const bool othr_allsame = other.diagAllSame();

    if(this_allsame && othr_allsame)
        {
        //In this case, scale_ represents
        //the value of each diag element
        scale_ += other.scale_;
        return *this;
        }

    //Check if this or other is effectively zero
    if(this->scale_.isRealZero())
        {
        *this = other;
        return *this;
        }

    if((other.scale_/scale_).isRealZero()) 
        { 
        return *this; 
        }

    //Determine a scale factor for the sum
    Real scalefac = 1;
    if(scale_.magnitudeLessThan(other.scale_)) 
        {
        this->scaleTo(other.scale_); 
        }
    else
        {
        scalefac = (other.scale_/scale_).real();
        }

    //Already checked both diagAllSame case
    if(this_allsame)
        {
        diag_ = other.diag_;
        diag_ = 1;
        diag_ += scalefac*other.diag_;
        }
    else
    if(othr_allsame)
        {
        diag_ += scalefac;
        }
    else
        {
        diag_ += scalefac*other.diag_;
        }
        
    return *this;
    }

ITSparse& ITSparse::
operator-=(const ITSparse& other)
    {
    if(this == &other) 
        { 
        scale_ = 0; 
        diag_.ReDimension(0); 
        return *this; 
        }
    operator*=(-1); 
    operator+=(other); 
    operator*=(-1);
    return *this;
    }

Real ITSparse::
norm() const
    {
    if(diagAllSame())
        return sqrt(diagSize())*fabs(scale_.real());
    else
        return fabs(Norm(diag_) * scale_.real());
    }

void ITSparse::
scaleOutNorm() const
	{
    if(diag_.Length() == 0) return;

    Real f = Norm(diag_);
    if(fabs(f-1) < 1E-12) return;
    //solo();
    if(f != 0) { diag_ = 1.0/f; scale_ *= f; }
	}

void ITSparse::
scaleTo(LogNumber newscale) const
	{
    if(diag_.Length() == 0) return;

    if(scale_ == newscale) return;
    //solo();
    if(newscale.isRealZero()) 
        { diag_ *= 0; }
    else 
        {
        scale_ /= newscale;
        if(scale_.isRealZero()) diag_ *= 0;
        else diag_ *= scale_.real();
        }
    scale_ = newscale;
	}

void ITSparse::
read(std::istream& s)
    {
    readVec(s,diag_);
    is_.read(s);
    scale_.read(s);
    }

void ITSparse::
write(std::ostream& s) const
    {
    writeVec(s,diag_);
    is_.write(s);
    scale_.write(s);
    }

void
product(const ITSparse& S, const ITensor& T, ITensor& res)
    {
    if(!S.isDiag()) 
        Error("product only implemented for diagonal ITSparses");

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
    int one = 1;
    array<int*,NMAX+1> ti,
                       ri; 
    ti.assign(&one);
    ri.assign(&one);

    //Index that will loop over 
    //the diagonal elems of S
    int diag_ind = 1;
    const int dsize = S.diagSize();

    //Create a Counter that only loops
    //over the free Indices of T
    Counter tc;
    tc.n[0] = 0;

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
    tcon.assign(0);
    scon.assign(0);
    int ncon = 0; //number contracted

    //Analyze contracted Indices
    for(int i = 1; i <= S.r(); ++i)
    for(int j = 1; j <= T.r(); ++j)
        if(S.index(i) == T.index(j))
            {
            scon[i] = j;
            tcon[j] = i;

            ++ncon;
            }

    //Put uncontracted m != 1 Indices
    //of S into res
    for(int i = 1; i <= S.rn(); ++i)
        if(scon[i] == 0)
            {
            res.is_.addindexn(S.index(i));
            alloc_size *= S.m(i);
            res_has_Sind = true;

            //Link ri pointer to diagonal of S
            ri[res.is_.r()] = &(diag_ind);
            }

    //Put uncontracted m != 1 Indices
    //of T into res
    for(int i = 1; i <= T.rn(); ++i)
        if(tcon[i] == 0)
            {
            res.is_.addindexn(T.index(i));
            alloc_size *= T.m(i);

            //Init appropriate elements
            //of Counter tc
            tc.n[++tc.rn_] = T.m(i);
            ++tc.r_;
            //Link up ti pointer
            //cerr << format("Linking ti[%d] to tc.i[%d] (tc.n[%d] = %d)\n") % i % tc.rn_ % tc.rn_ % (tc.n[tc.rn_]);
            ti[i] = &(tc.i[tc.rn_]);

            //Link ri pointer to free index of T
            //cerr << format("Linking ri[%d] to tc.i[%d] (tc.n[%d] = %d)\n") % res.is_.r() % tc.rn_ % tc.rn_ % (tc.n[tc.rn_]);
            ri[res.is_.r()] = &(tc.i[tc.rn_]);
            }
        else
            {
            //If contracted, will
            //be summed with diag of S
            ti[i] = &(diag_ind);
            }

    //Put uncontracted m == 1 Indices
    //of S into res
    for(int i = S.rn()+1; i <= S.r(); ++i)
        if(scon[i] == 0)
            {
            res.is_.addindex1(S.index(i));
            }

    //Put uncontracted m == 1 Indices
    //of T into res
    for(int i = T.rn()+1; i <= T.r(); ++i)
        if(tcon[i] == 0)
            {
            res.is_.addindex1(T.index(i));
            }

#ifdef DEBUG
    if(res.is_.r() != (S.r()+T.r() - 2*ncon))
        {
        Print(res.is_);
        cout << format("res.is_.r() = %d != (S.r()+T.r()-2*ncon) = %d")
            % res.is_.r() % (S.r()+T.r()-2*ncon) << endl;
        Error("Incorrect rank");
        }
#endif

    res.scale_ = S.scale_ * T.scale_;

    //If S has dimension 1
    //it is just a scalar.
    //res may have different m==1 
    //Indices than T, though.
    if(S.rn() == 0)
        {
        res.p = T.p;
        if(!S.diagAllSame())
            res *= S.diag_(1);
        return;
        }

    //Allocate a new dat for res if necessary
    if(res.isNull() || res.p->count() != 1) 
        { 
        res.p = new ITDat(alloc_size); 
        }
    else
        {
        res.p->v.ReDimension(alloc_size);
        res.p->v *= 0;
        }

    //Finish initting Counter tc
    for(int k = tc.rn_+1; k <= NMAX; ++k)
        {
        tc.n[k] = 1;
        }
    tc.reset(1);


    const Vector& Tdat = T.p->v;
    Vector& resdat = res.p->v;

    if(S.diagAllSame())
        {

        if(res_has_Sind)
            {
            //cout << "Doing allSame, res_has_Sind case" << endl;
            //cout << "Case I\n";
            for(; tc.notDone(); ++tc)
            for(diag_ind = 1; diag_ind <= dsize; ++diag_ind)
                {
                resdat(res._ind(*ri[1],*ri[2],
                                *ri[3],*ri[4],
                                *ri[5],*ri[6],
                                *ri[7],*ri[8]))
                 =  Tdat(T._ind(*ti[1],*ti[2],
                                *ti[3],*ti[4],
                                *ti[5],*ti[6],
                                *ti[7],*ti[8]));
                }
            }
        else
            {
            //cout << "Doing allSame, !res_has_Sind case" << endl;
            //cout << "Case II\n";
            for(; tc.notDone(); ++tc)
                {
                Real val = 0;
                for(diag_ind = 1; diag_ind <= dsize; ++diag_ind)
                    {
                    val +=
                    Tdat(T._ind(*ti[1],*ti[2],
                                *ti[3],*ti[4],
                                *ti[5],*ti[6],
                                *ti[7],*ti[8]));
                    }
                resdat(res._ind(*ri[1],*ri[2],
                                *ri[3],*ri[4],
                                *ri[5],*ri[6],
                                *ri[7],*ri[8]))
                = val;
                }
            }
        }
    else //!diag_allsame, use diag_ weights
        {

        if(res_has_Sind)
            {
            //cout << "Case III\n";
            for(; tc.notDone(); ++tc)
            for(diag_ind = 1; diag_ind <= dsize; ++diag_ind)
                {
                resdat(res._ind(*ri[1],*ri[2],
                                *ri[3],*ri[4],
                                *ri[5],*ri[6],
                                *ri[7],*ri[8]))
                 = S.diag_(diag_ind) 
                   * Tdat(T._ind(*ti[1],*ti[2],
                                 *ti[3],*ti[4],
                                 *ti[5],*ti[6],
                                 *ti[7],*ti[8]));
                }
            }
        else
            {
            //cout << "Case IV\n";
            for(; tc.notDone(); ++tc)
                {
                Real val = 0;
                for(diag_ind = 1; diag_ind <= dsize; ++diag_ind)
                    {
                    val +=
                    S.diag_(diag_ind) 
                    * Tdat(T._ind(*ti[1],*ti[2],
                                  *ti[3],*ti[4],
                                  *ti[5],*ti[6],
                                  *ti[7],*ti[8]));
                    }
                resdat(res._ind(*ri[1],*ri[2],
                                *ri[3],*ri[4],
                                *ri[5],*ri[6],
                                *ri[7],*ri[8]))
                = val;
                }
            }
        }


    } // void product(const ITSparse& S, const ITensor& T, ITensor& res)

ostream& 
operator<<(ostream & s, const ITSparse& t)
    {
    s << "log(scale)[incl in elems] = " << t.scale().logNum() 
      << ", r = " << t.r() << ": ";

    s << t.is_;

    if(t.scale_.isFiniteReal())
        {
        Real nrm = t.norm();
        if(nrm >= 1E-2)
            s << format(" (N=%.2f)\n") % nrm;
        else
            s << format(" (N=%.1E)\n") % nrm;
        }
    else
        {
        s << " (N=too big)\n";
        }

    if(Globals::printdat())
        {
        bool finite_scale = t.scale_.isFiniteReal();
        Real scale = (finite_scale ? t.scale_.real() : 1);

        if(t.diag_.Length() == 0)
            {
            if(finite_scale)
                s << format("Diag = %.5f (all same, %d elements)\n")
                     % scale % t.diagSize();
            else
                s << format("Diag = (scale factor too large) (all same, %d elements)\n")
                     % t.diagSize();
            }
        else
            {
            if(!finite_scale) s << "(omitting too large scale factor)\n";

            for(int j = 1; j <= t.diagSize(); ++j)
                {
                Real val = t.diag_(j)*scale;
                if(fabs(val) > Globals::printScale())
                    { s << j << " " << val << "\n"; }
                }
            }
        }
    else 
        {
        s << "\n";
        }

    return s;
    }
