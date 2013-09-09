//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_QN_H
#define __ITENSOR_QN_H

#define DEF_NMAX 2

//
// QN
//
// Quantum number label for IQIndexes.
// 
// For a QN "q", 
//
// q.sz() tracks the Sz eigenvalue mz in
// units of 1/2 (i.e. sz() == 1 is mz = 1/2,
// sz() == -1 is mz = -1/2, sz() == 2 is mz = 1, etc.)
//
// q.Nf() tracks the number of particles
//
// q.Nfp() tracks the fermion parity, which is
// just the particle number mod 2 (always 0 or 1).
// This is useful e.g. for superconductors which
// respect parity but do not conserve Nf.
// If Nmax is set > 2, then Nfp runs from 0 up to Nmax-1;
// this is useful for conserving charge mod Nmax as in a 
// Potts/parafermion model.
//

class QN
    {
    public:

    QN(int sz = 0, int Nf = 0);

    QN(int sz, int Nf, int Nfp);

    //sz is measured in units of 1/2
    //(sz==1 is 1/2, sz==2 is 1, etc.)
    int 
    sz() const { return sz_; }

    //Number of particles
    int 
    Nf() const { return Nf_; }

    //Fermion parity
    //Nfp is either 0 or 1
    //(or can range up to Nmax-1 if Nmax > 2)
    int 
    Nfp() const { return Nfp_; }

    QN& 
    operator+=(const QN &other);

    QN& 
    operator-=(const QN &other);

    QN 
    operator-() const;

    const QN&
    operator+() const { return *this; }
    
    // Multiplication by Arrows can change QN sign
    QN& 
    operator*=(Arrow dir);

    std::string 
    toString() const;

    void 
    write(std::ostream& s) const;

    void 
    read(std::istream& s);

    static int&
    Nmax()
        {
        static int Nmax_ = 2;
        return Nmax_;
        }

    private:

    int sz_, 
        Nf_, 
        Nfp_; //Nfp_ stands for fermion number parity, 
              //and tracks whether Nf is even or odd
              //(can't just calculate Nfp from Nf on the fly
              //because we may not be tracking Nf, i.e. Nf==0)
    };

inline QN::
QN(int sz, int Nf)
    :
    sz_(sz),
    Nf_(Nf),
    Nfp_(abs(Nf%DEF_NMAX))
    { }

inline QN::
QN(int sz, int Nf, int Nfp)
    : 
    sz_(sz), 
    Nf_(Nf), 
    Nfp_(abs(Nfp%Nmax()))
    { 
#ifdef DEBUG
    if(Nf_ != 0 && abs(Nf_%Nmax()) != Nfp_)
        {
        Print(Nmax());
        Error("Nfp should equal abs(Nf%Nmax)");
        }
#endif
    }

QN inline
operator+(QN A, const QN& B) { A += B; return A; }

QN inline
operator-(QN A, const QN& B) { A -= B; return A; }

inline
QN& QN::
operator+=(const QN &other)
    {
    sz_ += other.sz_; 
    Nf_ += other.Nf_; 
    Nfp_ = abs(Nfp_+other.Nfp_)%Nmax();
    return *this;
    }

inline
QN& QN::
operator-=(const QN &other)
    {
    sz_ -= other.sz_; 
    Nf_ -= other.Nf_; 
    Nfp_ = abs(Nfp_-other.Nfp_)%Nmax();
    return *this;
    }

QN inline QN::
operator-() const 
    { 
    return QN(-sz_,-Nf_,(Nmax()-Nfp_)%Nmax()); 
    }

inline
QN& QN::
operator*=(Arrow dir) 
    { 
    const int i = dir;
    sz_*=i; 
    Nf_*=i; 
    if(i == -1) Nfp_ = (Nmax()-Nfp_)%Nmax();
    return *this; 
    }

QN inline
operator*(QN q, Arrow dir) { q *= dir; return q; }

void inline QN::
write(std::ostream& s) const 
    { 
    s.write((char*)&sz_,sizeof(sz_)); 
    s.write((char*)&Nf_,sizeof(Nf_)); 
    s.write((char*)&Nfp_,sizeof(Nfp_)); 
    }

void inline QN::
read(std::istream& s) 
    { 
    s.read((char*)&sz_,sizeof(sz_)); 
    s.read((char*)&Nf_,sizeof(Nf_)); 
    s.read((char*)&Nfp_,sizeof(Nfp_)); 
    }

std::string inline QN::
toString() const
    { return (boost::format("(%+d:%d)")%sz_%Nf_).str(); }

inline std::ostream& 
operator<<(std::ostream &o, const QN &q)
    { 
    return o << boost::format("(sz=%d, Nf=%d, p=%d)") 
                % q.sz() 
                % q.Nf()
                % q.Nfp();
    }

bool inline
operator==(const QN &a,const QN &b)
    { 
    return a.sz() == b.sz() 
        && a.Nf() == b.Nf() 
        && a.Nfp() == b.Nfp(); 
    }

bool inline
operator!=(const QN &a,const QN &b)
    { 
    return a.sz() != b.sz() 
        || a.Nf() != b.Nf() 
        || a.Nfp() != b.Nfp(); 
    }

bool inline
operator<(const QN &a,const QN &b)
    { 
    return a.sz() < b.sz() 
        || (a.sz() == b.sz() && a.Nf() < b.Nf()) 
        || (a.sz() == b.sz() && a.Nf() == b.Nf() && a.Nfp() < b.Nfp()); 
    }

QN inline
operator*(Arrow dir, QN q)
    { 
    q *= dir;
    return q;
    }


#undef DEF_NMAX

#endif
