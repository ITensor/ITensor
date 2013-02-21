//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_QN_H
#define __ITENSOR_QN_H

//
// QN
//
// Quantum number label for IQIndex's.
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
// respect parity but do not conserve Nf
//

class QN
    {
    public:

    QN(int sz=0, int Nf=0);

    QN(int sz, int Nf, int Nfp);

    //sz is measured in units of 1/2
    //(sz==1 is 1/2, sz==2 is 1, etc.)
    int 
    sz() const { return _sz; }

    //Number of particles
    int 
    Nf() const { return _Nf; }

    //Fermion parity
    //Nfp is either 0 or 1
    int 
    Nfp() const { return _Nfp; }

    int 
    sign() const { return (_Nfp == 0 ? +1 : -1); }

    QN 
    operator+(const QN &other) const;

    QN 
    operator-(const QN &other) const;

    QN& 
    operator+=(const QN &other);

    QN& 
    operator-=(const QN &other);

    QN 
    operator-() const;
    
    // Multiplication by Arrows can change QN sign
    QN& 
    operator*=(Arrow dir);
    QN 
    operator*(Arrow dir) const;

    std::string 
    toString() const;

    void 
    print(std::string name = "") const;

    QN(std::istream& s) { read(s); }

    void 
    write(std::ostream& s) const;

    void 
    read(std::istream& s);

    private:

    int _sz, 
        _Nf, 
        _Nfp; //_Nfp stands for fermion number parity, 
              //and tracks whether Nf is even or odd
              //(can't just calculate Nfp from Nf on the fly
              //because we may not be tracking Nf, i.e. Nf==0)
    };

inline QN::
QN(int sz, int Nf)
    :
    _sz(sz),
    _Nf(Nf),
    _Nfp(abs(Nf%2))
    { }

inline QN::
QN(int sz, int Nf, int Nfp) 
    : 
    _sz(sz), 
    _Nf(Nf), 
    _Nfp(abs(Nfp%2))
    { 
#ifdef DEBUG
    if(_Nf != 0 && _Nf%2 != _Nfp)
        Error("Nfp should equal Nf%2");
#endif
    }

QN inline QN::
operator+(const QN &other) const
    { QN res(*this); res+=other; return res; }

QN inline QN::
operator-(const QN &other) const
    { QN res(*this); res-=other; return res; }

QN& inline QN::
operator+=(const QN &other)
    {
    _sz+=other._sz; 
    _Nf+=other._Nf; 
    _Nfp = abs(_Nfp+other._Nfp)%2;
    return *this;
    }

QN& inline QN::
operator-=(const QN &other)
    {
    _sz-=other._sz; 
    _Nf-=other._Nf; 
    _Nfp = abs(_Nfp-other._Nfp)%2;
    return *this;
    }

QN inline QN::
operator-() const  
    { return QN(-_sz,-_Nf,_Nfp); }

QN& inline QN::
operator*=(Arrow dir) 
    { 
    const int i = dir;
    _sz*=i; 
    _Nf*=i; 
    return *this; 
    }

QN inline QN::
operator*(Arrow dir) const 
    { QN res(*this); res*=dir; return res; }

void inline QN::
write(std::ostream& s) const 
    { 
    s.write((char*)&_sz,sizeof(_sz)); 
    s.write((char*)&_Nf,sizeof(_Nf)); 
    s.write((char*)&_Nfp,sizeof(_Nfp)); 
    }

void inline QN::
read(std::istream& s) 
    { 
    s.read((char*)&_sz,sizeof(_sz)); 
    s.read((char*)&_Nf,sizeof(_Nf)); 
    s.read((char*)&_Nfp,sizeof(_Nfp)); 
    }

std::string inline QN::
toString() const
    { return (boost::format("(%+d:%d)")%_sz%_Nf).str(); }

inline std::ostream& 
operator<<(std::ostream &o, const QN &q)
    { 
    return o << boost::format("sz = %d, Nf = %d, fp = %s") 
                % q.sz() % q.Nf() % (q.sign() < 0 ? "-" : "+"); 
    }

void inline QN::
print(std::string name) const
    { 
    std::cout << "\n" << name << " =\n" 
              << *this << std::endl; 
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
operator*(int i,const QN& a)
    { 
    return a*i; 
    }



#endif
