//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SPECTRUM_H
#define __ITENSOR_SPECTRUM_H

#define Cout std::cout
#define Endl std::endl
#define Format boost::format

namespace itensor {

//
// Spectrum
//
// Stores density matrix eigenvalue spectrum following a call to
// svd, denmatDecomp, etc.
//

class Spectrum
    {
    public:

    Spectrum(const OptSet& opts = Global::opts());

    Real 
    truncerr() const { return truncerr_; }
    void 
    truncerr(Real val) { truncerr_ = val; }

    const Vector& 
    eigsKept() const { return eigsKept_; }
    void 
    eigsKept(const Vector& val) { eigsKept_ = val; }

    int
    numEigsKept() const { return eigsKept_.Length(); }

    //
    // Other Methods
    //

    void 
    read(std::istream& s);
    void 
    write(std::ostream& s) const;

    private:

    /////////////////

    Vector eigsKept_;
    Real truncerr_;

    /////////////////

    }; //class Spectrum

inline Spectrum::
Spectrum(const OptSet& opts) 
    :
    truncerr_(NAN)
    { 
    }

void inline Spectrum::
read(std::istream& s)
    {
    s.read((char*)&truncerr_,sizeof(truncerr_));
    eigsKept_.read(s);
    }

void inline Spectrum::
write(std::ostream& s) const
    {
    s.write((char*)&truncerr_,sizeof(truncerr_));
    eigsKept_.write(s);
    }

inline std::ostream& 
operator<<(std::ostream & s, const Spectrum& spec)
    {
    const Vector& eigs = spec.eigsKept();
    const int N = eigs.Length();
    if(N > 0)
        {
        const int max_show = 20;
        int stop = min(N,max_show);
        s << "  Eigs kept: ";
        for(int j = 1; j <= stop; ++j)
            {
            s << Format(eigs(j) > 1E-3 ? ("%.3f") : ("%.3E")) 
                           % eigs(j);
            s << ((j != stop) ? ", " : "\n");
            }
        s << Format("  Trunc. error = %.3E")
             % spec.truncerr()
             << Endl;
        }
    return s;
    }

}; //namespace itensor

#undef Cout
#undef Format
#undef Endl

#endif
