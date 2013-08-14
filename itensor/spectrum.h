//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SPECTRUM_H
#define __ITENSOR_SPECTRUM_H

#define Cout std::cout
#define Endl std::endl
#define Format boost::format

//
// Spectrum
//
// Set parameters for truncating the spectrum of a density-matrix-like
// tensor. Used in conjunction with the methods in svdalgs.h.
//
// On return, stores the eigenvalue spectrum (aka the entanglement spectrum).
//

class Spectrum
    {
    public:

    //
    // Constructors
    //

    Spectrum(const OptSet& opts = Global::opts());

    Spectrum(Real cutoff, 
             int maxm = MAX_M, 
             int minm = 1, 
             Real noise = 0,
             const OptSet& opts = Global::opts());

    //
    // Accessor Methods
    //

    Real 
    cutoff() const { return cutoff_; }
    void 
    cutoff(Real val) { cutoff_ = val; }

    int 
    minm() const { return minm_; }
    void 
    minm(int val) { minm_ = val; }

    int 
    maxm() const { return maxm_; }
    void 
    maxm(int val) { maxm_ = val; }

    Real
    noise() const { return noise_; }
    void
    noise(Real val) { noise_ = val; }

    //Perform the SVD, but do not truncate.
    //Try to keep the original bond dimension 
    //as determined by the shared indices of 
    //the tensors holding the factorized
    //pieces of the original tensor.
    bool 
    useOrigM() const { return use_orig_m_; }
    void 
    useOrigM(bool val) { use_orig_m_ = val; }

    bool
    truncate() const { return truncate_; }
    void
    truncate(bool val) { truncate_ = val; }

    // If doRelCutoff_ is false,
    // refNorm_ defines an overall scale factored
    // out of the denmat before truncating.
    //
    // If doRelCutoff_ is true, refNorm_
    // is determined automatically.
    //
    // (Default is false.)
    //
    bool 
    doRelCutoff() const { return doRelCutoff_; }
    void 
    doRelCutoff(bool val) { doRelCutoff_ = val; }

    // If absoluteCutoff_ == true, do a strict limit
    // on the eigenvalues of the density matrix, 
    // not trying for a truncation error
    bool 
    absoluteCutoff() const { return absoluteCutoff_; }
    void 
    absoluteCutoff(bool val) { absoluteCutoff_ = val; }

    LogNumber 
    refNorm() const { return refNorm_; }
    void 
    refNorm(const LogNumber& val) 
        { 
        if(val.sign() == 0) Error("zero refNorm");
        refNorm_ = val; 
        }

    //
    // Entanglement spectrum analysis methods
    //

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

    void
    initOpts(const OptSet& opts);

    /////////////////
    //
    // Data Members
    //

    Real cutoff_;
    int maxm_;
    int minm_;
    Real noise_;

    bool use_orig_m_,
         doRelCutoff_,
         absoluteCutoff_,
         truncate_;

    LogNumber refNorm_;

    Vector eigsKept_;
    Real truncerr_;

    //
    /////////////////

    }; //class Spectrum

inline Spectrum::
Spectrum(const OptSet& opts) 
    : cutoff_(MIN_CUT),
      maxm_(MAX_M),
      minm_(1),
      noise_(0),
      refNorm_(1),
      truncerr_(NAN)
    { 
    initOpts(opts);
    }


inline Spectrum::
Spectrum(Real cutoff, int maxm, int minm, Real noise,
         const OptSet& opts)
    : cutoff_(cutoff), 
      maxm_(maxm),
      minm_(minm), 
      noise_(0),
      refNorm_(1),
      truncerr_(NAN)
    {  
    initOpts(opts);
    }

void inline Spectrum::
initOpts(const OptSet& opts)
    {
    absoluteCutoff_ = opts.getBool("AbsoluteCutoff",false);
    cutoff_ = opts.getReal("Cutoff",cutoff_);
    doRelCutoff_ = opts.getBool("DoRelCutoff",false);
    maxm_ = opts.getInt("Maxm",maxm_);
    minm_ = opts.getInt("Minm",minm_);
    noise_ = opts.getReal("Noise",noise_);
    truncate_ = opts.getBool("Truncate",true);
    use_orig_m_ = opts.getBool("UseOrigM",false);
    }

void inline Spectrum::
read(std::istream& s)
    {
    s.read((char*)&truncerr_,sizeof(truncerr_));
    s.read((char*)&cutoff_,sizeof(cutoff_));
    s.read((char*)&minm_,sizeof(minm_));
    s.read((char*)&maxm_,sizeof(maxm_));
    s.read((char*)&use_orig_m_,sizeof(use_orig_m_));
    s.read((char*)&doRelCutoff_,sizeof(doRelCutoff_));
    s.read((char*)&absoluteCutoff_,sizeof(absoluteCutoff_));
    s.read((char*)&refNorm_,sizeof(refNorm_));
    eigsKept_.read(s);
    }

void inline Spectrum::
write(std::ostream& s) const
    {
    s.write((char*)&truncerr_,sizeof(truncerr_));
    s.write((char*)&cutoff_,sizeof(cutoff_));
    s.write((char*)&minm_,sizeof(minm_));
    s.write((char*)&maxm_,sizeof(maxm_));
    s.write((char*)&use_orig_m_,sizeof(use_orig_m_));
    s.write((char*)&doRelCutoff_,sizeof(doRelCutoff_));
    s.write((char*)&absoluteCutoff_,sizeof(absoluteCutoff_));
    s.write((char*)&refNorm_,sizeof(refNorm_));
    eigsKept_.write(s);
    }

inline std::ostream& 
operator<<(std::ostream & s, const Spectrum& spec)
    {
    s << Format("  cutoff = %.2E\n  maxm = %d\n  minm = %d\n  noise = %.2E\n")
         % spec.cutoff()
         % spec.maxm()
         % spec.minm()
         % spec.noise();

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

#undef Cout
#undef Format
#undef Endl

#endif
