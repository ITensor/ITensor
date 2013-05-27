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

    const Vector& 
    eigsKept() const 
        { return eigsKept_; }
    int
    numEigsKept() const 
        { return eigsKept_.Length(); }

    int 
    maxEigsKept() const;

    Real 
    maxTruncerr() const;

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
    : 
    refNorm_(1)
    { 
    initOpts(opts);
    truncerr_ = NAN;
    }


inline Spectrum::
Spectrum(Real cutoff, int maxm, int minm, Real noise,
         const OptSet& opts = Global::opts())
    : cutoff_(cutoff), 
      maxm_(maxm),
      minm_(minm), 
      noise_(0),
      use_orig_m_(false), 
      doRelCutoff_(doRelCutoff),
      absoluteCutoff_(false), 
      truncate_(true),
      refNorm_(refNorm)
    {  
    initOpts(opts);
    }

void inline Spectrum::
initOpts(const OptSet& opts)
    {
    absoluteCutoff_ = opts.getBool("AbsoluteCutoff",false);
    cutoff_ = opts.getReal("Cutoff",MIN_CUT);
    doRelCutoff_ = opts.getBool("DoRelCutoff",false);
    maxm_ = opts.getInt("Maxm",MAX_M);
    minm_ = opts.getInt("Minm",1);
    noise_ = opts.getReal("Noise",0.);
    truncate_ = opts.getBool("Truncate",true);
    use_orig_m_ = opts.getBool("UseOrigM",false);
    }

int inline Spectrum::
maxEigsKept() const
    {
    int res = -1;
    Foreach(const Vector& eigs,eigsKept_)
        {
        res = max(res,eigs.Length());
        }
    return res;
    }

Real inline Spectrum::
maxTruncerr() const
    {
    Real res = -1;
    Foreach(const Real& te,truncerr_)
        {
        res = max(res,te);
        }
    return res;
    }

void inline Spectrum::
read(std::istream& s)
    {
    s.read((char*)&truncerr_,sizeof(truncerr_));
    s.read((char*)&cutoff_,sizeof(cutoff_));
    s.read((char*)&minm_,sizeof(minm_));
    s.read((char*)&maxm_,sizeof(maxm_));
    s.read((char*)&use_orig_m_,sizeof(use_orig_m_));
    s.read((char*)&showeigs_,sizeof(showeigs_));
    s.read((char*)&doRelCutoff_,sizeof(doRelCutoff_));
    s.read((char*)&absoluteCutoff_,sizeof(absoluteCutoff_));
    s.read((char*)&refNorm_,sizeof(refNorm_));
    //readVec(s,eigsKept_);
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
    s.write((char*)&showeigs_,sizeof(showeigs_));
    s.write((char*)&doRelCutoff_,sizeof(doRelCutoff_));
    s.write((char*)&absoluteCutoff_,sizeof(absoluteCutoff_));
    s.write((char*)&refNorm_,sizeof(refNorm_));
    //writeVec(s,eigsKept_);
    eigsKept_.write(s);
    }

#undef Cout
#undef Format
#undef Endl

#endif
