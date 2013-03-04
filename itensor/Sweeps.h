//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SWEEPS_HEADER_H
#define __ITENSOR_SWEEPS_HEADER_H
#include "global.h"
#include "input.h"

template <typename T>
class SweepSetter;

class Sweeps
    {
    public:
    enum Scheme {ramp_m, fixed_m, fixed_cutoff, exp_m, table};

    //Constructors --------------

    Sweeps();

    Sweeps(int nsw, int _minm = 1, int _maxm = 500, Real _cut = 1E-8);

    Sweeps(Scheme sch, int nsw, int _minm, int _maxm, Real _cut);

    Sweeps(Scheme sch, int nsw, int nwarm, int _minm, int _maxm, Real _cut);

    Sweeps(int nsw, InputGroup& sweep_table);
    
    //Accessor methods ----------

    Scheme 
    scheme() const { return scheme_; }

    int 
    minm(int sw) const { return Minm_.at(sw); }
    void 
    setMinm(int sw, int val) { Minm_.at(sw) = val; }

    //Use as sweeps.minm() = 20,20,10; (all remaining set to 10)
    SweepSetter<int> 
    minm();

    int 
    maxm(int sw) const { return Maxm_.at(sw); }
    void 
    setMaxm(int sw, int val) { Maxm_.at(sw) = val; }

    //Use as sweeps.maxm() = 50,50,100,100,500; (all remaining set to 500)
    SweepSetter<int> 
    maxm();

    Real 
    cutoff(int sw) const { return Cutoff_.at(sw); }
    void 
    setCutoff(int sw, Real val) { Cutoff_.at(sw) = val; }

    //Use as sweeps.cutoff() = 1E-8; (cutoff set to 1E-8 for all sweeps)
    //or as sweeps.cutoff() = 1E-4,1E-4,1E-5,1E-5,1E-10; (all remaining set to 1E-10)
    SweepSetter<Real> 
    cutoff();

    Real 
    noise(int sw) const { return Noise_.at(sw); }
    void 
    setNoise(int sw, Real val) { Noise_.at(sw) = val; }
    void 
    setNoise(Real val) { Noise_.assign(Nsweep_+1,val); }

    //Use as sweeps.noise() = 1E-10; (noise set to 1E-10 for all sweeps)
    //or as sweeps.noise() = 1E-8,1E-9,1E-10,0.0; (all remaining set to 0)
    SweepSetter<Real> 
    noise();

    int 
    nsweep() const { return Nsweep_; }
    void 
    nsweep(int val);

    int 
    nwarm() const { return Nwarm_; }

    int 
    niter(int sw) const { return Niter_.at(sw); }
    void 
    setNiter(int sw, int val) { Niter_.at(sw) = val; }
    void 
    setNiter(int val) { Niter_.assign(Nsweep_+1,val); }

    //Use as sweeps.niter() = 5,4,3,2; (all remaining set to 2)
    SweepSetter<int> 
    niter();

    int
    numSiteCenter() const { return num_site_center_; }
    void
    setNumSiteCenter(int val) { num_site_center_ = val; }

    Real 
    expFac() const { return exp_fac_; }
    void 
    setExpFac(Real val) { exp_fac_ = val; }

private:

    void 
    init(int _minm, int _maxm, Real _cut);

    void 
    tableInit(InputGroup& table);

    Scheme scheme_;
    std::vector<int> Maxm_,
                     Minm_,
                     Niter_;
    std::vector<Real> Cutoff_,
                      Noise_;
    int Nsweep_, Nwarm_;
    int num_site_center_;        // May not be implemented in some cases
    Real exp_fac_;

    };

//
// Helper class for Sweeps accessor methods.
// Accumulates a comma separated list of 
// values of type T, storing them in the
// vector v passed to its constructor.
//
template <typename T>
class SweepSetter
    {
    public:

    SweepSetter(std::vector<T>& v)
        :
        started_(false),
        v_(v),
        j_(1)
        { 
        last_val_ = v_[j_];
        }
    
    ~SweepSetter()
        {
        for(; j_ < int(v_.size()); ++j_)
            {
            v_[j_] = last_val_;
            }
        }

    SweepSetter& 
    operator=(T val)
        {
        started_ = true;
        return operator,(val);
        }

    SweepSetter& 
    operator,(T val)
        {
        if(!started_)
            {
            Error("SweepSetter notation is setter() = #, #, #, ... ;");
            }
        if(j_ >= int(v_.size())) return *this;
        v_[j_] = val;
        ++j_;
        last_val_ = val;
        return *this;
        }

    private:
    bool started_;
    std::vector<T>& v_;
    int j_;
    T last_val_;
    };

inline Sweeps::
Sweeps()
    :
    scheme_(table),
    Nsweep_(0),
    Nwarm_(0),
    num_site_center_(2), 
    exp_fac_(0.5)
    {
    }

inline Sweeps::
Sweeps(int nsw, int _minm, int _maxm, Real _cut)
    :
    scheme_(fixed_m),
    Nsweep_(nsw),
    Nwarm_(0),
    num_site_center_(2), 
    exp_fac_(0.5)
    {
    init(_minm,_maxm,_cut);
    }

inline Sweeps::
Sweeps(Scheme sch, int nsw, int _minm, int _maxm, Real _cut)
    : scheme_(sch), 
      Nsweep_(nsw), 
      Nwarm_(nsw-1), 
      num_site_center_(2), 
      exp_fac_(0.5)
    { 
    init(_minm,_maxm,_cut);
    }

inline Sweeps::
Sweeps(Scheme sch, int nsw, int nwm, int _minm, int _maxm, Real _cut)
    : scheme_(sch), 
      Nsweep_(nsw), 
      Nwarm_(nwm), 
      num_site_center_(2), 
      exp_fac_(0.5)
    { 
    init(_minm,_maxm,_cut);
    }

inline Sweeps::
Sweeps(int nsw, InputGroup& sweep_table)
    : scheme_(table),
      Nsweep_(nsw),
      Nwarm_(0), 
      num_site_center_(2),
      exp_fac_(0.5)
    {
    tableInit(sweep_table);
    }

SweepSetter<int> inline Sweeps::
minm() 
    { 
    return SweepSetter<int>(Minm_); 
    }

SweepSetter<int> inline Sweeps::
maxm() 
    { 
    return SweepSetter<int>(Maxm_); 
    }

SweepSetter<Real> inline Sweeps::
cutoff() 
    { 
    return SweepSetter<Real>(Cutoff_); 
    }

SweepSetter<Real> inline Sweeps::
noise() 
    { 
    return SweepSetter<Real>(Noise_); 
    }

SweepSetter<int> inline Sweeps::
niter() 
    { 
    return SweepSetter<int>(Niter_); 
    }

void inline Sweeps::
nsweep(int val)
    { 
    if(val > Nsweep_) 
        Error("Can't use nsweep accessor to increase number of sweeps.");
    Nsweep_ = val; 
    }


void inline Sweeps::
init(int _minm, int _maxm, Real _cut)
    {
    Minm_ = std::vector<int>(Nsweep_+1,_minm);
    Maxm_ = std::vector<int>(Nsweep_+1,_maxm);
    Niter_ = std::vector<int>(Nsweep_+1,2);
    Cutoff_ = std::vector<Real>(Nsweep_+1,_cut);
    Noise_ = std::vector<Real>(Nsweep_+1,0);

    //Don't want to start with m too low unless requested
    int start_m = (_maxm < 10 ? _minm : 10);

    if(scheme_ == ramp_m)
        {
        //Actual number of warmup sweeps to do
        int act_nwm = min(Nwarm_+1,Nsweep_);
        if(act_nwm > 1) 
        for(int s = 1; s <= act_nwm; ++s)
            Maxm_.at(s) = (int)(start_m + (s-1.0)/(act_nwm-1.0) * (_maxm - start_m)); 
        }
    else
    if(scheme_ == exp_m)
        {
        int act_nwm = min(Nwarm_+2,Nsweep_);
        if(act_nwm > 1)
        for(int s = 1; s <= act_nwm; ++s)
            {
            int p = (act_nwm-s)/2; //intentional integer division
            Maxm_.at(s) = (int)(start_m + pow(exp_fac_,p) * (_maxm - start_m)); 
            }
        }
    
    //Set number of Davidson iterations
    const int Max_Niter = 9;
    for(int s = 0; s <= min(Nwarm_,4); ++s)
        {
        int ni = Max_Niter-s;
        Niter_.at(1+s) = (ni > 2 ? ni : 2);
        }

    } //Sweeps::init

void inline Sweeps::
tableInit(InputGroup& table)
    {
    if(!table.GotoGroup()) 
        Error("Couldn't find table " + table.name);

    Minm_ = std::vector<int>(Nsweep_+1);
    Maxm_ = std::vector<int>(Nsweep_+1);
    Cutoff_ = std::vector<Real>(Nsweep_+1);
    Niter_ = std::vector<int>(Nsweep_+1);
    Noise_ = std::vector<Real>(Nsweep_+1);

    table.SkipLine(); //SkipLine so we can have a table key
    for(int i = 1; i <= Nsweep_; i++)
        {
        table.infile.file >> Maxm_[i] >> Minm_[i] >> Cutoff_[i] >> Niter_[i] >> Noise_[i];
        }

    } //Sweeps::tableInit

inline std::ostream&
operator<<(std::ostream& s, const Sweeps& swps)
    {
    s << "Sweeps:\n";
    for(int sw = 1; sw <= swps.nsweep(); ++sw)
        s << boost::format("%d  Maxm=%d, Minm=%d, Cutoff=%.1E, Niter=%d, Noise=%.1E\n")
             % sw % swps.maxm(sw) % swps.minm(sw) % swps.cutoff(sw) %swps.niter(sw) % swps.noise(sw);
    return s;
    }

void inline
sweepnext(int &b, int &ha, int N, int min_b = 1)
    {
    const int inc = (ha==1 ? +1 : -1);
    b += inc;
    if(b == (ha==1 ? N : min_b-1))
        {
        b -= inc;
        ++ha;
        }
    }



#endif //__ITENSOR_SWEEPS_HEADER_H
