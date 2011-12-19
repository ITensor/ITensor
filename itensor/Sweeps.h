#ifndef __ITENSOR_SWEEPS_HEADER_H
#define __ITENSOR_SWEEPS_HEADER_H
#include "types.h"

class Sweeps
{
public:
    enum Scheme {ramp_m, fixed_m, fixed_cutoff, exp_m};

    //Constructors --------------

    Sweeps();
    
    Sweeps(Scheme sch);

    Sweeps(Scheme sch, int nsw);

    Sweeps(Scheme sch, int nsw, int _minm, int _maxm, Real _cut);

    Sweeps(Scheme sch, int nsw, int nwm, int _minm, int _maxm, Real _cut);
    
    //Accessor methods ----------

    Scheme 
    scheme() const { return scheme_; }

    void 
    setScheme(Scheme val) { uninit(); scheme_ = val; }

    int 
    minm(int sw) const { return Minm_; }
    void 
    setMinm(int val) { uninit(); Minm_ = val; }

    int 
    maxm(int sw) const { init(); return Maxm_.at(sw); }
    void 
    setMaxm(int val) { uninit(); finalMaxm_ = val; }

    Real 
    cutoff(int sw) const { init(); return Cutoff_.at(sw); }
    void 
    setCutoff(Real val) { uninit(); finalCut_ = val; }

    int 
    nsweep() const { return Nsweep_; }
    void 
    setNsweep(int val) { uninit(); Nsweep_ = val; }

    int 
    nwarm() const { return Nwarm_; }
    void 
    setNwarm(int val) { uninit(); Nwarm_ = val; }

    int 
    niter(int sw) const { init(); return Niter_.at(sw); }

    int 
    maxNiter() const { return Max_Niter_; }
    void 
    setMaxNiter(int val) { uninit(); Max_Niter_ = val; }

    int
    numSiteCenter() const { return num_site_center_; }
    void
    setNumSiteCenter(int val) { num_site_center_ = val; }

private:

    void 
    uninit() { init_ = false; }

    void 
    init() const;

private:
    Scheme scheme_;
    int Minm_, finalMaxm_;
    Real finalCut_;
    mutable std::vector<int>  Maxm_, Niter_;
    mutable std::vector<Real> Cutoff_;
    int Nsweep_, Nwarm_, Max_Niter_;
    int num_site_center_;        // May not be implemented in some cases
    mutable bool init_;
};

inline Sweeps::
Sweeps()
    : scheme_(ramp_m),
      Minm_(1), 
      finalMaxm_(MAX_M), 
      finalCut_(MIN_CUT),
      Nsweep_(0), 
      Nwarm_(0), 
      Max_Niter_(9),
      num_site_center_(2)
    { }

inline Sweeps::
Sweeps(Scheme sch)
    : scheme_(sch),
      Minm_(1), 
      finalMaxm_(MAX_M), 
      finalCut_(MIN_CUT),
      Nsweep_(0), 
      Nwarm_(0), 
      Max_Niter_(9),
      num_site_center_(2)
    { }

inline Sweeps::
Sweeps(Scheme sch, int nsw, int _minm, int _maxm, Real _cut)
    : scheme_(sch), 
      Minm_(_minm), 
      finalMaxm_(_maxm), 
      finalCut_(_cut),
      Maxm_(nsw+1), 
      Niter_(nsw+1,2), 
      Cutoff_(nsw+1), 
      Nsweep_(nsw), 
      Nwarm_(nsw-1), 
      Max_Niter_(9),
      num_site_center_(2), 
      init_(false)
    { }

inline Sweeps::
Sweeps(Scheme sch, int nsw, int nwm, int _minm, int _maxm, Real _cut)
    : scheme_(sch), 
      Minm_(_minm), 
      finalMaxm_(_maxm), 
      finalCut_(_cut),
      Maxm_(nsw+1), 
      Niter_(nsw+1,2), 
      Cutoff_(nsw+1), 
      Nsweep_(nsw), 
      Nwarm_(nwm), 
      Max_Niter_(9),
      num_site_center_(2), 
      init_(false)
    { }

inline Sweeps::
Sweeps(Scheme sch, int nsw)
    : scheme_(sch),
      Minm_(1), 
      finalMaxm_(MAX_M), 
      finalCut_(MIN_CUT),
      Nsweep_(nsw), 
      Nwarm_(0), 
      Max_Niter_(9),
      num_site_center_(2)
    { }

inline void Sweeps::
init() const
    {
    if(init_) return;

    Maxm_.assign(Nsweep_+1,finalMaxm_);
    Niter_.assign(Nsweep_+1,2);
    Cutoff_.assign(Nsweep_+1,finalCut_);

    //Don't want to start with m too low unless requested
    int start_m = (finalMaxm_ < 10 ? Minm_ : 10);

    if(scheme_ == ramp_m)
        {
        //Actual number of warmup sweeps to do
        int act_nwm = min(Nwarm_+1,Nsweep_);
        if(act_nwm > 1) 
        for(int s = 1; s <= act_nwm; ++s)
            Maxm_.at(s) = (int)(start_m + (s-1.0)/(act_nwm-1.0) * (finalMaxm_ - start_m)); 
        }
    else
    if(scheme_ == exp_m)
        {
        int act_nwm = min(Nwarm_+2,Nsweep_);
        if(act_nwm > 1)
        for(int s = 1; s <= act_nwm; ++s)
            {
            int p = (act_nwm-s)/2; //intentional integer division
            Maxm_.at(s) = (int)(start_m + pow(0.5,p) * (finalMaxm_ - start_m)); 
            }
        }
    
    for(int s = 0; s <= min(Nwarm_,3); ++s)
        {
        int ni = Max_Niter_-s;
        Niter_.at(1+s) = (ni > 2 ? ni : 2);
        }

    init_ = true;
    } //Sweeps::init

inline void 
sweepnext(int &l, int &ha, int N, int min_l = 1)
    {
    if(ha == 1)
        {
        if(++l == N) 
            l = N-1, ha = 2;
        return;
        }
    if(l-- == min_l) ha = 3;
    }

inline std::ostream&
operator<<(std::ostream& s, const Sweeps& swps)
    {
    s << "Sweeps:\n";
    for(int sw = 1; sw <= swps.nsweep(); ++sw)
        s << boost::format("    Maxm(%d)=%d, Niter(%d)=%d, Cutoff(%d)=%.2E\n")
             % sw % swps.maxm(sw) % sw % swps.niter(sw) % sw % swps.cutoff(sw);
    return s;
    }


#endif //__ITENSOR_SWEEPS_HEADER_H
