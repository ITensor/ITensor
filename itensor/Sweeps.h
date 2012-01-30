#ifndef __ITENSOR_SWEEPS_HEADER_H
#define __ITENSOR_SWEEPS_HEADER_H
#include "types.h"
#include "input.h"

class Sweeps
    {
    public:
    enum Scheme {ramp_m, fixed_m, fixed_cutoff, exp_m, table};

    //Constructors --------------

    Sweeps(Scheme sch, int nsw, int _minm, int _maxm, Real _cut);

    Sweeps(Scheme sch, int nsw, int nwm, int _minm, int _maxm, Real _cut);

    Sweeps(int nsw, InputGroup& sweep_table);
    
    //Accessor methods ----------

    Scheme 
    scheme() const { return scheme_; }

    int 
    minm(int sw) const { return Minm_.at(sw); }

    int 
    maxm(int sw) const { return Maxm_.at(sw); }

    Real 
    cutoff(int sw) const { return Cutoff_.at(sw); }

    int 
    nsweep() const { return Nsweep_; }

    int 
    nwarm() const { return Nwarm_; }

    int 
    niter(int sw) const { return Niter_.at(sw); }

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
    std::vector<int> Maxm_;
    std::vector<int> Minm_;
    std::vector<int> Niter_;
    std::vector<Real> Cutoff_;
    int Nsweep_, Nwarm_;
    int num_site_center_;        // May not be implemented in some cases
    Real exp_fac_;
    };

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

inline void Sweeps::
init(int _minm, int _maxm, Real _cut)
    {

    Minm_ = std::vector<int>(Nsweep_+1,_minm);
    Maxm_ = std::vector<int>(Nsweep_+1,_maxm);
    Niter_ = std::vector<int>(Nsweep_+1,2);
    Cutoff_ = std::vector<Real>(Nsweep_+1,_cut);

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
    for(int s = 0; s <= min(Nwarm_,3); ++s)
        {
        int ni = Max_Niter-s;
        Niter_.at(1+s) = (ni > 2 ? ni : 2);
        }

    } //Sweeps::init

inline void Sweeps::
tableInit(InputGroup& table)
    {
    if(!table.GotoGroup()) 
        Error("Couldn't find table " + table.name);

    Minm_ = std::vector<int>(Nsweep_+1);
    Maxm_ = std::vector<int>(Nsweep_+1);
    Cutoff_ = std::vector<Real>(Nsweep_+1);
    Niter_ = std::vector<int>(Nsweep_+1);

    table.SkipLine(); //SkipLine so we can have a table key
    for(int i = 1; i <= Nsweep_; i++)
        {
        table.infile.file >> Maxm_[i] >> Minm_[i] >> Cutoff_[i] >> Niter_[i];
        }

    } //Sweeps::tableInit

inline std::ostream&
operator<<(std::ostream& s, const Sweeps& swps)
    {
    s << "Sweeps:\n";
    for(int sw = 1; sw <= swps.nsweep(); ++sw)
        s << boost::format("    Maxm(%d)=%d, Niter(%d)=%d, Cutoff(%d)=%.2E\n")
             % sw % swps.maxm(sw) % sw % swps.niter(sw) % sw % swps.cutoff(sw);
    return s;
    }

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


#endif //__ITENSOR_SWEEPS_HEADER_H
