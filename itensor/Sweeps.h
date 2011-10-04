#ifndef SWEEPS_HEADER_H
#define SWEEPS_HEADER_H

inline void sweepnext(int &l, int &ha, int N, int min_l = 1)
{
    if(ha == 1)
    {
        if(++l == N) 
            l = N-1, ha = 2;
        return;
    }
    if(l-- == min_l) ha = 3;
}

class Sweeps
{
public:
    enum Scheme {ramp_m, fixed_m, fixed_cutoff};
private:
    Scheme scheme_;
    int Minm, finalMaxm;
    Real finalCut;
    mutable std::vector<int>  Maxm, Niter;
    mutable std::vector<Real> Cutoff;
    int Nsweep;
    int num_site_center;        // May not be implemented in some cases
    mutable bool init_;

public:

    void setScheme(Scheme val) { uninit(); scheme_ = val; }
    void setMinm(int val) { uninit(); Minm = val; }
    void setMaxm(int val) { uninit(); finalMaxm = val; }
    void setCutoff(Real val) { uninit(); finalCut = val; }
    void setNsweep(int val) { uninit(); Nsweep = val; }

    Sweeps() 
    : scheme_(ramp_m),
      Minm(1), Maxm(MAX_M), Cutoff(MIN_CUT),
      Nsweep(0), num_site_center(2)
    { }
    
    Sweeps(Scheme sch)
    : scheme_(sch),
      Minm(1), Maxm(MAX_M), Cutoff(MIN_CUT),
      Nsweep(0), num_site_center(2)
    { }

    Sweeps(Scheme sch, int nsw)
    : scheme_(sch),
      Minm(1), Maxm(MAX_M), Cutoff(MIN_CUT),
      Nsweep(nsw), num_site_center(2)
    { }

    Sweeps(Scheme sch, int nsw, int _minm, int _maxm, Real _cut)
    : scheme_(sch), 
      Minm(_minm), finalMaxm(_maxm), finalCut(_cut),
      Maxm(nsw+1), Niter(nsw+1,4), Cutoff(nsw+1), Nsweep(nsw), 
      num_site_center(2), init_(false)
    { }

    Real cutoff(int sw) const { init(); return Cutoff.at(sw); }
    int minm(int sw) const {  return Minm; }
    int maxm(int sw) const { init(); return Maxm.at(sw); }
    int nsweep() const { return Nsweep; }
    int niter(int sw) const { init(); return Niter.at(sw); }

private:
    void uninit() { init_ = false; }

    void init() const
    {
        if(init_) return;

        Maxm.resize(Nsweep+1);
        Niter.resize(Nsweep+1);
        Cutoff.resize(Nsweep+1);

        if(scheme_ == ramp_m)
        {
            for(int s = 1; s <= Nsweep; s++)
            { Cutoff.at(s) = finalCut; Maxm.at(s) = (int)(Minm + (s-1.0)/Nsweep * (finalMaxm - Minm)); }
        }
        else if(scheme_ == fixed_m || scheme_ == fixed_cutoff)
        {
            for(int s = 1; s <= Nsweep; s++)
            { Cutoff.at(s) = finalCut; Maxm.at(s) = finalMaxm; }
        }
        
        for(int s = 1; s <= Nsweep; s++)
        { 
            if(s <= min(Nsweep,4))
                Niter.at(s) = 10-s;
            else
                Niter.at(s) = 4; 
        }

        init_ = true;
    }
};

#endif // SWEEPS_HEADER_h
