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
    Scheme scheme_;
    int Minm;
    std::vector<int>  Maxm, Niter;
    std::vector<Real> Cutoff;
    int Nsweep;
    int num_site_center;        // May not be implemented in some cases
    Sweeps(Scheme sch, int nsw, int _minm, int _maxm, Real _cut)
    : scheme_(sch), Minm(_minm), Maxm(nsw+1), Niter(nsw+1,4), Cutoff(nsw+1), Nsweep(nsw), num_site_center(2)
    {
        if(scheme_ == ramp_m)
        {
            for(int s = 1; s <= Nsweep; s++)
            { Cutoff.at(s) = _cut; Maxm.at(s) = (int)(_minm + (s-1.0)/nsw * (_maxm - _minm)); }
        }
        else if(scheme_ == fixed_m || scheme_ == fixed_cutoff)
        {
            for(int s = 1; s <= Nsweep; s++)
            { Cutoff.at(s) = _cut; Maxm.at(s) = _maxm; }
        }
        
        for(int s = 1; s <= min(Nsweep,4); s++)
        { Niter.at(s) = 10 - s; }
    }
    Real cutoff(int sw) const { return Cutoff.at(sw); }
    int minm(int sw) const { return Minm; }
    int maxm(int sw) const { return Maxm.at(sw); }
    int nsweep() const { return Nsweep; }
    int niter(int sw) const { return Niter.at(sw); }
};

#endif // SWEEPS_HEADER_h
