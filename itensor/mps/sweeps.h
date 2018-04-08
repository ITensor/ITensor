//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SWEEPS_HEADER_H
#define __ITENSOR_SWEEPS_HEADER_H
#include <functional>
#include "itensor/util/input.h"
#include "itensor/global.h"
#include "itensor/util/readwrite.h"

namespace itensor {

template <typename T>
class SweepSetter;

//
// To use the InputGroup / table constructor,
// the format required in the input file is:
//
// sweep_table_name
//      {
//      maxm   minm  cutoff  niter  noise
//      20     20    1E-8    4      1E-8
//      40     20    1E-8    3      1E-9
//      80     20    1E-10   2      1E-10
//      160    20    1E-12   2      0
//      240    20    1E-12   2      0
//      }
//
class Sweeps
    {
    public:

    //Constructors --------------

    Sweeps();

    Sweeps(int nsweeps, 
           int minm = 1, 
           int maxm = 500, 
           Real cutoff = 1E-8,
           Real noise = 0.);

    Sweeps(Args const& args);

    Sweeps(int nsweeps, InputGroup& sweep_table);
    
    //Accessor methods ----------

    int 
    nsweep() const { return nsweep_; }
    void 
    nsweep(int val);

    int 
    size() const { return nsweep_; }

    int 
    minm(int sw) const { return minm_.at(sw); }
    void 
    setminm(int sw, int val) { minm_.at(sw) = val; }

    //Use as sweeps.minm() = 20,20,10; (all remaining set to 10)
    SweepSetter<int> 
    minm();

    int 
    maxm(int sw) const { return maxm_.at(sw); }
    void 
    setmaxm(int sw, int val) { maxm_.at(sw) = val; }

    //Use as sweeps.maxm() = 50,50,100,100,500; (all remaining set to 500)
    SweepSetter<int> 
    maxm();

    Real 
    cutoff(int sw) const { return cutoff_.at(sw); }
    void 
    setcutoff(int sw, Real val) { cutoff_.at(sw) = val; }

    //Use as sweeps.cutoff() = 1E-8; (cutoff set to 1E-8 for all sweeps)
    //or as sweeps.cutoff() = 1E-4,1E-4,1E-5,1E-5,1E-10; (all remaining set to 1E-10)
    SweepSetter<Real> 
    cutoff();

    Real 
    noise(int sw) const { return noise_.at(sw); }
    void 
    setnoise(int sw, Real val) { noise_.at(sw) = val; }
    void 
    setnoise(Real val) { noise_.assign(nsweep_+1,val); }

    //Use as sweeps.noise() = 1E-10; (noise set to 1E-10 for all sweeps)
    //or as sweeps.noise() = 1E-8,1E-9,1E-10,0.0; (all remaining set to 0)
    SweepSetter<Real> 
    noise();

    int 
    niter(int sw) const { return niter_.at(sw); }
    void 
    setniter(int sw, int val) { niter_.at(sw) = val; }
    void 
    setniter(int val) { niter_.assign(nsweep_+1,val); }

    //Use as sweeps.niter() = 5,4,3,2; (all remaining set to 2)
    SweepSetter<int> 
    niter();

    void
    read(std::istream& s);

    void
    write(std::ostream& s) const;

    private:

    void 
    init(Args const& args);

    void 
    tableInit(InputGroup& table);

    std::vector<int> maxm_,
                     minm_,
                     niter_;
    std::vector<Real> cutoff_,
                      noise_;
    int nsweep_;
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
        nsweep_(int(v.size())-1),
        started_(false),
        v_(v),
        j_(1)
        { 
        last_val_ = v_[j_];
        }
    
    ~SweepSetter()
        {
        for(; j_ <= nsweep_; ++j_)
            {
            v_[j_] = last_val_;
            }
        }

    template <typename Arg>
    SweepSetter& 
    operator=(Arg val)
        {
        started_ = true;
        return operator,(val);
        }

    SweepSetter& 
    operator,(T val)
        {
        checkStarted();
        if(j_ > nsweep_) return *this;
        v_[j_] = val;
        ++j_;
        last_val_ = val;
        return *this;
        }

    SweepSetter& 
    operator,(const Args& args)
        { 
        checkStarted();
        if(args.defined("Repeat"))
            {
            if(j_ == 1) Error("No value to repeat");
            for(int n = 1; n < args.getInt("Repeat"); ++n, ++j_)
                {
                if(j_ > nsweep_) return *this;
                v_[j_] = last_val_;
                }
            }
        return *this;
        }

    using Func = std::function<T(int,int)>;

    //
    //Sets the remaining values of v_ by 
    //accepting (the address of a) function 
    //with signature:
    //
    //  T f(int sw, int nsweep) { ... }
    //
    //or an instance of a function object
    //which defines a method:
    //
    //  T operator()(int sw, int nsweep) const { ... } 
    //
    SweepSetter&
    operator,(Func f)
        {
        checkStarted();
        for(; j_ <= nsweep_ ; ++j_)
            {
            v_[j_] = f(j_,nsweep_);
            }
        last_val_ = v_.back();
        return *this;
        }

    private:

    const int nsweep_;
    bool started_;
    std::vector<T>& v_;
    int j_;
    T last_val_;

    void
    checkStarted() const
        {
        if(!started_) Error("SweepSetter notation is setter() = #,#,...;");
        }
    };

struct RampM
    {
    RampM(int start_m, int end_m,
          const Args& args = Global::args())
        :
        start_m_(start_m),
        end_m_(end_m),
        nwarm_(args.getInt("Warmup",-1))
        { }

    int
    operator()(int sw, int nsweep) const 
        { 
        const int actual_nwarm = (nwarm_ < 0 ? nsweep : std::min(nwarm_+1,nsweep));
        if(sw <= actual_nwarm)
            return (int) (start_m_ + (sw-1.)/(actual_nwarm-1.)*(end_m_-start_m_));
        else
            return end_m_;
        }

    private:
    const
    int start_m_,
        end_m_,
        nwarm_;
    };

struct ExpM
    {
    ExpM(int start_m, int end_m,
         const Args& args = Global::args())
        :
        start_m_(start_m),
        end_m_(end_m),
        exp_base_(args.getReal("ExpBase",2.))
        { }

    int
    operator()(int sw, int nsweep) const 
        { 
        int expm = int(start_m_*pow(exp_base_,sw-1));
        if(expm <= 0) return end_m_; //catch overflow
        return std::min(expm,end_m_);
        }

    private:
    const
    int start_m_,
        end_m_;
    const
    Real exp_base_;
    };

inline Sweeps::
Sweeps()
    :
    nsweep_(0)
    { }

inline Sweeps::
Sweeps(int nsw, 
       int min_m, 
       int max_m, 
       Real cut,
       Real noise)
  : nsweep_(nsw)
    {
    init({"Minm",min_m,
          "Maxm",max_m,
          "Cutoff",cut,
          "Noise",noise});
    }

inline Sweeps::
Sweeps(Args const& args)
    {
    nsweep_ = args.getInt("Nsweep");
    init(args);
    }

inline Sweeps::
Sweeps(int nsw, InputGroup& sweep_table)
    : 
    nsweep_(nsw)
    {
    tableInit(sweep_table);
    }

SweepSetter<int> inline Sweeps::
minm() { return SweepSetter<int>(minm_); }

SweepSetter<int> inline Sweeps::
maxm() { return SweepSetter<int>(maxm_); }

SweepSetter<Real> inline Sweeps::
cutoff() { return SweepSetter<Real>(cutoff_); }

SweepSetter<Real> inline Sweeps::
noise() { return SweepSetter<Real>(noise_); }

SweepSetter<int> inline Sweeps::
niter() { return SweepSetter<int>(niter_); }

void inline Sweeps::
nsweep(int val)
    { 
    if(val > nsweep_) 
        Error("Can't use nsweep accessor to increase number of sweeps.");
    nsweep_ = val; 
    }


void inline Sweeps::
init(Args const& args)
    {
    auto min_m = args.getInt("Minm",1);
    auto max_m = args.getInt("Maxm");
    auto cutoff = args.getReal("Cutoff");
    auto noise = args.getReal("Noise",0.);
    auto niter = args.getInt("Niter",2);

    minm_ = std::vector<int>(nsweep_+1,min_m);
    maxm_ = std::vector<int>(nsweep_+1,max_m);
    cutoff_ = std::vector<Real>(nsweep_+1,cutoff);
    niter_ = std::vector<int>(nsweep_+1,niter);
    noise_ = std::vector<Real>(nsweep_+1,noise);

    ////Set number of Davidson iterations
    //const int Max_niter = 9;
    //for(int s = 1; s <= std::min(4,nsweep_); ++s)
    //    {
    //    niter_.at(s) = std::max(Max_niter-s+1,2);
    //    }
    } //Sweeps::init

void inline Sweeps::
tableInit(InputGroup& table)
    {
    if(!table.GotoGroup()) 
        {
        Error("Couldn't find table " + table.name());
        }

    minm_ = std::vector<int>(nsweep_+1,0);
    maxm_ = std::vector<int>(nsweep_+1,0);
    cutoff_ = std::vector<Real>(nsweep_+1,0);
    niter_ = std::vector<int>(nsweep_+1,0);
    noise_ = std::vector<Real>(nsweep_+1,0);

    //printfln("Got nsweep_=%d",nsweep_);
    table.SkipLine(); //SkipLine so we can have a table key
    int n_last = nsweep_;
    for(int i = 1; i <= nsweep_; i++)
        {
        table.file() >> maxm_[i] >> minm_[i] >> cutoff_[i] >> niter_[i] >> noise_[i];
        //printfln("Line %d: %d %d %.2E %d %.2E",i,maxm_[i],minm_[i],cutoff_[i],niter_[i],noise_[i]);
        if(maxm_[i] == 0)
          {
          //printfln("Got maxm_[i]==0 at line number i=%d",i);
          n_last = i - 1 ;
          //printfln("Set n_last equal to n_last=%d",n_last);
          break ;
          }
        }
    //printfln("Filling n_last+1=%d to nsweep_=%d",n_last+1,nsweep_);
    for(int i = n_last + 1; i <= nsweep_; i++)
       {
       maxm_[i] = maxm_[n_last];
       minm_[i] = minm_[n_last];
       cutoff_[i] = cutoff_[n_last];
       niter_[i] = niter_[n_last];
       noise_[i] = noise_[n_last];
       }

    } //Sweeps::tableInit

void inline Sweeps::
write(std::ostream& s) const
    {
    itensor::write(s,maxm_);
    itensor::write(s,minm_);
    itensor::write(s,cutoff_);
    itensor::write(s,niter_);
    itensor::write(s,noise_);
    itensor::write(s,nsweep_);
    }

void inline Sweeps::
read(std::istream& s)
    {
    itensor::read(s,maxm_);
    itensor::read(s,minm_);
    itensor::read(s,cutoff_);
    itensor::read(s,niter_);
    itensor::read(s,noise_);
    itensor::read(s,nsweep_);
    }

inline std::ostream&
operator<<(std::ostream& s, const Sweeps& swps)
    {
    s << "Sweeps:\n";
    for(int sw = 1; sw <= swps.nsweep(); ++sw)
        {
        s << format("%d  Maxm=%d, Minm=%d, Cutoff=%.1E, Niter=%d, Noise=%.1E\n",
              sw,swps.maxm(sw),swps.minm(sw),swps.cutoff(sw),swps.niter(sw),swps.noise(sw));
        }
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

//one-site version of sweepnext
void inline
sweepnext1(int &b, int &ha, int N, int min_b = 1)
    {
    const int inc = (ha==1 ? +1 : -1);
    b += inc;
    if(b == (ha==1 ? N+1 : min_b-1))
        {
        b -= inc;
        ++ha;
        }
    }

} //namespace itensor


#endif //__ITENSOR_SWEEPS_HEADER_H
