//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
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
//      maxdim   mindim  cutoff  niter  noise
//      20       20      1E-8    4      1E-8
//      40       20      1E-8    3      1E-9
//      80       20      1E-10   2      1E-10
//      160      20      1E-12   2      0
//      240      20      1E-12   2      0
//      }
//
class Sweeps
    {
    public:

    //Constructors --------------

    Sweeps();

    Sweeps(int nsweeps, 
           int mindim = 1, 
           int maxdim = 500, 
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
    mindim(int sw) const { return mindim_.at(sw); }
    void 
    setmindim(int sw, int val) { mindim_.at(sw) = val; }

    //Use as sweeps.mindim() = 20,20,10; (all remaining set to 10)
    SweepSetter<int> 
    mindim();

    int 
    maxdim(int sw) const { return maxdim_.at(sw); }
    void 
    setmaxdim(int sw, int val) { maxdim_.at(sw) = val; }

    //Use as sweeps.maxdim() = 50,50,100,100,500; (all remaining set to 500)
    SweepSetter<int> 
    maxdim();

    // Deprecated, use maxdim() instead
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
    init(Args args);

    void 
    tableInit(InputGroup& table);

    std::vector<int> maxdim_,
                     mindim_,
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
          const Args& args = Args::global())
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
         Args const& args = Args::global())
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
       int min_dim, 
       int max_dim, 
       Real cut,
       Real noise)
  : nsweep_(nsw)
    {
    init({"MinDim",min_dim,
          "MaxDim",max_dim,
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
mindim() { return SweepSetter<int>(mindim_); }

SweepSetter<int> inline Sweeps::
maxdim() { return SweepSetter<int>(maxdim_); }

SweepSetter<int> inline Sweeps::
maxm()
  {
  Global::warnDeprecated("maxm() is deprecated, use maxdim() instead.");
  return SweepSetter<int>(maxdim_);
  }

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
init(Args args)
    {
    if( args.defined("Minm") )
      {
      if( args.defined("MinDim") )
        {
        Global::warnDeprecated("Args Minm and MinDim are both defined. Minm is deprecated in favor of MinDim, MinDim will be used.");
        }
      else
        {
        Global::warnDeprecated("Arg Minm is deprecated in favor of MinDim.");
        args.add("MinDim",args.getInt("Minm"));
        }
      }

    if( args.defined("Maxm") )
      {
      if( args.defined("MaxDim") )
        {
        Global::warnDeprecated("Args Maxm and MaxDim are both defined. Maxm is deprecated in favor of MaxDim, MaxDim will be used.");
        }
      else
        {
        Global::warnDeprecated("Arg Maxm is deprecated in favor of MaxDim.");
        args.add("MaxDim",args.getInt("Maxm"));
        }
      }

    auto min_dim = args.getInt("MinDim",1);
    auto max_dim = args.getInt("MaxDim");
    auto cutoff = args.getReal("Cutoff");
    auto noise = args.getReal("Noise",0.);
    auto niter = args.getInt("Niter",2);

    mindim_ = std::vector<int>(nsweep_+1,min_dim);
    maxdim_ = std::vector<int>(nsweep_+1,max_dim);
    cutoff_ = std::vector<Real>(nsweep_+1,cutoff);
    niter_ = std::vector<int>(nsweep_+1,niter);
    noise_ = std::vector<Real>(nsweep_+1,noise);
    } //Sweeps::init

void inline Sweeps::
tableInit(InputGroup& table)
    {
    if(!table.GotoGroup()) 
        {
        Error("Couldn't find table " + table.name());
        }

    mindim_ = std::vector<int>(nsweep_+1,0);
    maxdim_ = std::vector<int>(nsweep_+1,0);
    cutoff_ = std::vector<Real>(nsweep_+1,0);
    niter_ = std::vector<int>(nsweep_+1,0);
    noise_ = std::vector<Real>(nsweep_+1,0);

    //printfln("Got nsweep_=%d",nsweep_);
    table.SkipLine(); //SkipLine so we can have a table key
    int n_last = nsweep_;
    for(int i = 1; i <= nsweep_; i++)
        {
        table.file() >> maxdim_[i] >> mindim_[i] >> cutoff_[i] >> niter_[i] >> noise_[i];
        //printfln("Line %d: %d %d %.2E %d %.2E",i,maxdim_[i],mindim_[i],cutoff_[i],niter_[i],noise_[i]);
        if(maxdim_[i] == 0)
          {
          //printfln("Got maxdim_[i]==0 at line number i=%d",i);
          n_last = i - 1 ;
          //printfln("Set n_last equal to n_last=%d",n_last);
          break ;
          }
        }
    //printfln("Filling n_last+1=%d to nsweep_=%d",n_last+1,nsweep_);
    for(int i = n_last + 1; i <= nsweep_; i++)
       {
       maxdim_[i] = maxdim_[n_last];
       mindim_[i] = mindim_[n_last];
       cutoff_[i] = cutoff_[n_last];
       niter_[i] = niter_[n_last];
       noise_[i] = noise_[n_last];
       }

    } //Sweeps::tableInit

void inline Sweeps::
write(std::ostream& s) const
    {
    itensor::write(s,maxdim_);
    itensor::write(s,mindim_);
    itensor::write(s,cutoff_);
    itensor::write(s,niter_);
    itensor::write(s,noise_);
    itensor::write(s,nsweep_);
    }

void inline Sweeps::
read(std::istream& s)
    {
    itensor::read(s,maxdim_);
    itensor::read(s,mindim_);
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
        s << format("%d  MaxDim=%d, MinDim=%d, Cutoff=%.1E, Niter=%d, Noise=%.1E\n",
              sw,swps.maxdim(sw),swps.mindim(sw),swps.cutoff(sw),swps.niter(sw),swps.noise(sw));
        }
    return s;
    }

void inline
sweepnext(int &b, int &ha, int N, Args const& args = Args::global())
    {
    const int numCenter = args.getInt("NumCenter",2);
    const int min_b = args.getInt("Minb",1);
    const int inc = (ha==1 ? +1 : -1);
    b += inc;
    if(b == (ha==1 ? N+2-numCenter : min_b-1))
        {
        b -= inc;
        ++ha;
        }
    }

void inline
sweepnext(int &b, int &ha, int N, int min_b)
    {
    auto args = Args("NumCenter=",2,"Minb=",min_b);
    sweepnext(b,ha,N,args);
    }

void inline
sweepnext1(int &b, int &ha, int N, int min_b = 1)
    {
	auto args = Args("NumCenter=",1,"Minb=",min_b);
    sweepnext(b,ha,N,args);
    }

} //namespace itensor


#endif //__ITENSOR_SWEEPS_HEADER_H
