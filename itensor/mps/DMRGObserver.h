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
#ifndef __ITENSOR_DMRGOBSERVER_H
#define __ITENSOR_DMRGOBSERVER_H
#include "itensor/mps/mps.h"
#include "itensor/mps/observer.h"
#include "itensor/spectrum.h"

namespace itensor {

//
// Class for monitoring DMRG calculations.
// The measure and checkDone methods are virtual
// so that behavior can be customized in a
// derived class.
//

class DMRGObserver : public Observer
    {
    public:
    
    DMRGObserver(MPS const& psi, 
                 Args const& args = Args::global());

    virtual ~DMRGObserver() { }

    void virtual
    measure(Args const& args = Args::global());
    
    bool virtual
    checkDone(Args const& args = Args::global());

    void virtual
    lastSpectrum(Spectrum const& spec) { last_spec_ = spec; }

    MPS const& 
    psi() const { return psi_; }
    
    Spectrum const&
    spectrum() const { return last_spec_; }

    private:

    /////////////

    MPS const& psi_;
    Real energy_errgoal; //Stop DMRG once energy has converged to this precision
    bool printeigs;      //Print slowest decaying eigenvalues after every sweep
    int max_eigs;
    Real max_te;
    bool done_;
    Real last_energy_;
    Spectrum last_spec_;

    /////////////

    }; // class DMRGObserver

inline DMRGObserver::
DMRGObserver(MPS const& psi, Args const& args) 
    : 
    psi_(psi),
    energy_errgoal(args.getReal("EnergyErrgoal",-1)), 
    printeigs(args.getBool("PrintEigs",true)),
    max_eigs(-1),
    max_te(-1),
    done_(false),
    last_energy_(1000)
    //default_ops_(psi.sites().defaultOps())
    { 
    }

void inline DMRGObserver::
measure(Args const& args)
    {
    auto N = length(psi_);
    auto b = args.getInt("AtBond",1);
    auto sw = args.getInt("Sweep",0);
    auto nsweep = args.getInt("NSweep",0);
    auto ha = args.getInt("HalfSweep",0);
    auto energy = args.getReal("Energy",0);
    auto silent = args.getBool("Silent",false);

    //if(!args.getBool("Quiet",false) && !args.getBool("NoMeasure",false))
    //    {
    //    if(b < N && b > 0)
    //        {
    //        auto wfb = psi_(b)*psi_(b+1);
    //        //for(const std::string& opname : default_ops_)
    //        //    {
    //        //    auto sb = IndexT(psi_.sites()(b));
    //        //    auto z = (dag(prime(wfb,sb))*psi_.sites().op(opname,b)*wfb).eltC();
    //        //    //auto z = (prime(wfb,psi_.sites()(b))*psi_.sites().op(opname,b)*wfb).eltC();
    //        //    if(std::fabs(z.imag()) < 1E-14)
    //        //        printfln("<%s>(%d) = %.10E",opname,b,z.real());
    //        //    else
    //        //        printfln("<%s>(%d) = (%.10E,%.10E)",opname,b,z.real(),z.imag());
    //        //    }
    //        }
    //    }

    if(!silent && printeigs)
        {
        if(b == N/2 && ha == 2)
            {
            println();
            auto center_eigs = last_spec_.eigsKept();
            // Normalize eigs
            Real norm_eigs = 0;
            for(auto& p : center_eigs)
              norm_eigs += p;
            center_eigs /= norm_eigs;
            // Calculate entropy
            Real S = 0;
            for(auto& p : center_eigs)
                {
                if(p > 1E-13) S += p*log(p);
                }
            S *= -1;
            printfln("    vN Entropy at center bond b=%d = %.12f",N/2,S);
            printf(  "    Eigs at center bond b=%d: ",N/2);
            auto ten = decltype(center_eigs.size())(10);
            for(auto j : range(std::min(center_eigs.size(),ten)))
                {
                auto eig = center_eigs(j);
                if(eig < 1E-3) break;
                printf("%.4f ",eig);
                }
            println();
            }
        }

    max_eigs = std::max(max_eigs,last_spec_.numEigsKept());
    max_te = std::max(max_te,last_spec_.truncerr());
    if(!silent)
        {
        if(b == 1 && ha == 2) 
            {
            if(!printeigs) println();
            auto swstr = (nsweep>0) ? format("%d/%d",sw,nsweep) 
                                    : format("%d",sw);
            println("    Largest link dim during sweep ",swstr," was ",(max_eigs > 1 ? max_eigs : 1));
            max_eigs = -1;
            println("    Largest truncation error: ",(max_te > 0 ? max_te : 0.));
            max_te = -1;
            printfln("    Energy after sweep %s is %.12f",swstr,energy);
            }
        }

    }


bool inline DMRGObserver::
checkDone(Args const& args)
    {
    const int sw = args.getInt("Sweep",0);
    const Real energy = args.getReal("Energy",0);
    
    if(sw == 1) last_energy_ = 1000;
    if(energy_errgoal > 0 && sw%2 == 0)
        {
        Real dE = std::fabs(energy-last_energy_);
        if(dE < energy_errgoal)
            {
            printfln("    Energy error goal met (dE = %.3E < %.3E); returning after %d sweeps.",
                      dE, energy_errgoal, sw);
            last_energy_ = 1000;
            return true;
            }
        }
    last_energy_ = energy;

    //If STOP_DMRG found, will return true (i.e. done) once, but 
    //outer calling using same Observer may continue running e.g. infinite dmrg calling finite dmrg.
    if(fileExists("STOP_DMRG"))
        {
        println("File STOP_DMRG found: stopping this DMRG run after sweep ",sw);
        std::remove("STOP_DMRG");
        return true;
        }

    //Set done_ flag to true so any outer callers using this Observer will also terminate.
    if(fileExists("STOP_DMRG_ALL"))
        {
        println("File STOP_DMRG_ALL found: stopping this run after sweep ",sw);
        std::remove("STOP_DMRG_ALL");
        done_ = true;
        return done_;
        }
    
    return done_;
    }

} //namespace itensor

#endif // __ITENSOR_DMRGOBSERVER_H
