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
#ifndef __ITENSOR_DMRG_H
#define __ITENSOR_DMRG_H

#include "itensor/iterativesolvers.h"
#include "itensor/mps/localmposet.h"
#include "itensor/mps/localmpo_mps.h"
#include "itensor/mps/sweeps.h"
#include "itensor/mps/DMRGObserver.h"
#include "itensor/util/cputime.h"


namespace itensor {

template<class LocalOpT>
Real
DMRGWorker(MPS & psi,
           LocalOpT & PH,
           Sweeps const& sweeps,
           Args const& args = Args::global());

template<class LocalOpT>
Real
DMRGWorker(MPS & psi,
           LocalOpT & PH,
           Sweeps const& sweeps,
           DMRGObserver & obs,
           Args args = Args::global());

//
// Available DMRG methods:
//

//
//DMRG with an MPO
//
Real inline
dmrg(MPS & psi, 
     MPO const& H, 
     Sweeps const& sweeps,
     Args const& args = Args::global())
    {
    LocalMPO PH(H,args);
    Real energy = DMRGWorker(psi,PH,sweeps,args);
    return energy;
    }

//
//DMRG with an MPO
//Version that takes a starting guess MPS
//and returns the optimized MPS
//
std::tuple<Real,MPS> inline
dmrg(MPO const& H,
     MPS const& psi0,
     Sweeps const& sweeps,
     Args const& args = Args::global())
    {
    auto psi = psi0;
    auto energy = dmrg(psi,H,sweeps,args);
    return std::tuple<Real,MPS>(energy,psi);
    }

//
//DMRG with an MPO and custom DMRGObserver
//
Real inline
dmrg(MPS& psi, 
     MPO const& H, 
     Sweeps const& sweeps, 
     DMRGObserver & obs,
     Args const& args = Args::global())
    {
    LocalMPO PH(H,args);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,args);
    return energy;
    }

//
//DMRG with an MPO and custom DMRGObserver
//Version that takes a starting guess MPS
//and returns the optimized MPS
//
std::tuple<Real,MPS> inline
dmrg(MPO const& H,
     MPS const& psi0,
     Sweeps const& sweeps,
     DMRGObserver & obs,
     Args const& args = Args::global())
    {
    auto psi = psi0;
    auto energy = dmrg(psi,H,sweeps,obs,args);
    return std::tuple<Real,MPS>(energy,psi);
    }

//
//DMRG with an MPO and boundary tensors LH, RH
// LH - H1 - H2 - ... - HN - RH
//(ok if one or both of LH, RH default constructed)
//
Real inline
dmrg(MPS& psi, 
     MPO const& H, 
     ITensor const& LH, 
     ITensor const& RH,
     Sweeps const& sweeps,
     Args const& args = Args::global())
    {
    LocalMPO PH(H,LH,RH,args);
    Real energy = DMRGWorker(psi,PH,sweeps,args);
    return energy;
    }

//
//DMRG with an MPO and boundary tensors LH, RH
// LH - H1 - H2 - ... - HN - RH
//(ok if one or both of LH, RH default constructed)
//Version that takes a starting guess MPS
//and returns the optimized MPS
//
std::tuple<Real,MPS> inline
dmrg(MPO const& H,
     ITensor const& LH,
     ITensor const& RH,
     MPS const& psi0,
     Sweeps const& sweeps,
     Args const& args = Args::global())
    {
    auto psi = psi0;
    auto energy = dmrg(psi,H,LH,RH,sweeps,args);
    return std::tuple<Real,MPS>(energy,psi);
    }

//
//DMRG with an MPO and boundary tensors LH, RH
//and a custom observer
//
Real inline
dmrg(MPS& psi, 
     MPO const& H, 
     ITensor const& LH, 
     ITensor const& RH,
     Sweeps const& sweeps, 
     DMRGObserver& obs,
     Args const& args = Args::global())
    {
    LocalMPO PH(H,LH,RH,args);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,args);
    return energy;
    }

std::tuple<Real,MPS> inline
dmrg(MPO const& H,
     ITensor const& LH,
     ITensor const& RH,
     MPS const& psi0,
     Sweeps const& sweeps,
     DMRGObserver& obs,
     Args const& args = Args::global())
    {
    auto psi = psi0;
    auto energy = dmrg(psi,H,LH,RH,sweeps,obs,args);
    return std::tuple<Real,MPS>(energy,psi);
    }

//
//DMRG with a set of MPOs (lazily summed)
//(H vector is 0-indexed)
//
Real inline
dmrg(MPS& psi, 
     std::vector<MPO> const& Hset, 
     Sweeps const& sweeps,
     Args const& args = Args::global())
    {
    LocalMPOSet PH(Hset,args);
    Real energy = DMRGWorker(psi,PH,sweeps,args);
    return energy;
    }

std::tuple<Real,MPS> inline
dmrg(std::vector<MPO> const& Hset,
     MPS const& psi0,
     Sweeps const& sweeps,
     Args const& args = Args::global())
    {
    auto psi = psi0;
    auto energy = dmrg(psi,Hset,sweeps,args);
    return std::tuple<Real,MPS>(energy,psi);
    }

//
//DMRG with a set of MPOs and a custom DMRGObserver
//(H vector is 0-indexed)
//
Real inline
dmrg(MPS& psi, 
     std::vector<MPO> const& Hset, 
     Sweeps const& sweeps, 
     DMRGObserver& obs,
     Args const& args = Args::global())
    {
    LocalMPOSet PH(Hset,args);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,args);
    return energy;
    }

std::tuple<Real,MPS> inline
dmrg(std::vector<MPO> const& Hset,
     MPS const& psi0,
     Sweeps const& sweeps,
     DMRGObserver& obs,
     Args const& args = Args::global())
    {
    auto psi = psi0;
    auto energy = dmrg(psi,Hset,sweeps,obs,args);
    return std::tuple<Real,MPS>(energy,psi);
    }

//
//DMRG with a single Hamiltonian MPO and a set of 
//MPS to orthogonalize against
//(psis vector is 0-indexed)
//Named Args recognized:
// Weight - real number w > 0; calling dmrg(psi,H,psis,sweeps,Args("Weight",w))
//          sets the effective Hamiltonian to be
//          H + w * (|0><0| + |1><1| + ...) where |0> = psis[0], |1> = psis[1]
//          etc.
//
Real inline
dmrg(MPS& psi, 
     MPO const& H, 
     std::vector<MPS> const& psis, 
     Sweeps const& sweeps, 
     Args const& args = Args::global())
    {
    if(hasQNs(psi))
        {
        auto psi_qn = totalQN(psi);
        for(auto n : range(psis))
            {
            auto qn_n = totalQN(psis[n]);
            if(qn_n != psi_qn)
                {
                printfln("totalQN of initial state:  %s",psi_qn);
                printfln("totalQN of state n=%d (n is 0-indexed): %s",n,qn_n);
                Error("Excited-state DMRG intended for states with same totalQN");
                }
            }
        }
    LocalMPO_MPS PH(H,psis,args);
    Real energy = DMRGWorker(psi,PH,sweeps,args);
    return energy;
    }

std::tuple<Real,MPS> inline
dmrg(MPO const& H,
     std::vector<MPS> const& psis,
     MPS const& psi0,
     Sweeps const& sweeps,
     Args const& args = Args::global())
    {
    auto psi = psi0;
    auto energy = dmrg(psi,H,psis,sweeps,args);
    return std::tuple<Real,MPS>(energy,psi);
    }

//
//DMRG with a single Hamiltonian MPO, 
//a set of MPS to orthogonalize against, 
//and a custom DMRGObserver.
//(psis vector is 0-indexed)
//Named Args recognized:
// Weight - real number w > 0; calling dmrg(psi,H,psis,sweeps,Args("Weight",w))
//          sets the effective Hamiltonian to be
//          H + w * (|0><0| + |1><1| + ...) where |0> = psis[0], |1> = psis[1]
//          etc.
//
Real inline
dmrg(MPS & psi, 
     MPO const& H, 
     std::vector<MPS> const& psis, 
     Sweeps const& sweeps, 
     DMRGObserver& obs, 
     Args const& args = Args::global())
    {
    LocalMPO_MPS PH(H,psis,args);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,args);
    return energy;
    }

std::tuple<Real,MPS> inline
dmrg(MPO const& H,
     std::vector<MPS> const& psis,
     MPS const& psi0,
     Sweeps const& sweeps,
     DMRGObserver& obs, 
     Args const& args = Args::global())
    {
    auto psi = psi0;
    auto energy = dmrg(psi,H,psis,sweeps,obs,args);
    return std::tuple<Real,MPS>(energy,psi);
    }


//
// DMRGWorker
//

template<class LocalOpT>
Real
DMRGWorker(MPS & psi,
           LocalOpT & PH,
           Sweeps const& sweeps,
           Args const& args)
    {
    DMRGObserver obs(psi,args);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,args);
    return energy;
    }

template<class LocalOpT>
Real
DMRGWorker(MPS & psi,
           LocalOpT & PH,
           Sweeps const& sweeps,
           DMRGObserver & obs,
           Args args)
    {
    if( args.defined("WriteM") )
      {
      if( args.defined("WriteDim") )
        {
        Global::warnDeprecated("Args WirteM and WriteDim are both defined. WriteM is deprecated in favor of WriteDim, WriteDim will be used.");
        }
      else
        {
        Global::warnDeprecated("Arg WriteM is deprecated in favor of WriteDim.");
        args.add("WriteDim",args.getInt("WriteM"));
        }
      }

    // Truncate blocks of degenerate singular values (or not)
    args.add("RespectDegenerate",args.getBool("RespectDegenerate",true));

    const bool silent = args.getBool("Silent",false);
    if(silent)
        {
        args.add("Quiet",true);
        args.add("PrintEigs",false);
        args.add("NoMeasure",true);
        args.add("DebugLevel",0);
        }
    const bool quiet = args.getBool("Quiet",false);
    const int debug_level = args.getInt("DebugLevel",(quiet ? 0 : 1));

    const int N = length(psi);
    Real energy = NAN;

    psi.position(1);

    args.add("DebugLevel",debug_level);
    args.add("DoNormalize",true);
    
    for(int sw = 1; sw <= sweeps.nsweep(); ++sw)
        {
        cpu_time sw_time;
        args.add("Sweep",sw);
        args.add("NSweep",sweeps.nsweep());
        args.add("Cutoff",sweeps.cutoff(sw));
        args.add("MinDim",sweeps.mindim(sw));
        args.add("MaxDim",sweeps.maxdim(sw));
        args.add("Noise",sweeps.noise(sw));
        args.add("MaxIter",sweeps.niter(sw));

        if(!PH.doWrite()
           && args.defined("WriteDim")
           && sweeps.maxdim(sw) >= args.getInt("WriteDim"))
            {
            if(!quiet)
                {
                println("\nTurning on write to disk, write_dir = ",
                        args.getString("WriteDir","./"));
                }

            //psi.doWrite(true);
            PH.doWrite(true,args);
            }

        for(int b = 1, ha = 1; ha <= 2; sweepnext(b,ha,N))
            {
            if(!quiet)
                {
                printfln("Sweep=%d, HS=%d, Bond=%d/%d",sw,ha,b,(N-1));
                }

TIMER_START(1);
            PH.position(b,psi);
TIMER_STOP(1);

TIMER_START(2);
            auto phi = psi(b)*psi(b+1);
TIMER_STOP(2);

TIMER_START(3);
            energy = davidson(PH,phi,args);
TIMER_STOP(3);
            
TIMER_START(4);
            auto spec = psi.svdBond(b,phi,(ha==1?Fromleft:Fromright),PH,args);
TIMER_STOP(4);

            if(!quiet)
                { 
                printfln("    Truncated to Cutoff=%.1E, Min_dim=%d, Max_dim=%d",
                          sweeps.cutoff(sw),
                          sweeps.mindim(sw), 
                          sweeps.maxdim(sw) );
                printfln("    Trunc. err=%.1E, States kept: %s",
                         spec.truncerr(),
                         showDim(linkIndex(psi,b)) );
                }

            obs.lastSpectrum(spec);

            args.add("AtBond",b);
            args.add("HalfSweep",ha);
            args.add("Energy",energy); 
            args.add("Truncerr",spec.truncerr()); 

            obs.measure(args);

            } //for loop over b

        if(!silent)
            {
            auto sm = sw_time.sincemark();
            printfln("    Sweep %d/%d CPU time = %s (Wall time = %s)",
                      sw,sweeps.nsweep(),showtime(sm.time),showtime(sm.wall));
#ifdef COLLECT_TIMES
            println(timers());
            timers().reset();
#endif
            }

        if(obs.checkDone(args)) break;
    
        } //for loop over sw
    
    psi.normalize();

    return energy;
    }

} //namespace itensor


#endif
