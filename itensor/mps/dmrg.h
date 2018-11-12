//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
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

//
// Available DMRG methods:
//

//
//DMRG with an MPO
//
template <class Tensor>
Real
dmrg(MPSt<Tensor>& psi, 
     const MPOt<Tensor>& H, 
     const Sweeps& sweeps,
     const Args& args = Global::args())
    {
    LocalMPO<Tensor> PH(H,args);
    Real energy = DMRGWorker(psi,PH,sweeps,args);
    return energy;
    }

//
//DMRG with an MPO and custom DMRGObserver
//
template <class Tensor>
Real
dmrg(MPSt<Tensor>& psi, 
     const MPOt<Tensor>& H, 
     const Sweeps& sweeps, 
     DMRGObserver<Tensor>& obs,
     const Args& args = Global::args())
    {
    LocalMPO<Tensor> PH(H,args);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,args);
    return energy;
    }

//
//DMRG with an MPO and boundary tensors LH, RH
// LH - H1 - H2 - ... - HN - RH
//(ok if one or both of LH, RH default constructed)
//
template <class Tensor>
Real
dmrg(MPSt<Tensor>& psi, 
     const MPOt<Tensor>& H, 
     const Tensor& LH, const Tensor& RH,
     const Sweeps& sweeps,
     const Args& args = Global::args())
    {
    LocalMPO<Tensor> PH(H,LH,RH,args);
    Real energy = DMRGWorker(psi,PH,sweeps,args);
    return energy;
    }

//
//DMRG with an MPO and boundary tensors LH, RH
//and a custom observer
//
template <class Tensor>
Real
dmrg(MPSt<Tensor>& psi, 
     const MPOt<Tensor>& H, 
     const Tensor& LH, const Tensor& RH,
     const Sweeps& sweeps, 
     DMRGObserver<Tensor>& obs,
     const Args& args = Global::args())
    {
    LocalMPO<Tensor> PH(H,LH,RH,args);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,args);
    return energy;
    }

//
//DMRG with a set of MPOs (lazily summed)
//(H vector is 0-indexed)
//
template <class Tensor>
Real
dmrg(MPSt<Tensor>& psi, 
     const std::vector<MPOt<Tensor> >& Hset, 
     const Sweeps& sweeps,
     const Args& args = Global::args())
    {
    LocalMPOSet<Tensor> PH(Hset,args);
    Real energy = DMRGWorker(psi,PH,sweeps,args);
    return energy;
    }

//
//DMRG with a set of MPOs and a custom DMRGObserver
//(H vector is 0-indexed)
//
template <class Tensor>
Real 
dmrg(MPSt<Tensor>& psi, 
     const std::vector<MPOt<Tensor> >& Hset, 
     const Sweeps& sweeps, 
     DMRGObserver<Tensor>& obs,
     const Args& args = Global::args())
    {
    LocalMPOSet<Tensor> PH(Hset,args);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,args);
    return energy;
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
template <class Tensor>
Real
dmrg(MPSt<Tensor>& psi, 
     const MPOt<Tensor>& H, 
     const std::vector<MPSt<Tensor> >& psis, 
     const Sweeps& sweeps, 
     const Args& args = Global::args())
    {
    LocalMPO_MPS<Tensor> PH(H,psis,args);
    Real energy = DMRGWorker(psi,PH,sweeps,args);
    return energy;
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
template <class Tensor>
Real
dmrg(MPSt<Tensor>& psi, 
     const MPOt<Tensor>& H, 
     const std::vector<MPSt<Tensor> >& psis, 
     const Sweeps& sweeps, 
     DMRGObserver<Tensor>& obs, 
     const Args& args = Global::args())
    {
    LocalMPO_MPS<Tensor> PH(H,psis,args);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,args);
    return energy;
    }



//
// DMRGWorker
//

template <class Tensor, class LocalOpT>
Real inline
DMRGWorker(MPSt<Tensor>& psi,
           LocalOpT& PH,
           const Sweeps& sweeps,
           const Args& args = Global::args())
    {
    DMRGObserver<Tensor> obs(psi,args);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,args);
    return energy;
    }

template <class Tensor, class LocalOpT>
Real
DMRGWorker(MPSt<Tensor>& psi,
           LocalOpT& PH,
           const Sweeps& sweeps,
           DMRGObserver<Tensor>& obs,
           Args args = Global::args())
    {
    const bool quiet = args.getBool("Quiet",false);
    const int debug_level = args.getInt("DebugLevel",(quiet ? 0 : 1));
    const int numCenter = args.getInt("NumCenter",2);

    const int N = psi.N();
    Real energy = NAN;

    psi.position(1);

    args.add("DebugLevel",debug_level);
    args.add("DoNormalize",true);


    if(numCenter == 2)
		{	
        for(int sw = 1; sw <= sweeps.nsweep(); ++sw)
            {
            cpu_time sw_time;
            args.add("Sweep",sw);
            args.add("NSweep",sweeps.nsweep());
            args.add("Cutoff",sweeps.cutoff(sw));
            args.add("Minm",sweeps.minm(sw));
            args.add("Maxm",sweeps.maxm(sw));
            args.add("Noise",sweeps.noise(sw));
            args.add("MaxIter",sweeps.niter(sw));
    
            if(!PH.doWrite()
               && args.defined("WriteM")
               && sweeps.maxm(sw) >= args.getInt("WriteM"))
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
    
                PH.position(b,psi);
    
                auto phi = psi.A(b)*psi.A(b+1);
    
                energy = davidson(PH,phi,args);
                
                auto spec = psi.svdBond(b,phi,(ha==1?Fromleft:Fromright),PH,args);
    
    
                if(!quiet)
                    { 
                    printfln("    Truncated to Cutoff=%.1E, Min_m=%d, Max_m=%d",
                              sweeps.cutoff(sw),
                              sweeps.minm(sw), 
                              sweeps.maxm(sw) );
                    printfln("    Trunc. err=%.1E, States kept: %s",
                             spec.truncerr(),
                             showm(linkInd(psi,b)) );
                    }
    
                obs.lastSpectrum(spec);
    
                args.add("AtBond",b);
                args.add("HalfSweep",ha);
                args.add("Energy",energy); 
                args.add("Truncerr",spec.truncerr()); 
    
                obs.measure(args);
    
                } //for loop over b
    
            auto sm = sw_time.sincemark();
            printfln("    Sweep %d/%d CPU time = %s (Wall time = %s)",
                      sw,sweeps.nsweep(),showtime(sm.time),showtime(sm.wall));
    
            if(obs.checkDone(args)) break;
        
            } //for loop over sw
        }
    else
        {	
        if(numCenter == 3)
            {
            for(int sw = 1; sw <= sweeps.nsweep(); ++sw)
            	{
                cpu_time sw_time;
                args.add("Sweep",sw);
                args.add("Cutoff",sweeps.cutoff(sw));
                args.add("Minm",sweeps.minm(sw));
                args.add("Maxm",sweeps.maxm(sw));
                args.add("Noise",sweeps.noise(sw));
                args.add("MaxIter",sweeps.niter(sw));
                    
                if(!PH.doWrite()
                    && args.defined("WriteM")
                    && sweeps.maxm(sw) >= args.getInt("WriteM"))
                	{
                    if(!quiet)
                    	{
                        println("\nTurning on write to disk, write_dir = ",
                                    args.getString("WriteDir","./"));
                    	}
                        
                    //psi.doWrite(true);
                    PH.doWrite(true);//default do_write_ is false, "true" argument to triger the initWrite, do_write_ = true after done. This will only be called at sw==1. 
                	}
                
                Tensor temp;
                for(int b = 1, ha = 1; ha <= 2; sweepnext3(b,ha,N))
                	{
                    if(!quiet)
                    	{
                        printfln("Sweep=%d, HS=%d, Bond=(%d,%d,%d)",sw,ha,b,(b+1),(b+2));
                    	}
                        
                    PH.position(b,psi);//make the MPO represented in the respective local states
                    
                    Tensor phi;
                    if(ha == 1)
        				{
                        if(b == 1)  
        					phi = psi.A(b)*psi.A(b+1)*psi.A(b+2);
                        else        
        					phi = temp*psi.A(b+2);
                    	}
                    else
        				{
                        if(b == N-2)    
        					phi = psi.A(b)*psi.A(b+1)*psi.A(b+2);
                        else            
        					phi = psi.A(b)*temp;
                    	}
        				
                    energy = davidson(PH,phi,args);
                    //since the svd for a three-site tensor is only used here,
                    //so we do not define a new function to do it.
        	    //now define a new svdBond3 and setBond3, setBond3=setBond.
                    //Now only deal with the case that "noise=0 && cutoff<1E-12", i.e. using svd but not density matrix, since denmatdecomp need commonindex of A[b] and A[b+1] and modify the combiner function.
                    //and forget about the donormalize args. --> now deal with the donormalize args.
        	    //Through setbond to move OC in MPS (not in MPO, which was done before by PH.position) and write and read from disk
        	    ////////////however,l_orth_lim_ etc. setting by Aref is behind the actual position during the sweep,though it's correct in the end. Will this cause problem?
        			auto spec = psi.svdBond3(b,phi,(ha==1?Fromleft:Fromright),PH,temp,args);
        
                    if(!quiet)
                    	{
                        printfln("    Truncated to Cutoff=%.1E, Min_m=%d, Max_m=%d",
                                     sweeps.cutoff(sw),
                                     sweeps.minm(sw),
                                     sweeps.maxm(sw) );
                        printfln("    Trunc. err=%.1E, States kept: %s",
                                   spec.truncerr() );
                    	}
                        
                    obs.lastSpectrum(spec);
                        
                    args.add("AtBond",b);
                    args.add("HalfSweep",ha);
                    args.add("Energy",energy);
                        
                    obs.measure(args);
        			            
                } //for loop over b
                
                    
                auto sm = sw_time.sincemark();
                printfln("    Sweep %d CPU time = %s (Wall time = %s)",
                            sw,showtime(sm.time),showtime(sm.wall));
        
                
                if(obs.checkDone(args)) break;
                    
            	} //for loop over sw        
            }
        else
            Error("Only support 2 and 3 sites algorithm presently.");
        }
    
    psi.normalize();

    return energy;
    }

} //namespace itensor


#endif
