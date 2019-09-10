#ifndef __ITENSOR_TDVP_H
#define __ITENSOR_TDVP_H

#include "itensor/iterativesolvers.h"
#include "itensor/mps/localmposet.h"
#include "itensor/mps/sweeps.h"
#include "itensor/mps/DMRGObserver.h"
#include "itensor/util/cputime.h"


namespace itensor {

template <class LocalOpT>
Real
TDVPWorker(MPS & psi,
           LocalOpT& PH,
		   Real t,
           const Sweeps& sweeps,
           const Args& args = Args::global());

template <class LocalOpT>
Real
TDVPWorker(MPS & psi,
           LocalOpT& PH,
		   Real t,
           const Sweeps& sweeps,
           DMRGObserver & obs,
           Args args = Args::global());

//
// Available TDVP methods:
// second order integrator: sweep left-to-right and right-to-left
//

//
//TDVP with an MPO
//
Real inline
tdvp(MPS & psi, 
     MPO const& H,
	 Real t, 
     const Sweeps& sweeps,
     const Args& args = Args::global())
    {
    LocalMPO PH(H,args);
    Real energy = TDVPWorker(psi,PH,t,sweeps,args);
    return energy;
    }

//
//TDVP with an MPO and custom DMRGObserver
//
Real inline
tdvp(MPS & psi, 
     MPO const& H, 
	 Real t,
     const Sweeps& sweeps, 
     DMRGObserver & obs,
     const Args& args = Args::global())
    {
    LocalMPO PH(H,args);
    Real energy = TDVPWorker(psi,PH,t,sweeps,obs,args);
    return energy;
    }

//
//TDVP with an MPO and boundary tensors LH, RH
// LH - H1 - H2 - ... - HN - RH
//(ok if one or both of LH, RH default constructed)
//
Real inline
tdvp(MPS & psi, 
     MPO const& H, 
	 Real t,
     ITensor const& LH, 
	 ITensor const& RH,
     const Sweeps& sweeps,
     const Args& args = Args::global())
    {
    LocalMPO PH(H,LH,RH,args);
    Real energy = TDVPWorker(psi,PH,t,sweeps,args);
    return energy;
    }

//
//TDVP with an MPO and boundary tensors LH, RH
//and a custom observer
//
Real inline
tdvp(MPS & psi, 
     MPO const& H, 
	 Real t,
     ITensor const& LH, 
	 ITensor const& RH,
     const Sweeps& sweeps, 
     DMRGObserver& obs,
     const Args& args = Args::global())
    {
    LocalMPO PH(H,LH,RH,args);
    Real energy = TDVPWorker(psi,PH,t,sweeps,obs,args);
    return energy;
    }

//
//TDVP with a set of MPOs (lazily summed)
//(H vector is 0-indexed)
//
Real inline
tdvp(MPS& psi, 
     std::vector<MPO> const& Hset,
	 Real t, 
     const Sweeps& sweeps,
     const Args& args = Args::global())
    {
    LocalMPOSet PH(Hset,args);
    Real energy = TDVPWorker(psi,PH,t,sweeps,args);
    return energy;
    }

//
//TDVP with a set of MPOs and a custom DMRGObserver
//(H vector is 0-indexed)
//
Real inline
tdvp(MPS & psi, 
     std::vector<MPO> const& Hset, 
	 Real t,
     const Sweeps& sweeps, 
     DMRGObserver& obs,
     const Args& args = Args::global())
    {
    LocalMPOSet PH(Hset,args);
    Real energy = TDVPWorker(psi,PH,t,sweeps,obs,args);
    return energy;
    }


//
// TDVPWorker
//

template <class LocalOpT>
Real
TDVPWorker(MPS & psi,
           LocalOpT& PH,
		   Real t,
           Sweeps const& sweeps,
           Args const& args)
    {
    DMRGObserver obs(psi,args);
    Real energy = TDVPWorker(psi,PH,t,sweeps,obs,args);
    return energy;
    }

template <class LocalOpT>
Real
TDVPWorker(MPS & psi,
           LocalOpT& PH,
		   Real t,
           Sweeps const& sweeps,
           DMRGObserver& obs,
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
    const int numCenter = args.getInt("NumCenter",2);

    const int N = length(psi);
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
            args.add("MinDim",sweeps.mindim(sw));
            args.add("MaxDim",sweeps.maxdim(sw));
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
    
                PH.positionA(b,psi,2);//position2

                auto phi = psi.A(b)*psi.A(b+1);

				energy = krylov(PH,phi,t/2,args);

				if(args.getBool("DoNormalize",true))
					{
					phi/=norm(phi);
					}
				
				auto spec = psi.svdBond(b,phi,(ha==1?Fromleft:Fromright),PH,args);

				if((ha == 1 && b+1 != N) || (ha == 2 && b != 1))
					{	
					PH.positionA((ha == 1? b+1: b),psi,1);//position1: fromleft: b+1,fromright: b
					auto& M = (ha == 1? psi.Aref(b+1):psi.Aref(b));
					energy = krylov(PH,M,-t/2,args);
					if(args.getBool("DoNormalize",true))
						{
						M/=norm(M);
						}
					}

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
    			}
            
			if(obs.checkDone(args)) break;
        
            } //for loop over sw
        }
    else
        {	
        if(numCenter == 1)
            {
			for(int sw = 1; sw <= sweeps.nsweep(); ++sw)
            	{
            	cpu_time sw_time;
            	args.add("Sweep",sw);
            	args.add("NSweep",sweeps.nsweep());
            	args.add("Cutoff",sweeps.cutoff(sw));
            	args.add("MinDim",sweeps.mindim(sw));
            	args.add("MaxDim",sweeps.maxdim(sw));
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
	
            	ITensor M;
				for(int b = 1, ha = 1; ha <= 2; sweepnext1(b,ha,N))
            	    {
            	    if(!quiet)
            	        {
            	        printfln("Sweep=%d, HS=%d, Bond=%d/%d",sw,ha,b,(N-1));
            	        }

            	    PH.positionA(b,psi,1);//position1

            	    ITensor phi;
				  	if((ha == 1 && b != 1) || (ha == 2 && b != 	N))	phi = M*psi.A(b);//const reference
					else	   phi = psi.A(b);

					energy = krylov(PH,phi,t/2,args);
					if(args.getBool("DoNormalize",true))
						{
						phi/=norm(phi);
						}
  	    
					Spectrum spec;
					if((ha == 1 && b != N) || (ha == 2 && b != 1))
						{
						ITensor U,V,S;
						if(ha == 1)	V = ITensor(commonIndex(psi.A(b),psi.A(b+1),"Link"));
						else		V = ITensor(commonIndex(psi.A(b-1),psi.A(b),"Link"));
            	    	spec = svd(phi,U,S,V,args);// QR, oc tensor to be returned
						psi.Aref(b) = U;
						M = S*V;

						PH.positionA((ha == 1? b+1: b),psi,0);//position0

						energy = krylov(PH,M,-t/2,args);
						if(args.getBool("DoNormalize",true))
							{
							M/=norm(M);
							}
						}
					else
						{
						psi.Aref(b) = phi;
						}

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
    				}
            	
				if(obs.checkDone(args)) break;
				} //for loop over sw
			}
        else
            Error("Only support 1 and 2 sites algorithm presently.");
        }
    
	if(args.getBool("DoNormalize",true))
		{
		psi.normalize();
		}

    return energy;
    }

} //namespace itensor


#endif
