//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IDMRG_H
#define __ITENSOR_IDMRG_H

#include "itensor/mps/dmrg.h"

namespace itensor {

struct idmrgRVal
    {
    Real energy;
    ITensor HL;
    ITensor HR;
    ITensor IL;
    ITensor V;
    };

void
read(std::istream & s, idmrgRVal & rval);

void
write(std::ostream & s, idmrgRVal const& rval);

idmrgRVal
idmrg(MPS      & psi, 
      MPO const& H, 
      Sweeps       const& sweeps, 
      Args         const& args = Args::global());

idmrgRVal
idmrg(MPS      & psi, 
      MPO const& H,
      Sweeps       const& sweeps,
      DMRGObserver & obs,
      Args         const& args = Args::global());

//For restarting idmrg calculations
//from a previous run (creates a new DMRGObserver automatically)
idmrgRVal
idmrg(MPS      & psi, 
      MPO const& H, 
      idmrgRVal const& last_res,
      Sweeps       const& sweeps, 
      Args         const& args);

idmrgRVal
idmrg(MPS & psi, 
      MPO H,
      idmrgRVal last_rval,
      Sweeps const& sweeps,
      DMRGObserver & obs,
      Args args = Args::global());


//Given an MPS (or MPO) A1 A2 A3 | A4 A5 A6,
//modifies it to A4 A5 A6 | A1 A2 A3
//(treating the left and right halves
//each as one unit cell)
//Very efficient: swap method only uses pointer swaps internally.
template <class MPSType>
void
swapUnitCells(MPSType & psi)
    {
    auto Nuc = length(psi)/2;
    for(auto n : range1(Nuc))
        {
        psi.ref(n).swap(psi.ref(Nuc+n));
        }
    }

//
// Implementations
//


namespace detail {

struct PseudoInvert
    {
    Real cutoff = 0.;
    PseudoInvert(Real cut) : cutoff(cut) { }

    Real
    operator()(Real x) const
        {
        return (x > cutoff) ? 1./x : 0.;
        }
    };

} //namespace detail


idmrgRVal inline
idmrg(MPS & psi, 
      MPO H,        //Copies H since algorithm swaps tensors in-place
      idmrgRVal last_rval,
      Sweeps const& sweeps,
      DMRGObserver & obs,
      Args args)
    {
    auto olevel = args.getInt("OutputLevel",0);
    auto quiet = args.getBool("Quiet",olevel == 0);
    auto nucsweeps = args.getInt("NUCSweeps",1);
    auto nucs_decr = args.getInt("NUCSweepsDecrement",0);
    auto do_randomize = args.getBool("Randomize",false);
    auto show_overlap = args.getBool("ShowOverlap",false);
    //inverse_cut is cutoff for computing pseudo inverse 
    auto inverse_cut = args.getReal("InverseCut",1E-8);
    auto actual_nucsweeps = nucsweeps;

    int N0 = length(psi); //Number of sites in center
    int Nuc = N0/2;   //Number of sites in unit cell
    int N = N0;       //Current system size

    if(N0 == 2) args.add("CombineMPO",false);

    Real energy = NAN;

    auto lastV = last_rval.V;
    ITensor D;

    if(psi(0))
        {
        lastV = dag(psi(0));
        lastV /= norm(lastV);
        lastV.apply(detail::PseudoInvert(0));
        }

    ITensor HL(last_rval.HL),
            HR(last_rval.HR);

    auto IL = last_rval.IL;

    //If last_rval is trivial,
    //get edge tensors from MPO
    if(not HL) HL = H(0);
    if(not HR) HR = H(N0+1);

    int sw = 1;

    //Start with two unit cells
        { 
        if(!quiet)
            {
            printfln("\niDMRG Step = %d, N=%d sites",sw,N);
            }

        auto ucsweeps = Sweeps(actual_nucsweeps);
        ucsweeps.mindim() = sweeps.mindim(sw);
        ucsweeps.maxdim() = sweeps.maxdim(sw);
        ucsweeps.cutoff() = sweeps.cutoff(sw);
        ucsweeps.noise() = sweeps.noise(sw);
        ucsweeps.niter() = sweeps.niter(sw);
        print(ucsweeps);

        auto extra_args = Args("Quiet",olevel < 2,
                               "iDMRG_Step",sw,
                               "NSweep",ucsweeps.nsweep());
        energy = dmrg(psi,H,HL,HR,ucsweeps,obs,args + extra_args);

        if(do_randomize)
            {
            println("Randomizing psi");
            for(int j = 1; j <= length(psi); ++j)
                {
                psi.ref(j).randomize();
                }
            psi.normalize();
            }

        printfln("\n    Energy per site = %.14f\n",energy/N0);

        psi.position(Nuc);

        args.add("Sweep",sw);
        args.add("AtBond",Nuc);
        args.add("Energy",energy);
        obs.measure(args+Args("AtCenter",true,"NoMeasure",true));

        svd(psi(Nuc)*psi(Nuc+1),psi.ref(Nuc),D,psi.ref(Nuc+1));
        D /= norm(D);
        
        //Prepare MPO for next step
        for(int j = 1; j <= Nuc; ++j)
            {
            HL *= psi(j);
            HL *= H(j);
            HL *= dag(prime(psi(j)));
            IL *= psi(j);
            IL *= H(j);
            IL *= dag(prime(psi(j)));

            HR *= psi(N0-j+1);
            HR *= H(N0-j+1);
            HR *= dag(prime(psi(N0-j+1)));
            }
        //H = HG(sw);
        swapUnitCells(H);

        HL += -energy*IL;

        //Prepare MPS for next step
        swapUnitCells(psi);
        if(lastV) psi.ref(Nuc+1) *= lastV;
        psi.ref(1) *= D;
        psi.ref(N0) *= D;
        psi.position(1);

        ++sw;
        }

    Spectrum spec;

    for(; sw <= sweeps.nsweep(); ++sw)
        {
        auto ucsweeps = Sweeps(actual_nucsweeps);
        ucsweeps.mindim() = sweeps.mindim(sw);
        ucsweeps.maxdim() = sweeps.maxdim(sw);
        ucsweeps.cutoff() = sweeps.cutoff(sw);
        ucsweeps.noise() = sweeps.noise(sw);
        ucsweeps.niter() = sweeps.niter(sw);
        args.add("MaxDim",sweeps.maxdim(sw));

        print(ucsweeps);

        if(actual_nucsweeps > 1) actual_nucsweeps -= nucs_decr;

        N += N0;

        if(!quiet) printfln("\niDMRG Step = %d, N=%d sites",sw,N);

        auto initPsi = psi;

        auto PH = LocalMPO(H,HL,HR,args);

        if(olevel >= 1)
            {
            auto ien = overlap(psi,H,HL,HR,psi);
            printfln("Initial energy = %.20f",ien);
            }
        
        auto extra_args = Args("Quiet",olevel<2,"NoMeasure",sw%2==0,"iDMRG_Step",sw,"NSweep",ucsweeps.nsweep());
        energy = DMRGWorker(psi,PH,ucsweeps,obs,args + extra_args);

        if(show_overlap || olevel >= 1)
            {
            Real ovrlap, im;
            overlap(initPsi,psi,ovrlap,im);
            print("\n    Overlap of initial and final psi = ");
            printfln((std::fabs(ovrlap) > 1E-4 ? "%.10f" : "%.10E"),std::fabs(ovrlap));
            print("\n    1-Overlap of initial and final psi = ");
            printfln((1-std::fabs(ovrlap) > 1E-4 ? "%.10f" : "%.10E"),1-std::fabs(ovrlap));
            }

        printfln("    Energy per site = %.14f",energy/N0);


        //Save last center matrix
        lastV = dag(D);
        lastV /= norm(lastV);
        lastV.apply(detail::PseudoInvert(inverse_cut));

        //Calculate new center matrix
        psi.position(Nuc);

        args.add("Sweep",sw);
        args.add("AtBond",Nuc);
        args.add("Energy",energy);
        obs.measure(args+Args("AtCenter",true,"NoMeasure",true));

        D = ITensor();
        svd(psi(Nuc)*psi(Nuc+1),psi.ref(Nuc),D,psi.ref(Nuc+1),args);
        D /= norm(D);

        //Prepare MPO for next step
        for(int j = 1; j <= Nuc; ++j)
            {
            HL *= psi(j);
            HL *= H(j);
            HL *= dag(prime(psi(j)));
            IL *= psi(j);
            IL *= H(j);
            IL *= dag(prime(psi(j)));

            HR *= psi(N0-j+1);
            HR *= H(N0-j+1);
            HR *= dag(prime(psi(N0-j+1)));
            }
        swapUnitCells(H);

        HL += -energy*IL;

        //Prepare MPS for next step
        swapUnitCells(psi);

        psi.ref(N0) *= D;

        if((obs.checkDone(args) && sw%2==0)
           || sw == sweeps.nsweep()) 
            {
            //Convert A's (left-ortho) to B's by moving D (center matrix)
            //through until last V*A_j*D == B_j
            for(int b = N0-1; b >= Nuc+1; --b)
                {
                ITensor d;
                svd(psi(b)*psi(b+1),psi.ref(b),d,psi.ref(b+1));
                psi.ref(b) *= d;
                }
            psi.ref(Nuc+1) *= lastV;

            psi.ref(0) = D;

            break;
            }

        if(fileExists("WRITE_WF") && sw%2==0)
            {
            println("File WRITE_WF found: writing out wavefunction after step",sw);
            system("rm -f WRITE_WF");
            auto wpsi = psi;
            for(int b = N0-1; b >= Nuc+1; --b)
                {
                ITensor d;
                svd(wpsi(b)*wpsi(b+1),wpsi.ref(b),d,wpsi.ref(b+1));
                wpsi.ref(b) *= d;
                }
            wpsi.ref(Nuc+1) *= lastV;
            wpsi.ref(0) = D;
            writeToFile(format("psi_%d",sw),wpsi);
            //writeToFile("sites",wpsi.sites());
            }

        psi.ref(Nuc+1) *= lastV;
        psi.ref(1) *= D;

        psi.orthogonalize();
        psi.normalize();

        } //for loop over sw
    
    auto res = idmrgRVal();
    res.energy = energy;
    res.HL = HL;
    res.HR = HR;
    res.IL = IL;
    res.V = lastV;

    return res;
    }

idmrgRVal inline
idmrg(MPS      & psi, 
      MPO const& H,
      Sweeps       const& sweeps,
      DMRGObserver & obs,
      Args         const& args)
    {
    //Assumes H(N+1) contains vector
    //picking out ending state of MPO
    //automaton:
    auto lval = idmrgRVal();
    lval.IL = ITensor(dag(H(length(H)+1)));
    return idmrg(psi,H,lval,sweeps,obs,args);
    }

idmrgRVal inline
idmrg(MPS      & psi, 
      MPO const& H, 
      Sweeps       const& sweeps, 
      Args         const& args)
    {
    auto obs = DMRGObserver(psi);
    return idmrg(psi,H,sweeps,obs,args);
    }

idmrgRVal inline
idmrg(MPS      & psi, 
      MPO const& H, 
      idmrgRVal const& last_res,
      Sweeps       const& sweeps, 
      Args         const& args)
    {
    auto obs = DMRGObserver(psi);
    return idmrg(psi,H,last_res,sweeps,obs,args);
    }

void inline
read(std::istream & s, idmrgRVal & R)
    {
    itensor::read(s,R.energy);
    itensor::read(s,R.HL);
    itensor::read(s,R.HR);
    itensor::read(s,R.IL);
    itensor::read(s,R.V);
    }

void inline
write(std::ostream & s, idmrgRVal const& R)
    {
    itensor::write(s,R.energy);
    itensor::write(s,R.HL);
    itensor::write(s,R.HR);
    itensor::write(s,R.IL);
    itensor::write(s,R.V);
    }

} //namespace itensor


#endif
