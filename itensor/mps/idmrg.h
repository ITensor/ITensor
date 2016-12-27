//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IDMRG_H
#define __ITENSOR_IDMRG_H

#include "itensor/mps/dmrg.h"

namespace itensor {

template <class Tensor>
struct idmrgRVal
    {
    Real energy;
    Tensor HL;
    Tensor HR;
    Tensor IL;
    Tensor V;
    };

template<typename Tensor>
void
read(std::istream & s, idmrgRVal<Tensor> & rval);

template<typename Tensor>
void
write(std::ostream & s, idmrgRVal<Tensor> const& rval);

template <class Tensor>
idmrgRVal<Tensor>
idmrg(MPSt<Tensor>      & psi, 
      MPOt<Tensor> const& H, 
      Sweeps       const& sweeps, 
      Args         const& args = Global::args());

template <class Tensor>
idmrgRVal<Tensor>
idmrg(MPSt<Tensor>      & psi, 
      MPOt<Tensor> const& H,
      Sweeps       const& sweeps,
      DMRGObserver<Tensor> & obs,
      Args         const& args = Global::args());

//For restarting idmrg calculations
//from a previous run (creates a new DMRGObserver automatically)
template <class Tensor>
idmrgRVal<Tensor>
idmrg(MPSt<Tensor>      & psi, 
      MPOt<Tensor> const& H, 
      idmrgRVal<Tensor> const& last_res,
      Sweeps       const& sweeps, 
      Args         const& args);

template <class Tensor>
idmrgRVal<Tensor>
idmrg(MPSt<Tensor> & psi, 
      MPOt<Tensor> H,
      idmrgRVal<Tensor> last_rval,
      Sweeps const& sweeps,
      DMRGObserver<Tensor> & obs,
      Args args = Global::args());


//Given an MPS (or MPO) A1 A2 A3 | A4 A5 A6,
//modifies it to A4 A5 A6 | A1 A2 A3
//(treating the left and right halves
//each as one unit cell)
//Very efficient: swap method only uses pointer swaps internally.
template <class MPSType>
void
swapUnitCells(MPSType & psi)
    {
    auto Nuc = psi.N()/2;
    for(auto n : range1(Nuc))
        {
        psi.Aref(n).swap(psi.Aref(Nuc+n));
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


template <class Tensor>
idmrgRVal<Tensor>
idmrg(MPSt<Tensor> & psi, 
      MPOt<Tensor> H,        //Copies H since algorithm swaps tensors in-place
      idmrgRVal<Tensor> last_rval,
      Sweeps const& sweeps,
      DMRGObserver<Tensor> & obs,
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

    int N0 = psi.N(); //Number of sites in center
    int Nuc = N0/2;   //Number of sites in unit cell
    int N = N0;       //Current system size

    if(N0 == 2) args.add("CombineMPO",false);

    Real energy = NAN;

    auto lastV = last_rval.V;
    Tensor D;

    if(psi.A(0))
        {
        lastV = dag(psi.A(0));
        lastV /= norm(lastV);
        lastV.apply(detail::PseudoInvert(0));
        }

    Tensor HL(last_rval.HL),
           HR(last_rval.HR);

    auto IL = last_rval.IL;

    //If last_rval is trivial,
    //get edge tensors from MPO
    if(not HL) HL = H.A(0);
    if(not HR) HR = H.A(N0+1);

    int sw = 1;

    //Start with two unit cells
        { 
        if(!quiet)
            {
            printfln("\niDMRG Step = %d, N=%d sites",sw,N);
            }

        auto ucsweeps = Sweeps(actual_nucsweeps);
        ucsweeps.minm() = sweeps.minm(sw);
        ucsweeps.maxm() = sweeps.maxm(sw);
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
            for(int j = 1; j <= psi.N(); ++j)
                {
                randomize(psi.Aref(j));
                }
            psi.normalize();
            }

        printfln("\n    Energy per site = %.14f\n",energy/N0);

        psi.position(Nuc);

        args.add("Sweep",sw);
        args.add("AtBond",Nuc);
        args.add("Energy",energy);
        obs.measure(args+Args("AtCenter",true,"NoMeasure",true));

        svd(psi.A(Nuc)*psi.A(Nuc+1),psi.Aref(Nuc),D,psi.Aref(Nuc+1));
        D /= norm(D);
        
        //Prepare MPO for next step
        for(int j = 1; j <= Nuc; ++j)
            {
            HL *= psi.A(j);
            HL *= H.A(j);
            HL *= dag(prime(psi.A(j)));
            IL *= psi.A(j);
            IL *= H.A(j);
            IL *= dag(prime(psi.A(j)));

            HR *= psi.A(N0-j+1);
            HR *= H.A(N0-j+1);
            HR *= dag(prime(psi.A(N0-j+1)));
            }
        //H = HG(sw);
        swapUnitCells(H);

        HL += -energy*IL;

        //Prepare MPS for next step
        swapUnitCells(psi);
        if(lastV) psi.Aref(Nuc+1) *= lastV;
        psi.Aref(1) *= D;
        psi.Aref(N0) *= D;
        psi.position(1);

        ++sw;
        }

    Spectrum spec;

    for(; sw <= sweeps.nsweep(); ++sw)
        {
        auto ucsweeps = Sweeps(actual_nucsweeps);
        ucsweeps.minm() = sweeps.minm(sw);
        ucsweeps.maxm() = sweeps.maxm(sw);
        ucsweeps.cutoff() = sweeps.cutoff(sw);
        ucsweeps.noise() = sweeps.noise(sw);
        ucsweeps.niter() = sweeps.niter(sw);
        args.add("Maxm",sweeps.maxm(sw));

        print(ucsweeps);

        if(actual_nucsweeps > 1) actual_nucsweeps -= nucs_decr;

        N += N0;

        if(!quiet) printfln("\niDMRG Step = %d, N=%d sites",sw,N);

        auto initPsi = psi;

        auto PH = LocalMPO<Tensor>(H,HL,HR,args);

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
            psiphi(initPsi,psi,ovrlap,im);
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

        D = Tensor();
        svd(psi.A(Nuc)*psi.A(Nuc+1),psi.Aref(Nuc),D,psi.Aref(Nuc+1),args);
        D /= norm(D);

        //Prepare MPO for next step
        for(int j = 1; j <= Nuc; ++j)
            {
            HL *= psi.A(j);
            HL *= H.A(j);
            HL *= dag(prime(psi.A(j)));
            IL *= psi.A(j);
            IL *= H.A(j);
            IL *= dag(prime(psi.A(j)));

            HR *= psi.A(N0-j+1);
            HR *= H.A(N0-j+1);
            HR *= dag(prime(psi.A(N0-j+1)));
            }
        swapUnitCells(H);

        HL += -energy*IL;

        //Prepare MPS for next step
        swapUnitCells(psi);

        psi.Aref(N0) *= D;

        if((obs.checkDone(args) && sw%2==0)
           || sw == sweeps.nsweep()) 
            {
            //Convert A's (left-ortho) to B's by moving D (center matrix)
            //through until last V*A_j*D == B_j
            for(int b = N0-1; b >= Nuc+1; --b)
                {
                Tensor d;
                svd(psi.A(b)*psi.A(b+1),psi.Aref(b),d,psi.Aref(b+1));
                psi.Aref(b) *= d;
                }
            psi.Aref(Nuc+1) *= lastV;

            psi.Aref(0) = D;

            break;
            }

        if(fileExists("WRITE_WF") && sw%2==0)
            {
            println("File WRITE_WF found: writing out wavefunction after step",sw);
            system("rm -f WRITE_WF");
            auto wpsi = psi;
            for(int b = N0-1; b >= Nuc+1; --b)
                {
                Tensor d;
                svd(wpsi.A(b)*wpsi.A(b+1),wpsi.Aref(b),d,wpsi.Aref(b+1));
                wpsi.Aref(b) *= d;
                }
            wpsi.Aref(Nuc+1) *= lastV;
            wpsi.Aref(0) = D;
            writeToFile(format("psi_%d",sw),wpsi);
            writeToFile("sites",wpsi.sites());
            }

        psi.Aref(Nuc+1) *= lastV;
        psi.Aref(1) *= D;

        psi.orthogonalize();
        psi.normalize();

        } //for loop over sw
    
    auto res = idmrgRVal<Tensor>();
    res.energy = energy;
    res.HL = HL;
    res.HR = HR;
    res.IL = IL;
    res.V = lastV;

    return res;
    }

template <class Tensor>
idmrgRVal<Tensor>
idmrg(MPSt<Tensor>      & psi, 
      MPOt<Tensor> const& H,
      Sweeps       const& sweeps,
      DMRGObserver<Tensor> & obs,
      Args         const& args)
    {
    //Assumes H.A(N+1) contains vector
    //picking out ending state of MPO
    //automaton:
    auto lval = idmrgRVal<Tensor>();
    lval.IL = Tensor(dag(H.A(H.N()+1)));
    return idmrg(psi,H,lval,sweeps,obs,args);
    }

template <class Tensor>
idmrgRVal<Tensor>
idmrg(MPSt<Tensor>      & psi, 
      MPOt<Tensor> const& H, 
      Sweeps       const& sweeps, 
      Args         const& args)
    {
    auto obs = DMRGObserver<Tensor>(psi);
    return idmrg(psi,H,sweeps,obs,args);
    }

template <class Tensor>
idmrgRVal<Tensor>
idmrg(MPSt<Tensor>      & psi, 
      MPOt<Tensor> const& H, 
      idmrgRVal<Tensor> const& last_res,
      Sweeps       const& sweeps, 
      Args         const& args)
    {
    auto obs = DMRGObserver<Tensor>(psi);
    return idmrg(psi,H,last_res,sweeps,obs,args);
    }

template<typename Tensor>
void
read(std::istream & s, idmrgRVal<Tensor> & R)
    {
    itensor::read(s,R.energy);
    itensor::read(s,R.HL);
    itensor::read(s,R.HR);
    itensor::read(s,R.IL);
    itensor::read(s,R.V);
    }

template<typename Tensor>
void
write(std::ostream & s, idmrgRVal<Tensor> const& R)
    {
    itensor::write(s,R.energy);
    itensor::write(s,R.HL);
    itensor::write(s,R.HR);
    itensor::write(s,R.IL);
    itensor::write(s,R.V);
    }

} //namespace itensor


#endif
