//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IDMRG_H
#define __ITENSOR_IDMRG_H

#include "dmrg.h"

namespace itensor {

template <class Tensor>
struct idmrgRVal
    {
    Real energy;
    Tensor HL;
    Tensor HR;
    Tensor V;
    };

template <class Tensor>
idmrgRVal<Tensor>
idmrg(MPSt<Tensor>& psi, 
      const MPOt<Tensor>& H, 
      const Sweeps& sweeps, 
      const OptSet& opts = Global::opts());

template <class Tensor>
idmrgRVal<Tensor>
idmrg(MPSt<Tensor>& psi, 
      const MPOt<Tensor>& H,
      const Sweeps& sweeps,
      DMRGObserver<Tensor>& obs,
      const OptSet& opts = Global::opts());

template <class Tensor>
idmrgRVal<Tensor>
idmrg(MPSt<Tensor>& psi, 
      MPOt<Tensor> H,        //Copies H since algorithm swaps tensors in-place
      Tensor IL,
      const Sweeps& sweeps,
      DMRGObserver<Tensor>& obs,
      OptSet opts = Global::opts());


//Given an MPS (or MPO) A1 A2 A3 | A4 A5 A6,
//modifies it to A4 A5 A6 | A1 A2 A3
//(treating the left and right halves
//each as one unit cell)
//Very efficient: swap method only uses pointer swaps internally.
template <class MPSType>
void
swapUnitCells(MPSType& psi)
    {
    const int Nuc = psi.N()/2;
    for(int n = 1; n <= Nuc; ++n)
        {
        psi.Anc(n).swap(psi.Anc(Nuc+n));
        }
    }

//
// Implementations
//



template <class Tensor>
idmrgRVal<Tensor>
idmrg(MPSt<Tensor>& psi, 
      MPOt<Tensor> H,        //Copies H since algorithm swaps tensors in-place
      Tensor IL,
      const Sweeps& sweeps,
      DMRGObserver<Tensor>& obs,
      OptSet opts)
    {
    typedef typename Tensor::IndexT
    IndexT;

    const int olevel = opts.getInt("OutputLevel",0);
    const bool quiet = opts.getBool("Quiet",olevel == 0);
    const int nucsweeps = opts.getInt("NUCSweeps",1);
    const int nucs_decr = opts.getInt("NUCSweepsDecrement",0);
    const bool randomize = opts.getBool("Randomize",false);
    int actual_nucsweeps = nucsweeps;

    const int N0 = psi.N(); //Number of sites in center
    const int Nuc = N0/2;   //Number of sites in unit cell
    int N = N0;             //Current system size

    if(N0 == 2) opts.add(Opt("CombineMPO",false));

    Real energy = NAN,
         lastenergy = 0;

    Tensor lastV,
           D;

    if(psi.A(0))
        {
        lastV = dag(psi.A(0));
        lastV /= lastV.norm();
        lastV.pseudoInvert(0);
        }

    Tensor HL(H.A(0)),
           HR(H.A(N0+1));

    int sw = 1;

    //Start with two unit cells
        { 
        if(!quiet)
            {
            printfln("\niDMRG Step = %d, N=%d sites",sw,N);
            }

        Sweeps ucsweeps(actual_nucsweeps);
        ucsweeps.minm() = sweeps.minm(sw);
        ucsweeps.maxm() = sweeps.maxm(sw);
        ucsweeps.cutoff() = sweeps.cutoff(sw);
        ucsweeps.noise() = sweeps.noise(sw);
        ucsweeps.niter() = sweeps.niter(sw);
        print(ucsweeps);

        energy = dmrg(psi,H,HL,HR,ucsweeps,obs,opts & Opt("Quiet",olevel < 3)
                                                    & Opt("iDMRG_Step",sw)
                                                    & Opt("NSweep",ucsweeps.nsweep()));

        if(randomize)
            {
            println("Randomizing psi");
            for(int j = 1; j <= psi.N(); ++j)
                {
                psi.Anc(j).randomize();
                }
            psi.normalize();
            }

        printfln("\n    Energy per site = %.14f\n",energy/N0);

        psi.position(Nuc);
        svd(psi.A(Nuc)*psi.A(Nuc+1),psi.Anc(Nuc),D,psi.Anc(Nuc+1));
        D /= D.norm();
        
        //Prepare MPO for next step
        for(int j = 1; j <= Nuc; ++j)
            {
            HL *= psi.A(j);
            HL *= H.A(j);
            HL *= dag(primed(psi.A(j)));
            IL *= psi.A(j);
            IL *= H.A(j);
            IL *= dag(primed(psi.A(j)));

            HR *= psi.A(N0-j+1);
            HR *= H.A(N0-j+1);
            HR *= dag(primed(psi.A(N0-j+1)));
            }
        //H = HG(sw);
        swapUnitCells(H);

        HL += -energy*IL;

        //Prepare MPS for next step
        swapUnitCells(psi);
        if(lastV) psi.Anc(Nuc+1) *= lastV;
        psi.Anc(1) *= D;
        psi.Anc(N0) *= D;
        psi.position(1);

        ++sw;
        }

    Spectrum spec;

    for(; sw <= sweeps.nsweep(); ++sw)
        {
        Sweeps ucsweeps(actual_nucsweeps);
        ucsweeps.minm() = sweeps.minm(sw);
        ucsweeps.maxm() = sweeps.maxm(sw);
        ucsweeps.cutoff() = sweeps.cutoff(sw);
        ucsweeps.noise() = sweeps.noise(sw);
        ucsweeps.niter() = sweeps.niter(sw);
        opts.add("Maxm",sweeps.maxm(sw));

        print(ucsweeps);

        if(actual_nucsweeps > 1) actual_nucsweeps -= nucs_decr;

        N += N0;

        if(!quiet)
            {
            printfln("\niDMRG Step = %d, N=%d sites",sw,N);
            }

        const MPSt<Tensor> initPsi(psi);

        lastenergy = energy;
        LocalMPO<Tensor> PH(H,HL,HR,opts);
        
        energy = DMRGWorker(psi,PH,ucsweeps,obs,opts & Opt("Quiet",olevel < 3) 
                                                     & Opt("NoMeasure",sw%2==0)
                                                     & Opt("iDMRG_Step",sw)
                                                     & Opt("NSweep",ucsweeps.nsweep()));

        Real ovrlap, im;
        psiphi(initPsi,psi,ovrlap,im);
        print("\n    Overlap of initial and final psi = ");
        printfln((fabs(ovrlap) > 1E-4 ? "%.10f" : "%.10E"),fabs(ovrlap));
        print("\n    1-Overlap of initial and final psi = ");
        printfln((1-fabs(ovrlap) > 1E-4 ? "%.10f" : "%.10E"),1-fabs(ovrlap));

        printfln("    Energy per site = %.14f",energy/N0);


        //Save last center matrix
        lastV = dag(D);
        lastV /= lastV.norm();
        lastV.pseudoInvert(0);

        //Calculate new center matrix
        psi.position(Nuc);

        opts.add("Sweep",sw);
        opts.add("AtBond",Nuc);
        opts.add("Energy",energy);
        obs.measure(opts&Opt("AtCenter")&Opt("NoMeasure"));

        D = Tensor();
        svd(psi.A(Nuc)*psi.A(Nuc+1),psi.Anc(Nuc),D,psi.Anc(Nuc+1),opts);
        D /= D.norm();


        //Prepare MPO for next step
        for(int j = 1; j <= Nuc; ++j)
            {
            HL *= psi.A(j);
            HL *= H.A(j);
            HL *= dag(primed(psi.A(j)));
            IL *= psi.A(j);
            IL *= H.A(j);
            IL *= dag(primed(psi.A(j)));

            HR *= psi.A(N0-j+1);
            HR *= H.A(N0-j+1);
            HR *= dag(primed(psi.A(N0-j+1)));
            }
        swapUnitCells(H);

        HL += -energy*IL;

        //Prepare MPS for next step
        swapUnitCells(psi);

        psi.Anc(N0) *= D;

        if((obs.checkDone(opts) && sw%2==0)
           || sw == sweeps.nsweep()) 
            {
            //Convert A's (left-ortho) to B's by moving D (center matrix)
            //through until last V*A_j*D == B_j
            for(int b = N0-1; b >= Nuc+1; --b)
                {
                Tensor d;
                svd(psi.A(b)*psi.A(b+1),psi.Anc(b),d,psi.Anc(b+1));
                psi.Anc(b) *= d;
                }
            psi.Anc(Nuc+1) *= lastV;

            psi.Anc(0) = D;

            break;
            }

        if(fileExists("WRITE_WF") && sw%2==0)
            {
            println("File WRITE_WF found: writing out wavefunction after step",sw);
            system("rm -f WRITE_WF");
            MPSt<Tensor> wpsi(psi);
            for(int b = N0-1; b >= Nuc+1; --b)
                {
                Tensor d;
                svd(wpsi.A(b)*wpsi.A(b+1),wpsi.Anc(b),d,wpsi.Anc(b+1));
                wpsi.Anc(b) *= d;
                }
            wpsi.Anc(Nuc+1) *= lastV;
            wpsi.Anc(0) = D;
            writeToFile(format("psi_%d",sw),wpsi);
            writeToFile("sites",wpsi.model());
            }

        psi.Anc(Nuc+1) *= lastV;
        psi.Anc(1) *= D;

        psi.orthogonalize();
        psi.normalize();

        } //for loop over sw
    
    idmrgRVal<Tensor> res;
    res.energy = energy;
    res.HL = HL;
    res.HR = HR;
    res.V = lastV;

    return res;
    }

template <class Tensor>
idmrgRVal<Tensor>
idmrg(MPSt<Tensor>& psi, 
      const MPOt<Tensor>& H,
      const Sweeps& sweeps,
      DMRGObserver<Tensor>& obs,
      const OptSet& opts)
    {
    Tensor IL(dag(H.A(H.N()+1)));
    return idmrg(psi,H,IL,sweeps,obs,opts);
    }

template <class Tensor>
idmrgRVal<Tensor>
idmrg(MPSt<Tensor>& psi, 
      const MPOt<Tensor>& H, 
      const Sweeps& sweeps, 
      const OptSet& opts)
    {
    DMRGObserver<Tensor> obs(psi);
    return idmrg(psi,H,sweeps,obs,opts);
    }

}; //namespace itensor


#endif
