//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IDMRG_H
#define __ITENSOR_IDMRG_H

#include "dmrg.h"

#define Cout std::cout
#define Endl std::endl
#define Format boost::format

template <class Tensor>
Real
idmrg(MPSt<Tensor>& psi, 
      const MPOt<Tensor>& H, 
      const Sweeps& sweeps, 
      const OptSet& opts = Global::opts());

template <class Tensor>
Real
idmrg(MPSt<Tensor>& psi, 
      MPOt<Tensor> H,
      const Sweeps& sweeps,
      Observer& obs,
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
Real
idmrg(MPSt<Tensor>& psi, 
      MPOt<Tensor> H,        //Copies H since algorithm swaps tensors in-place
      const Sweeps& sweeps,
      Observer& obs,
      OptSet opts)
    {
    typedef typename Tensor::IndexT
    IndexT;

    const int olevel = opts.getInt("OutputLevel",0);
    const bool quiet = opts.getBool("Quiet",olevel == 0);
    const bool measure_xi = opts.getBool("MeasureCorrLen",false);

    const Real orig_cutoff = psi.cutoff(),
               orig_noise  = psi.noise();
    const int orig_minm = psi.minm(), 
              orig_maxm = psi.maxm();

    const int N0 = psi.N(); //Number of sites in center
    const int Nuc = N0/2;   //Number of sites in unit cell
    int N = N0;             //Current system size

    if(N0 == 2)
        opts.add(Opt("CombineMPO",false));

    Real energy = NAN,
         lastenergy = 0;

    Tensor lastV,
           D;

    const int Nteig = 2;
    std::vector<Tensor> vv(Nteig);
    std::vector<Real>   ee(Nteig,1.);

    Tensor HL(H.A(0)),
           HR(H.A(N0+1));

    int sw = 1;

    //Start with two unit cells
        { 
        if(!quiet)
            {
            Cout << Format("\niDMRG Step = %d, N=%d sites") % sw % N << Endl;
            }

        Sweeps ucsweeps(sweeps.niter(sw));
        ucsweeps.minm() = sweeps.minm(sw);
        ucsweeps.maxm() = sweeps.maxm(sw);
        ucsweeps.cutoff() = sweeps.cutoff(sw);
        ucsweeps.noise() = sweeps.noise(sw);
        ucsweeps.niter() = 2;
        Cout << ucsweeps;

        const Real fac = sqrt(1./N0);
        HL *= fac;
        HR *= fac;

        energy = dmrg(psi,H,HL,HR,ucsweeps,obs,opts & Quiet(olevel < 3));

        psi.position(Nuc);
        svd(psi.A(Nuc)*psi.A(Nuc+1),psi.Anc(Nuc),D,psi.Anc(Nuc+1));
        D /= D.norm();

        //Prepare MPO for next step
        for(int j = 1; j <= Nuc; ++j)
            {
            HL *= psi.A(j);
            HL *= H.A(j);
            HL *= conj(primed(psi.A(j)));

            HR *= psi.A(N0-j+1);
            HR *= H.A(N0-j+1);
            HR *= conj(primed(psi.A(N0-j+1)));
            }
        swapUnitCells(H);

        if(measure_xi)
            {
            vv[0] = psi.A(1)*conj(primed(psi.A(1),Link));
            for(int j = 2; j <= Nuc; ++j)
                {
                vv[0] *= psi.A(j);
                vv[0] *= conj(primed(psi.A(j),Link));
                }
            }

        //Prepare MPS for next step
        swapUnitCells(psi);
        psi.Anc(1) *= D;
        psi.Anc(N0) *= D;
        psi.position(1);

        ++sw;
        }


    //Orthonormalize vv tensors
    if(measure_xi)
        {
        vv[0] /= vv[0].norm();
        for(int j = 1; j < Nteig; ++j)
            {
            vv[j] = vv[0];
            vv[j].randomize();
            for(int k = 0; k < j; ++k)
                {
                vv[j] += vv[k]*(-BraKet(vv[k],vv[j]));
                }
            vv[j] /= vv[j].norm();
            }
        }


    Real sub_en_per_site = energy/N0;

    Spectrum spec;

    for(; sw <= sweeps.nsweep(); ++sw)
        {
        Sweeps ucsweeps(sweeps.niter(sw));
        ucsweeps.minm() = sweeps.minm(sw);
        ucsweeps.maxm() = sweeps.maxm(sw);
        ucsweeps.cutoff() = sweeps.cutoff(sw);
        ucsweeps.noise() = sweeps.noise(sw);
        ucsweeps.niter() = 2;
        spec.maxm(sweeps.maxm(sw));

        Cout << ucsweeps;

        N += N0;

        if(!quiet)
            {
            Cout << Format("\niDMRG Step = %d, N=%d sites") % sw % N << Endl;
            }

        const Real fac = sqrt((1.*N-N0)/N);
        HL *= fac;
        HR *= fac;

        const MPSt<Tensor> initPsi(psi);

        lastenergy = energy;
        LocalMPO<Tensor> PH(H,HL,HR,opts);
        energy = DMRGWorker(psi,PH,ucsweeps,obs,opts & Quiet(olevel < 3) & Opt("NoMeasure",sw%2==0));


        Real ovrlap, im;
        psiphi(initPsi,psi,ovrlap,im);
        Cout << "\n    Overlap of initial and final psi = " << Format(fabs(ovrlap) > 1E-4 ? "%.10f" : "%.10E") % fabs(ovrlap) << Endl;
        Cout << "\n    1-Overlap of initial and final psi = " << Format(1-fabs(ovrlap) > 1E-4 ? "%.10f" : "%.10E") % (1-fabs(ovrlap)) << Endl;


        sub_en_per_site = (energy*N-lastenergy*(N-N0))/N0;

        Cout << Format("    Energy per site = %.10f\n") % energy << Endl;
        Cout << Format("    Subtracted Energy per site = %.14f\n") % sub_en_per_site << Endl;


        //Save last center matrix
        lastV = conj(D);
        lastV /= lastV.norm();
        lastV.pseudoInvert(0);

        //Calculate new center matrix
        psi.position(Nuc);
        obs.measure(N,sw,1,Nuc,spec,sub_en_per_site,opts&Opt("AtCenter")&Opt("NoMeasure"));

        D = Tensor();
        svd(psi.A(Nuc)*psi.A(Nuc+1),psi.Anc(Nuc),D,psi.Anc(Nuc+1),spec);
        D /= D.norm();


        //Prepare MPO for next step

        PH.position(Nuc,psi);

        HR = PH.R();
        HR *= psi.A(Nuc+1);
        HR *= H.A(Nuc+1);
        HR *= conj(primed(psi.A(Nuc+1)));

        HL = PH.L();
        HL *= psi.A(Nuc);
        HL *= H.A(Nuc);
        HL *= conj(primed(psi.A(Nuc)));

        swapUnitCells(H);

        if(measure_xi)
            {
            //Update correlation length estimate
            Real xi = NAN;
            for(int j = 0; j < Nteig; ++j)
                {
                for(int i = 1; i <= Nuc; ++i)
                    {
                    vv[j] *= psi.A(i);
                    vv[j] *= conj(primed(psi.A(i),Link));
                    }
                for(int k = 0; k < j; ++k)
                    {
                    vv[j] += vv[k]*(-BraKet(vv[k],vv[j]));
                    }
                ee[j] = vv[j].norm();
                const Real eig = ee[j];
                if(eig == 0)
                    {
                    Error("Zero norm of transfer matrix eigenvector");
                    }
                vv[j] /= eig;
//#ifdef DEBUG
                Cout << Format("    T eig(%d) = %.14f\n") % (1+j) % eig << Endl;
                if(j == 0 && fabs(eig-1.) > 1E-4)
                    {
                    Cout << Format("    Leading transfer eigenvalue = %.14f\n") % eig << Endl;
                    }
//#endif
                if(j > 0 && eig > 1E-12)
                    {
                    xi = -1.*Nuc/log(eig);
                    }
                }
            Cout << Format("    Correlation length = %.14f\n") % xi << Endl;
            }

        //Prepare MPS for next step
        swapUnitCells(psi);

        psi.Anc(N0) *= D;

        if((obs.checkDone(sw,sub_en_per_site) && sw%2==0)
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

        psi.Anc(Nuc+1) *= lastV;
        psi.Anc(1) *= D;

        psi.orthogonalize();
        psi.normalize();

        } //for loop over sw
    
    psi.cutoff(orig_cutoff); 
    psi.minm(orig_minm); 
    psi.maxm(orig_maxm);
    psi.noise(orig_noise); 

    return sub_en_per_site;
    }

template <class Tensor>
Real
idmrg(MPSt<Tensor>& psi, const MPOt<Tensor>& H, 
      const Sweeps& sweeps, 
      const OptSet& opts)
    {
    DMRGObserver<Tensor> obs(psi);
    return idmrg(psi,H,sweeps,obs,opts);
    }

#undef Cout
#undef Format
#undef Endl

#endif
