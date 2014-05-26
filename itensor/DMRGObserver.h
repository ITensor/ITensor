//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_DMRGOBSERVER_H
#define __ITENSOR_DMRGOBSERVER_H
#include "observer.h"

namespace itensor {

//
// Class for monitoring DMRG calculations.
// The measure and checkDone methods are virtual
// so that behavior can be customized in a
// derived class.
//

template<class Tensor>
class DMRGObserver : public Observer
    {
    public:
    
    DMRGObserver(const MPSt<Tensor>& psi, 
                 const OptSet& opts = Global::opts());

    virtual ~DMRGObserver() { }

    void virtual
    measure(const OptSet& opts = Global::opts());
    
    bool virtual
    checkDone(const OptSet& opts = Global::opts());

    void virtual
    lastSpectrum(const Spectrum& spec) { last_spec_ = spec; }

    const MPSt<Tensor>& 
    psi() const { return psi_; }
    
    private:

    /////////////
    //
    // Data Members

    const MPSt<Tensor>& psi_;

    Vector center_eigs;
    Real energy_errgoal; //Stop DMRG once energy has converged to this precision
    bool printeigs;      //Print slowest decaying eigenvalues after every sweep
    int max_eigs;
    Real max_te;
    bool done_;
    Real last_energy_;
    Spectrum last_spec_;

    Model::DefaultOpsT default_ops_;

    //
    /////////////

    }; // class DMRGObserver

template<class Tensor>
inline DMRGObserver<Tensor>::
DMRGObserver(const MPSt<Tensor>& psi, const OptSet& opts) 
    : 
    psi_(psi),
    energy_errgoal(opts.getReal("EnergyErrgoal",-1)), 
    printeigs(opts.getBool("PrintEigs",true)),
    max_eigs(-1),
    max_te(-1),
    done_(false),
    last_energy_(1000),
    default_ops_(psi.sites().defaultOps())
    { 
    }

template<class Tensor>
void inline DMRGObserver<Tensor>::
measure(const OptSet& opts)
    {
    const int N = psi_.N();
    const int b = opts.getInt("AtBond",1);
    const int sw = opts.getInt("Sweep",0);
    const int ha = opts.getInt("HalfSweep",0);
    const Real energy = opts.getReal("Energy",0);

    if(!opts.getBool("Quiet",false) && !opts.getBool("NoMeasure",false))
        {
        if(b < N && b > 0)
            {
            const Tensor wfb = psi_.A(b)*psi_.A(b+1);
            Foreach(const std::string& opname, default_ops_)
                {
                Complex z = 
                    BraKet(primed(wfb,psi_.sites()(b)),psi_.sites().op(opname,b)*wfb);
                if(fabs(z.imag()) < 1E-14)
                    printfln("<%s>(%d) = %.10E",opname,b,z.real());
                else
                    printfln("<%s>(%d) = (%.10E,%.10E)",opname,b,z.real(),z.imag());
                }
            }
        }

    if(printeigs)
        {
        if(b == N/2 && ha == 2)
            {
            println();
            Vector center_eigs = last_spec_.eigsKept();
            Real S = 0;
            for(int j = 1; j <= center_eigs.Length(); ++j) 
                {
                S -= center_eigs(j)*log(fabs(center_eigs(j)));
                }
            printfln("    vN Entropy at center bond b=%d = %.12f",N/2,S);
            printf("    Eigs at center bond b=%d: ",N/2);
            for(int j = 1; j <= min(center_eigs.Length(),10); ++j) 
                {
                const Real eig = center_eigs(j);
                if(eig < 1E-3) break;
                printf("%.4f ",eig);
                }
            println();
            }
        }

    max_eigs = max(max_eigs,last_spec_.numEigsKept());
    max_te = max(max_te,last_spec_.truncerr());
    if(b == 1 && ha == 2) 
        {
        if(!printeigs) println();
        println("    Largest m during sweep ",sw," was ",(max_eigs > 1 ? max_eigs : 1));
        max_eigs = -1;
        println("    Largest truncation error: ",(max_te > 0 ? max_te : 0.));
        max_te = -1;
        printfln("    Energy after sweep %d is %.12f",sw,energy);
        }

    }


template<class Tensor>
bool inline DMRGObserver<Tensor>::
checkDone(const OptSet& opts)
    {
    const int sw = opts.getInt("Sweep",0);
    const Real energy = opts.getReal("Energy",0);
    
    if(sw == 1) last_energy_ = 1000;
    if(energy_errgoal > 0 && sw%2 == 0)
        {
        Real dE = fabs(energy-last_energy_);
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
        system("rm -f STOP_DMRG");
        return true;
        }

    //Set done_ flag to true so any outer callers using this Observer will also terminate.
    if(fileExists("STOP_DMRG_ALL"))
        {
        println("File STOP_DMRG_ALL found: stopping this run after sweep ",sw);
        system("rm -f STOP_DMRG_ALL");
        done_ = true;
        return done_;
        }
    
    return done_;
    }

}; //namespace itensor

#endif // __ITENSOR_DMRGOBSERVER_H
