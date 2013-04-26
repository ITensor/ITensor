//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_DMRGOBSERVER_H
#define __ITENSOR_DMRGOBSERVER_H
#include "observer.h"
#include "svdworker.h"

#define Cout std::cout
#define Endl std::endl
#define Format boost::format

//
// Class for monitoring DMRG calculations.
// The measure and checkDone methods are virtual
// so that behavior can be customized in a
// derived class.
//

class DMRGObserver : public Observer
    {
    public:
    
    DMRGObserver();

    virtual ~DMRGObserver() { }

    void virtual
    measure(int sw, int ha, int b, const SVDWorker& svd, Real energy,
            const OptSet& opts = Global::opts());
    
    bool virtual
    checkDone(int sw, const SVDWorker& svd, Real energy,
              const OptSet& opts = Global::opts());

    Real 
    energyErrgoal() const { return energy_errgoal; }
    void 
    energyErrgoal(Real val) { energy_errgoal = val; }
    
    Real 
    orthWeight() const { return orth_weight; }
    void 
    orthWeight(Real val) { orth_weight = val; }
    
    bool 
    printEigs() const { return printeigs; }
    void 
    printEigs(bool val) { printeigs = val; }
    
    private:

    /////////////
    //
    // Data Members

    Vector center_eigs;
    Real energy_errgoal; //Stop DMRG once energy has converged to this precision
    Real orth_weight;    //How much to penalize non-orthogonality in multiple-state DMRG
    bool printeigs;      //Print slowest decaying eigenvalues after every sweep

    //
    /////////////

    }; // class DMRGObserver

inline DMRGObserver::
DMRGObserver() 
    : 
    energy_errgoal(-1), 
    orth_weight(1),
    printeigs(true)
    { }


void inline DMRGObserver::
measure(int sw, int ha, int b, const SVDWorker& svd, Real energy,
        const OptSet& opts)
    {
    if(printeigs)
        {
        if(b == 1 && ha == 2) 
            {
            Cout << "\n    Largest m during sweep " << sw << " was " << svd.maxEigsKept() << "\n";
            Cout << "    Largest truncation error: " << svd.maxTruncerr() << Endl;
            Vector center_eigs = svd.eigsKept(svd.N()/2);
            Cout << "    Eigs at center bond: ";
            for(int j = 1; j <= min(center_eigs.Length(),10); ++j) 
                {
                Cout << Format(center_eigs(j) > 1E-2 ? ("%.2f") : ("%.2E")) % center_eigs(j);
                Cout << ((j != min(center_eigs.Length(),10)) ? ", " : "");
                }
            Cout << std::endl;
            Cout << Format("    Energy after sweep %d is %f") % sw % energy << Endl;
            }
        }
    }


bool inline DMRGObserver::
checkDone(int sw, const SVDWorker& svd, Real energy,
          const OptSet& opts)
    {
    static Real last_energy;
    
    if(sw == 1) last_energy = 1000;
    if(energy_errgoal > 0 && sw%2 == 0)
        {
        Real dE = fabs(energy-last_energy);
        if(dE < energy_errgoal)
            {
            Cout << Format("    Energy error goal met (dE = %E); returning after %d sweeps.\n") % dE % sw;
            return true;
            }
        }
    last_energy = energy;

    if(fileExists("STOP_DMRG"))
        {
        Cout << "File STOP_DMRG found: stopping this DMRG run after sweep " << sw << Endl;
        system("rm -f STOP_DMRG");
        return true;
        }
    
    return false;
    }

#undef Cout
#undef Endl
#undef Format

#endif // __ITENSOR_DMRGOBSERVER_H
