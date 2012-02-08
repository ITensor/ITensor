#ifndef __ITENSOR_DMRGOPTS_H
#define __ITENSOR_DMRGOPTS_H
#include "BaseDMRGOpts.h"

//
// Class for fine-tuning DMRG algorithms.
// The measure and checkDone methods are virtual
// so that behavior can be customized in a
// derived class.
//

class DMRGOpts : public BaseDMRGOpts
    {
    public:

    bool 
    quiet() const { return quiet_; }
    void 
    quiet(bool val) { quiet_ = val; }
    
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
    
    DMRGOpts();

    virtual ~DMRGOpts() { }
    
    virtual void 
    measure(int sw, int ha, int b, const SVDWorker& svd, Real energy);
    
    virtual bool 
    checkDone(int sw, Real energy);

    private:

    Vector center_eigs;
    Real energy_errgoal; //Stop DMRG once energy has converged to this precision
    Real orth_weight;    //How much to penalize non-orthogonality in multiple-state DMRG
    bool printeigs;      //Print slowest decaying eigenvalues after every sweep
    bool quiet_;         //Show/don't show info after every step

    }; // class DMRGOpts

inline DMRGOpts::
DMRGOpts() 
    : energy_errgoal(-1), 
      orth_weight(1),
      printeigs(true), 
      quiet_(true)
    { }

inline
void DMRGOpts::
measure(int sw, int ha, int b, const SVDWorker& svd, Real energy)
    {
    if(printeigs)
        {
        if(b == 1 && ha == 2) 
            {
            std::cout << "\n    Largest m during sweep " << sw << " was " << svd.maxEigsKept() << "\n";
            std::cout << "    Largest truncation error: " << svd.maxTruncerr() << "\n";
            Vector center_eigs = svd.eigsKept(svd.NN()/2);
            std::cout << "    Eigs at center bond: ";
            for(int j = 1; j <= min(center_eigs.Length(),10); ++j) 
                {
                std::cout << boost::format(center_eigs(j) > 1E-2 ? ("%.2f") : ("%.2E")) % center_eigs(j);
                std::cout << ((j != min(center_eigs.Length(),10)) ? ", " : "\n");
                }
            std::cout << boost::format("    Energy after sweep %d is %f\n") % sw % energy;
            }
        }
    }

inline
bool DMRGOpts::
checkDone(int sw, Real energy)
    {
    static Real last_energy;
    
    if(sw == 1) last_energy = 1000;
    if(energy_errgoal > 0 && sw%2 == 0)
        {
        Real dE = fabs(energy-last_energy);
        if(dE < energy_errgoal)
            {
            std::cout << boost::format("    Energy error goal met (dE = %E); returning after %d sweeps.\n") % dE % sw;
            return true;
            }
        }
    last_energy = energy;
    
    return false;
    }

#endif // __ITENSOR_DMRGOPTS_H
