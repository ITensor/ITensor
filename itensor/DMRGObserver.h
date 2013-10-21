//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_DMRGOBSERVER_H
#define __ITENSOR_DMRGOBSERVER_H
#include "observer.h"

#define Cout std::cout
#define Endl std::endl
#define Format boost::format

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
    measure(int N, int sw, int ha, int b, const Spectrum& spec, Real energy,
            const OptSet& opts = Global::opts());
    
    bool virtual
    checkDone(int sw, Real energy,
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

    const MPSt<Tensor>& 
    psi() const { return psi_; }
    
    private:

    /////////////
    //
    // Data Members

    const MPSt<Tensor>& psi_;

    Vector center_eigs;
    Real energy_errgoal; //Stop DMRG once energy has converged to this precision
    Real orth_weight;    //How much to penalize non-orthogonality in multiple-state DMRG
    bool printeigs;      //Print slowest decaying eigenvalues after every sweep
    int max_eigs;
    Real max_te;
    bool done_;

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
    orth_weight(opts.getReal("OrthWeight",1)),
    printeigs(opts.getBool("PrintEigs",true)),
    max_eigs(-1),
    max_te(-1),
    done_(false),
    default_ops_(psi.model().defaultOps())
    { 
    }


template<class Tensor>
void inline DMRGObserver<Tensor>::
measure(int N, int sw, int ha, int b, const Spectrum& spec, Real energy,
        const OptSet& opts)
    {
    if(!opts.getBool("Quiet",false) && !opts.getBool("NoMeasure",false))
        {
        for(size_t j = 0; j < default_ops_.size(); ++j)
            {
            const std::string opname = default_ops_.at(j);
            Complex z = 
                BraKet(primed(psi_.A(b),Site),psi_.model().op(opname,b)*psi_.A(b));
            if(fabs(z.imag()) < 1E-14)
                Cout << Format("<%s>(%d) = %.10E") % opname % b % z.real() << Endl;
            else
                Cout << Format("<%s>(%d) = (%.10E,%.10E)") % opname % b % z.real() % z.imag() << Endl;
            }
        }

    if(printeigs)
        {
        if(b == N/2 && ha == 2)
            {
            Cout << Endl;
            Vector center_eigs = spec.eigsKept();
            Real S = 0;
            for(int j = 1; j <= center_eigs.Length(); ++j) 
                {
                S -= center_eigs(j)*log(fabs(center_eigs(j)));
                }
            Cout << Format("    vN Entropy at center bond b=%d = %.10f") % (N/2) % S << Endl;
            Cout << Format("    Eigs at center bond b=%d: ") % (N/2);
            for(int j = 1; j <= min(center_eigs.Length(),10); ++j) 
                {
                const Real eig = center_eigs(j);
                if(eig < 1E-3) break;
                Cout << Format("%.4f ") % eig;
                }
            Cout << Endl;
            }
        }

    max_eigs = max(max_eigs,spec.numEigsKept());
    max_te = max(max_te,spec.truncerr());
    if(b == 1 && ha == 2) 
        {
        if(!printeigs) Cout << Endl;
        if(max_eigs > 0)
            {
            Cout << "    Largest m during sweep " << sw << " was " << max_eigs << "\n";
            max_eigs = -1;
            }
        if(max_te > 0)
            {
            Cout << "    Largest truncation error: " << max_te << Endl;
            max_te = -1;
            }
        Cout << Format("    Energy after sweep %d is %f") % sw % energy << Endl;
        }

    }


template<class Tensor>
bool inline DMRGObserver<Tensor>::
checkDone(int sw, Real energy,
          const OptSet& opts)
    {
    static Real last_energy;
    
    if(sw == 1) last_energy = 1000;
    if(energy_errgoal > 0 && sw%2 == 0)
        {
        Real dE = fabs(energy-last_energy);
        if(dE < energy_errgoal)
            {
            Cout << Format("    Energy error goal met (dE = %.3E < %.3E); returning after %d sweeps.") 
                    % dE
                    % energy_errgoal
                    % sw
                 << Endl;
            done_ = true;
            return done_;
            }
        }
    last_energy = energy;

    //If STOP_DMRG found, will return true (i.e. done) once, but 
    //outer calling using same Observer may continue running e.g. infinite dmrg calling finite dmrg.
    if(fileExists("STOP_DMRG"))
        {
        Cout << "File STOP_DMRG found: stopping this DMRG run after sweep " << sw << Endl;
        system("rm -f STOP_DMRG");
        return true;
        }

    //Set done_ flag to true so any outer callers using this Observer will also terminate.
    if(fileExists("STOP_DMRG_ALL"))
        {
        Cout << "File STOP_DMRG_ALL found: stopping this run after sweep " << sw << Endl;
        system("rm -f STOP_DMRG_ALL");
        done_ = true;
        return done_;
        }
    
    return done_;
    }

#undef Cout
#undef Endl
#undef Format

#endif // __ITENSOR_DMRGOBSERVER_H
