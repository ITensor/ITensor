#ifndef DMRG_WORKER_H
#define DMRG_WORKER_H

#include "BaseDMRGOpts.h" //<-- abstract base class 
#include "DMRGOpts.h" // <-- default implementation
#include "itensor.h"
#include "Sweeps.h"

template <class Tensor=ITensor,class MPSType=MPS, class MPOType=MPO>
class DMRGWorker {
public:
    typedef DMRGOpts DMRGOptsType; // choose what type to use here
    
    DMRGWorker(const Sweeps& sweeps)
    : sweeps_(sweeps),opts_(),energy_(0) { }
    
    void run(MPSType& psi, const MPOType& H)
    {
        runInternal(psi,H);
    }

    const Real& energy() const { return energy_; } 
    
private:
   
    void runInternal(MPSType& psi, const MPOType& H)
    {
        const Real orig_cutoff = psi.cutoff(); 
        const int orig_minm = psi.minm(), orig_maxm = psi.maxm();
        int debuglevel = (opts_.quiet ? 0 : 1);
        int N = psi.NN();
        energy_ = 0;
        
        psi.position(1);
        //if(H.is_complex()) psi.AAnc(1) *= Complex_1;
        
        std::vector<Tensor> PH(N+1);
        for(int l = N-1; l >= 2; --l) psi.projectOp(l+1,Fromright,PH[l+1],H.AA(l+1),PH[l]);
        
        for(int sw = 1; sw <= sweeps_.nsweep(); ++sw)
        {
            psi.cutoff(sweeps_.cutoff(sw)); psi.minm(sweeps_.minm(sw)); psi.maxm(sweeps_.maxm(sw));
            for(int b = 1, ha = 1; ha != 3; sweepnext(b,ha,N))
            {
                if(!opts_.quiet) std::cout << boost::format("Sweep=%d, HS=%d, Bond=(%d,%d)\n") % sw % ha % b % (b+1);
                
                energy_ = psi.bondDavidson(b,H.bondTensor(b),PH[b],PH[b+1],
                                          sweeps_.niter(sw),debuglevel,(ha==1?Fromleft:Fromright));
                
                if(!opts_.quiet) 
                { 
                    std::cout << boost::format("    Truncated to Cutoff=%.1E, Min_m=%d, Max_m=%d\n") 
                    % sweeps_.cutoff(sw) % sweeps_.minm(sw) % sweeps_.maxm(sw);
                    std::cout << boost::format("    Trunc. err=%.1E, States kept=%s\n")
                    % psi.svd().truncerr(b) % psi.LinkInd(b).showm();
                }
                
                opts_.measure(sw,ha,b,psi.svd(),energy_);
                
                if(ha == 1 && b != N-1) psi.projectOp(b,Fromleft,PH[b],H.AA(b),PH[b+1]);
                if(ha == 2 && b != 1)   psi.projectOp(b+1,Fromright,PH[b+1],H.AA(b+1),PH[b]);
            } //for loop over b
            
        if(opts_.checkDone(sw,energy_))
        {
            psi.cutoff(orig_cutoff); 
            psi.minm(orig_minm); 
            psi.maxm(orig_maxm);
            return;
        }
        
        } //for loop over sw
        
        psi.cutoff(orig_cutoff); 
        psi.minm(orig_minm); 
        psi.maxm(orig_maxm);
        }

    Sweeps sweeps_;
    Real energy_;
    DMRGOptsType opts_;
}; // class DMRGWorker

#endif // DMRG_WORKER_H
