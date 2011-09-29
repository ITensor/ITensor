#ifndef __ITENSOR_DMRG_WORKER_H
#define __ITENSOR_DMRG_WORKER_H

#include "BaseDMRGWorker.h"
#include "itensor.h"

template <class MPSType, class MPOType, class DefaultDMRGOptsType=DMRGOpts>
class DMRGWorker : public BaseDMRGWorker<MPSType,MPOType,DefaultDMRGOptsType>
{
public:
    typedef BaseDMRGWorker<MPSType,MPOType,DefaultDMRGOptsType> Parent;

    typedef typename MPSType::TensorT Tensor;
    
    DMRGWorker(const Sweeps& sweeps)
    : Parent(sweeps), energy_(0) 
    { }

    DMRGWorker(const Sweeps& sweeps, BaseDMRGOpts& opts)
    : Parent(sweeps, opts), energy_(0) 
    { }

    virtual ~DMRGWorker() { }

    using Parent::sweeps;
    using Parent::opts;

    void run(const MPOType& H, MPSType& psi)
    {
        runInternal(H,psi);
    }

    Real energy() const { return energy_; } 

private:

    Real energy_;
   
    void runInternal(const MPOType& H, MPSType& psi)
    {
        const Real orig_cutoff = psi.cutoff(); 
        const int orig_minm = psi.minm(), 
                  orig_maxm = psi.maxm();
        int debuglevel = (opts().quiet() ? 0 : 1);
        int N = psi.NN();
        energy_ = 0;
        
        psi.position(1);
        //if(H.is_complex()) psi.AAnc(1) *= Complex_1;
        
        std::vector<Tensor> PH(N+1);
        for(int l = N-1; l >= 2; --l) 
            psi.projectOp(l+1,Fromright,PH[l+1],H.AA(l+1),PH[l]);
        
        for(int sw = 1; sw <= sweeps().nsweep(); ++sw)
        {
            psi.cutoff(sweeps().cutoff(sw)); 
            psi.minm(sweeps().minm(sw)); 
            psi.maxm(sweeps().maxm(sw));

            for(int b = 1, ha = 1; ha != 3; sweepnext(b,ha,N))
            {
                if(!opts().quiet()) 
                {
                    std::cout << 
                        boost::format("Sweep=%d, HS=%d, Bond=(%d,%d)\n") 
                        % sw % ha % b % (b+1);
                }
                
                energy_ = psi.bondDavidson(b,H.bondTensor(b),PH[b],PH[b+1],
                                           sweeps().niter(sw),debuglevel,
                                           (ha==1?Fromleft:Fromright));
                
                if(!opts().quiet()) 
                { 
                    std::cout << boost::format("    Truncated to Cutoff=%.1E, Min_m=%d, Max_m=%d\n") 
                    % sweeps().cutoff(sw) % sweeps().minm(sw) % sweeps().maxm(sw);
                    std::cout << boost::format("    Trunc. err=%.1E, States kept=%s\n")
                    % psi.svd().truncerr(b) % psi.LinkInd(b).showm();
                }
                
                opts().measure(sw,ha,b,psi.svd(),energy_);
                
                if(ha == 1 && b != N-1) 
                    psi.projectOp(b,Fromleft,PH[b],H.AA(b),PH[b+1]);
                if(ha == 2 && b != 1)   
                    psi.projectOp(b+1,Fromright,PH[b+1],H.AA(b+1),PH[b]);
            } //for loop over b
            
            if(opts().checkDone(sw,energy_)) break;
        
        } //for loop over sw
        
        psi.cutoff(orig_cutoff); 
        psi.minm(orig_minm); 
        psi.maxm(orig_maxm);
    }

}; // class DMRGWorker

template <class MPSType, class MPOType>
inline Real dmrg(MPSType& psi, const MPOType& H, const Sweeps& sweeps)
{
    DMRGWorker<MPSType,MPOType> worker(sweeps);
    worker.run(H,psi);
    return worker.energy();
}

template <class MPSType, class MPOType>
inline Real dmrg(MPSType& psi, const MPOType& H, const Sweeps& sweeps, 
                 BaseDMRGOpts& opts)
{
    DMRGWorker<MPSType,MPOType> worker(sweeps,opts);
    worker.run(H,psi);
    return worker.energy();
}

#endif // __ITENSOR_DMRG_WORKER_H
