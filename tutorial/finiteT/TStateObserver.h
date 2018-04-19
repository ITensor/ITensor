//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TSTATEOBSERVER_H
#define __ITENSOR_TSTATEOBSERVER_H
#include "itensor/mps/TEvolObserver.h"

//
// Class for monitoring time evolution calculations.
// The measure and checkDone methods are virtual
// so that behavior can be customized in a
// derived class.
//

namespace itensor {

template<class Tensor>
class TStateObserver : public TEvolObserver
    {
    public:

    using MPST = MPSt<Tensor>;
    
    TStateObserver(const MPST& psi,
                   const Args& args = Global::args());

    virtual ~TStateObserver() { }

    void virtual
    measure(const Args& args = Global::args());
    
    private:

    /////////////
    //
    // Data Members

    const MPST& psi_;
    bool show_maxm_;

    //
    /////////////

    }; // class TStateObserver

template<class Tensor>
inline TStateObserver<Tensor>::
TStateObserver(const MPST& psi,
               const Args& args) 
    :
    psi_(psi)
    { 
    show_maxm_ = args.getBool("ShowMaxm",true);
    }


template<class Tensor>
void inline TStateObserver<Tensor>::
measure(const Args& args)
    {
    const auto t = args.getReal("Time");
    if(show_maxm_)
        {
        const auto ttotal = args.getReal("TotalTime");
        const Real percentdone = (100.*t)/ttotal;
        long maxm = 0;
        for(int b = 1; b < psi_.N(); ++b)
            {
            maxm = std::max(maxm,linkInd(psi_,b).m());
            }
        printfln("%2.f%%:%d ",percentdone,maxm);
        }
    }

}


#endif // __ITENSOR_TSTATEOBSERVER_H
