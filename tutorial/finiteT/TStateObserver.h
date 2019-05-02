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

class TStateObserver : public TEvolObserver
    {
    public:

    TStateObserver(const MPS& psi,
                   const Args& args = Args::global());

    virtual ~TStateObserver() { }

    void virtual
    measure(const Args& args = Args::global());
    
    private:

    /////////////
    //
    // Data Members

    const MPS& psi_;
    bool show_maxdim_;

    //
    /////////////

    }; // class TStateObserver

inline TStateObserver::
TStateObserver(const MPS& psi,
               const Args& args) 
    :
    psi_(psi)
    { 
    show_maxdim_ = args.getBool("ShowMaxDim",true);
    }


void inline TStateObserver::
measure(const Args& args)
    {
    const auto t = args.getReal("Time");
    if(show_maxdim_)
        {
        const auto ttotal = args.getReal("TotalTime");
        const Real percentdone = (100.*t)/ttotal;
        long maxdim = 0;
        for(int b = 1; b < psi_.length(); ++b)
            {
            maxdim = std::max(maxdim,dim(linkIndex(psi_,b)));
            }
        printfln("%2.f%%:%d ",percentdone,maxdim);
        }
    }

}


#endif // __ITENSOR_TSTATEOBSERVER_H
