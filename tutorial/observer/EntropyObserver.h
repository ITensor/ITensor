//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ENTROPYOBSERVER_H
#define __ITENSOR_ENTROPYOBSERVER_H
#include "observer.h"

namespace itensor {

//
// Class for monitoring DMRG calculations.
// The measure and checkDone methods are virtual
// so that behavior can be customized in a
// derived class.
//

class EntropyObserver : public DMRGObserver<ITensor>
    {
    public:
    
    using Parent = DMRGObserver;
    
    EntropyObserver(const MPS& psi, 
                    const Args& args = Global::args());

    void virtual
    measure(const Args& args = Global::args());

    };

inline EntropyObserver::
EntropyObserver(const MPS& psi, 
                const Args& args) 
  : Parent(psi,args)
    { 
    }

void inline EntropyObserver::
measure(const Args& args)
    {
    auto& psi = Parent::psi();
    auto N = psi.N();
    auto b = args.getInt("AtBond");
    auto sw = args.getInt("Sweep");
    auto ha = args.getInt("HalfSweep");

    auto wf = psi.A(b)*psi.A(b+1);

    auto U = psi.A(b);
    auto V = psi.A(b+1);
    ITensor D;
    auto spectrum = svd(wf,U,D,V);

    //
    // Change the code below to compute
    // the von Neumann or Renyi entropy
    // and print this out
    // (along with the bond number b defined above)
    //
    print("Eigs: ");
    for(auto eig : spectrum.eigsKept())
        {
        print(eig," ");
        }
    println();
    println();
    PAUSE
    }


} //namespace itensor

#endif // __ITENSOR_DMRGOBSERVER_H
