#include "itensor/all.h"

using namespace itensor;

int main()
    {
    //
    // Single-site wavefunction
    //
    
    auto s = Index("s",2);

    auto psi = ITensor(s); //initialized to zero

    psi.set(s(1),1.);

    // TODO try changing above wavefunction
    //      to be Sx eigenstate

    PrintData(psi);
    
    //
    // Operators 
    //

    auto Sz = ITensor(s,prime(s));
    auto Sx = ITensor(s,prime(s));

    Sz.set(s(1),prime(s)(1),+0.5);
    Sz.set(s(2),prime(s)(2),-0.5);

    Sx.set(s(1),prime(s)(2),+0.5);
    Sx.set(s(2),prime(s)(1),+0.5);

    PrintData(Sz);
    PrintData(Sx);

    //
    // Product Sx * psi 
    //

    //
    // TODO 
    //
    // 1. Compute |phi> = Sx |psi> using
    //    the Sx and psi ITensors above
    //
    // 2. Next compute: auto olap = <psi|phi>;
    //    using the * operator and .real() method
    //
    // 3. Try normalizing |phi> before computing
    //    the overlap: phi /= norm(phi);
    //



    return 0;
    }
