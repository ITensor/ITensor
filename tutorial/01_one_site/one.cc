//
// 2019 Many Electron Collaboration Summer School
// ITensor Tutorial
//
#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

int main()
    {
    //
    // Define our Index 
    // (the space we are working with)
    //

    auto s = Index(2,"s");

    //
    // Operators 
    //

    auto Sx = ITensor(s,prime(s));

    Sx.set(s=1,prime(s)=2,+0.5);
    Sx.set(s=2,prime(s)=1,+0.5);

    PrintData(Sx);

    //
    // Single-site wavefunction
    //
    
    auto psi = ITensor(s); //initialized to zero

    //
    // TODO 
    //
    // 1. make the above wavefunction
    //    the (normalized) positive Sx eigenstate
    //    HINT: use psi.set(...)
    //

    /* Your code here */

    PrintData(psi);
    
    //
    // TODO
    //
    // 2. Compute |phi> = Sx |psi> using
    //    the Sx and psi ITensors above
    //    AND
    //    compute: auto olap = <psi|phi>
    //    using the * operator and elt(...) method.
    //    Print the result with PrintData(...).
    //

    /* Your code here */

    //
    // TODO
    //
    // 3. Try normalizing |phi> and recompute
    //    the inner product <psi|phi>
    //    Print the result with PrintData(...).
    //    HINT: use phi /= norm(phi)) to normalize.
    //

    /* Your code here */

    return 0;
    }
