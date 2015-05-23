#include "itensor/core.h"

using namespace itensor;

int
main(int argc, char* argv[])
    {
    //
    // Single-site wavefunction
    //
    
    //Make a dimension 2 Index
    Index s("s",2);

    //Construct an ITensor
    ITensor psi(s); //default initialized to zero


    //
    // Initialize up spin
    //

    //Set first element to 1.
    psi.set(1,s(1));

    PrintData(psi);
    
    //exit(0); //uncomment to exit here

    //
    // Operators 
    //

    ITensor Sz(s,prime(s)),
            Sx(s,prime(s));
    Sz.set(+0.5,s(1),prime(s)(1));
    Sz.set(-0.5,s(1),prime(s)(1));

    Sx.set(+0.5,s(1),prime(s)(2));
    Sx.set(+0.5,s(2),prime(s)(1));

    PrintData(Sz);
    PrintData(Sx);

    //exit(0); //uncomment to exit here

    //
    // Product Sx * phi 
    //

    ITensor phi = Sx * psi;

    phi.noprime();

    PrintData(phi);

    //exit(0); //uncomment to exit here

    //
    // 45* angle spin
    //

    Real theta = Pi/4;

    //Extra factors of two come from S=1/2 representation
    psi.set(cos(theta/2),s(1));
    psi.set(sin(theta/2),s(2));

    PrintData(psi);

    //exit(0); //uncomment to exit here

    //
    // Expectation values
    //

    ITensor cpsi = dag(prime(psi));

    Real zz = (cpsi * Sz * psi).real();
    Real xx = (cpsi * Sx * psi).real();

    println("<Sz> = ", zz);
    println("<Sx> = ", xx);

    return 0;
    }
