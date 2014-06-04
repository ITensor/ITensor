#include "core.h"

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
    psi(s(1)) = 1;

    PrintData(psi);
    
    //exit(0); //uncomment to exit here

    //
    // Operators 
    //

    ITensor Sz(s,primed(s)),
            Sx(s,primed(s));

    commaInit(Sz,s,primed(s)) = 0.5, 0.0,
                                0.0,-0.5;

    commaInit(Sx,s,primed(s)) = 0.0, 0.5,
                                0.5, 0.0;

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

    const Real theta = Pi/4;

    //Extra factors of two come from S=1/2 representation
    psi(s(1)) = cos(theta/2.);
    psi(s(2)) = sin(theta/2.);

    PrintData(psi);

    //exit(0); //uncomment to exit here

    //
    // Expectation values
    //

    ITensor cpsi = dag(primed(psi));

    Real zz = (cpsi * Sz * psi).toReal();
    Real xx = (cpsi * Sx * psi).toReal();

    println("<Sz> = ", zz);
    println("<Sx> = ", xx);

    return 0;
    }
