#include "itensor/all.h"

using namespace itensor;

int main()
    {
    //
    // Single-site wavefunction
    //
    
    //Make a dimension 2 Index
    auto s = Index("s",2);

    //Construct an ITensor
    auto psi = ITensor(s); //default initialized to zero

    //
    // Initialize up spin
    //

    //Set first element to 1.
    psi.set(s(1),1);

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
    // Product Sx * phi 
    //

    ITensor phi = Sx * psi;

    phi.noprime();

    PrintData(phi);

    //
    // 45* angle spin
    //

    Real theta = Pi/4;

    //Extra factors of two come from S=1/2 representation
    psi.set(s(1),cos(theta/2));
    psi.set(s(2),sin(theta/2));

    PrintData(psi);

    //
    // Expectation values
    //

    auto cpsi = dag(prime(psi));

    Real zz = (cpsi * Sz * psi).real();
    Real xx = (cpsi * Sx * psi).real();

    println("<Sz> = ", zz);
    println("<Sx> = ", xx);

    return 0;
    }
