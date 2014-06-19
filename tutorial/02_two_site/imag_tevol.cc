#include "core.h"

using namespace itensor;

ITensor
makeId(const Index& s)
    {
    ITensor Id(s,prime(s));
    Id(s(1),prime(s)(1)) = 1;
    Id(s(2),prime(s)(2)) = 1;
    return Id;
    }

ITensor
makeSp(const Index& s)
    {
    ITensor Sp(s,prime(s));
    Sp(s(2),prime(s)(1)) = 1;
    return Sp;
    }

ITensor
makeSm(const Index& s)
    {
    ITensor Sm(s,prime(s));
    Sm(s(1),prime(s)(2)) = 1;
    return Sm;
    }

ITensor
makeSz(const Index& s)
    {
    ITensor Sz(s,prime(s));
    Sz(s(1),prime(s)(1)) =  0.5;
    Sz(s(2),prime(s)(2)) = -0.5;
    return Sz;
    }

int
main(int argc, char* argv[])
    {

    //
    // Two-site wavefunction
    // initialized to a singlet
    //
    
    Index s1("s1",2,Site),
          s2("s2",2,Site);

    ITensor psi(s1,s2); //default initialized to zero

    psi.randomize();
    psi *= 1./psi.norm();

    //PrintData(psi);

    //
    // Single-site operators
    //

    ITensor Sz1 = makeSz(s1),
            Sz2 = makeSz(s2),
            Sp1 = makeSp(s1),
            Sp2 = makeSp(s2),
            Sm1 = makeSm(s1),
            Sm2 = makeSm(s2),
            Id1 = makeId(s1),
            Id2 = makeId(s2);

    //
    // Two-site Heisenberg Hamiltonian
    //

    ITensor H = Sz1*Sz2 + 0.5*Sp1*Sm2 + 0.5*Sm1*Sp2;


    //
    // Energy expectation value
    //

    ITensor cpsi = dag(prime(psi));
    Real initEn = (cpsi * H * psi).toReal();

    printfln("\nInitial energy = %.10f",initEn);

    //
    // Exponentiate H to form exp(-beta*H/2)
    // 
    // Use the formula:
    // exp(x) = 1 + x + x^2/2! + x^3/3! + ...
    //        = 1 + x * (1 + x/2 * (1 + x/3 * (...)))
    //        ~ ((x/3 + 1) * x/2 + 1) * x + 1
    //
    // to build up exp(x) in the for loop below.

    const Real beta = 10;
    const int max_order = 100;

    const ITensor I = Id1*Id2;

    ITensor x = H*(-beta);

    //Make expH = (x/M + 1), similar to the innermost
    //parentheses in the formula above:
    ITensor expH = x/max_order + I;

    x.mapprime(1,2);
    x.mapprime(0,1);

    //
    // The tensor x now looks like:
    //
    //    s1''   s2''
    //    |      |
    //    ========
    //    |      |
    //    s1'    s2'
    //

    for(int ord = max_order-1; ord >= 1; --ord)
        {
        //
        //
        // Fill in code here.
        //
        //
        }

    ITensor psi_beta = expH*psi;
    psi_beta.noprime();
    psi_beta *= 1./psi_beta.norm();

    Real En = (dag(prime(psi_beta)) * H * psi_beta).toReal();
    printfln("Energy at beta = %.3f: %.10f",beta,En);

    ITensor A(s1),B(s2),D;

    svd(psi_beta,A,D,B);

    PrintData(D);


    return 0;
    }
