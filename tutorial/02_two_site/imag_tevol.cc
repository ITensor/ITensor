#include "itensor/all.h"

using namespace itensor;

ITensor
makeId(Index const& s)
    {
    auto Id = ITensor(s,prime(s));
    Id.set(s(1),prime(s)(1),1);
    Id.set(s(2),prime(s)(2),1);
    return Id;
    }

ITensor
makeSp(Index const& s)
    {
    auto Sp = ITensor(s,prime(s));
    Sp.set(s(2),prime(s)(1), 1);
    return Sp;
    }

ITensor
makeSm(Index const& s)
    {
    auto Sm = ITensor(s,prime(s));
    Sm.set(s(1),prime(s)(2),1);
    return Sm;
    }

ITensor
makeSz(Index const& s)
    {
    auto Sz = ITensor(s,prime(s));
    Sz.set(s(1),prime(s)(1), 0.5);
    Sz.set(s(2),prime(s)(2),-0.5);
    return Sz;
    }


int main()
    {

    //
    // Two-site wavefunction
    // initialized to a singlet
    //
    
    auto s1 = Index("s1",2,Site);
    auto s2 = Index("s2",2,Site);

    auto psi = ITensor(s1,s2); //default initialized to zero

    randomize(psi);
    psi *= 1./norm(psi);

    //PrintData(psi);

    //
    // Single-site operators
    //

    auto Sz1 = makeSz(s1);
    auto Sz2 = makeSz(s2);
    auto Sp1 = makeSp(s1);
    auto Sp2 = makeSp(s2);
    auto Sm1 = makeSm(s1);
    auto Sm2 = makeSm(s2);
    auto Id1 = makeId(s1);
    auto Id2 = makeId(s2);

    //
    // Two-site Heisenberg Hamiltonian
    //

    ITensor H = Sz1*Sz2 + 0.5*Sp1*Sm2 + 0.5*Sm1*Sp2;


    //
    // Energy expectation value
    //

    auto cpsi = dag(prime(psi));
    Real initEn = (cpsi * H * psi).real();

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
        // TODO: ADD CODE here to compute
        // exponential of H, saving result
        // in tensor expH
        //
        }

    ITensor psi_beta = expH*psi;
    psi_beta.noprime();
    psi_beta *= 1./norm(psi_beta);

    Real En = (dag(prime(psi_beta)) * H * psi_beta).real();
    printfln("Energy at beta = %.3f: %.10f",beta,En);

    ITensor A(s1),B(s2),D;

    svd(psi_beta,A,D,B);

    PrintData(D);


    return 0;
    }
