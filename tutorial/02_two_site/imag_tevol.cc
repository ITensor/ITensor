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
    // Random initial wavefunction
    //
    auto s1 = Index("s1",2,Site);
    auto s2 = Index("s2",2,Site);
    auto psi = ITensor(s1,s2);
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

    // Initial energy expectation value
    Real initEn = (dag(prime(psi)) * H * psi).real();
    printfln("\nInitial energy = %.10f",initEn);


    //
    // Pieces needed to build exp(-beta*H)
    //
    Real beta = 10;
    int max_order = 100;
    ITensor Id = Id1*Id2;
    ITensor bH = (-beta)*H;

    //Make expH = (bH/M + 1), similar to the innermost
    //parentheses in the formula below:
    ITensor expH = bH/max_order + Id;

    bH.mapprime(1,2);
    bH.mapprime(0,1);

    //
    // The tensor bH now looks like:
    //
    //    s1''   s2''
    //    |      |
    //    ========
    //    |      |
    //    s1'    s2'
    //
    //
    // Exponentiate H to form exp(-beta*H)
    // 
    // Use the formula:
    // exp(bH) = 1 + bH + bH^2/2! + bH^3/3! + ...
    //        = 1 + bH * (1 + bH/2 * (1 + bH/3 * (...)))
    //        ~ ((bH/3 + 1) * bH/2 + 1) * bH + 1
    //
    // to build up exp(bH) in the for loop below.

    for(int ord = max_order-1; ord >= 1; --ord)
        {
        //
        // TODO: ADD CODE here to compute
        // exponential of H, saving result
        // in tensor expH
        //
        //Steps to do:
        //1. Multiply expH times bH
        //2. Divide expH by ord
        //3. Restore primelevels of expH
        //   to their original values
        //   Hint: use expH.mapprime(2,1)
        //4. Add Id to expH
        //
        }

    ITensor psi_beta = expH*psi;
    psi_beta.noprime();
    psi_beta /= norm(psi_beta);

    Real En = (dag(prime(psi_beta)) * H * psi_beta).real();
    printfln("Energy at beta = %.3f: %.10f",beta,En);

    //
    // Inspect entanglement spectrum
    //
    auto A = ITensor(s1);
    ITensor D,B;
    svd(psi_beta,A,D,B);

    PrintData(D);


    return 0;
    }
