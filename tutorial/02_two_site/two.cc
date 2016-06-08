#include "itensor/all.h"

using namespace itensor;

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

    psi.set(s1(1),s2(2), 1./sqrt(2));
    psi.set(s1(2),s2(1),-1./sqrt(2));

    PrintData(psi);

    //EXIT //uncomment to exit here

    //
    // Single-site operators
    //

    auto Sz1 = makeSz(s1);
    auto Sz2 = makeSz(s2);
    auto Sp1 = makeSp(s1);
    auto Sp2 = makeSp(s2);
    auto Sm1 = makeSm(s1);
    auto Sm2 = makeSm(s2);

    PrintData(Sz1);
    PrintData(Sp1);
    PrintData(Sm1);

    //EXIT //uncomment to exit here

    //
    // Two-site Heisenberg Hamiltonian
    //

    auto H = Sz1*Sz2 + 0.5*Sp1*Sm2 + 0.5*Sm1*Sp2;

    //
    // Energy expectation value
    //

    auto cpsi = dag(prime(psi));
    Real E = (cpsi * H * psi).real();

    Print(E);

    return 0;
    }
