#include "core.h"

using namespace std;
using boost::format;
using namespace itensor;

ITensor
makeSp(const Index& s)
    {
    ITensor Sp(s,primed(s));
    Sp(s(2),primed(s)(1)) = 1;
    return Sp;
    }

ITensor
makeSm(const Index& s)
    {
    ITensor Sm(s,primed(s));
    Sm(s(1),primed(s)(2)) = 1;
    return Sm;
    }

ITensor
makeSz(const Index& s)
    {
    ITensor Sz(s,primed(s));
    Sz(s(1),primed(s)(1)) =  0.5;
    Sz(s(2),primed(s)(2)) = -0.5;
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

    psi(s1(1),s2(2)) =  1./sqrt(2);
    psi(s1(2),s2(1)) = -1./sqrt(2);

    PrintData(psi);

    //exit(0);

    //
    // Single-site operators
    //

    ITensor Sz1 = makeSz(s1),
            Sz2 = makeSz(s2),
            Sp1 = makeSp(s1),
            Sp2 = makeSp(s2),
            Sm1 = makeSm(s1),
            Sm2 = makeSm(s2);

    PrintData(Sz1);
    PrintData(Sp1);
    PrintData(Sm1);

    //exit(0);

    //
    // Two-site Heisenberg Hamiltonian
    //

    ITensor H = Sz1*Sz2 + 0.5*Sp1*Sm2 + 0.5*Sm1*Sp2;


    //
    // Energy expectation value
    //

    ITensor cpsi = conj(primed(psi));
    Real E = (cpsi * H * psi).toReal();

    Print(E);

    return 0;
    }
