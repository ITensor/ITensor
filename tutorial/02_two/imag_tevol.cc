#include "core.h"

using namespace std;
using boost::format;

ITensor
makeId(const Index& s)
    {
    ITensor Id(s,primed(s));
    Id(s(1),primed(s)(1)) = 1;
    Id(s(2),primed(s)(2)) = 1;
    return Id;
    }

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

    psi.randomize();

    cout << "psi.norm() = " << psi.norm() << endl;
    
    psi *= 1./psi.norm();

    cout << "psi.norm() = " << psi.norm() << endl;

    PrintDat(psi);

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

    ITensor cpsi = conj(primed(psi));
    Real initEn = (cpsi * H * psi).toReal();

    cout << format("Initial energy = %.10f")
            % initEn
            << endl;

    //
    // Exponentiate H to form exp(-beta*H/2)
    // 
    // Use formula:
    // exp(x) = 1 + x + x^2/2! + x^3/3! + ...
    //        = 1 + x * (1 + x/2 * (1 + x/3 * (...)))
    //        ~ ((x/3 + 1) * x/2 + 1) * x + 1

    const Real beta = 3;
    const int max_order = 100;

    const ITensor I = Id1*Id2;

    ITensor x = H*(-beta/2.);

    ITensor expH = x;

    x.mapprime(1,2);
    x.mapprime(0,1);

    for(int ord = max_order; ord >= 1; --ord)
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

    Real En = (conj(primed(psi_beta)) * H * psi_beta).toReal();
    cout << format("Energy at beta = %.3f: %.10f")
            % beta
            % En
            << endl;

    ITensor A(s1),B(s2),D;

    svd(psi_beta,A,D,B);

    PrintDat(D);


    return 0;
    }
