//
// 2019 Many Electron Collaboration Summer School
// ITensor Tutorial
//
#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

ITensor
makeSp(Index const& s)
    {
    auto Sp = ITensor(s,prime(s));
    Sp.set(s=2,prime(s)=1, 1);
    return Sp;
    }

ITensor
makeSm(Index const& s)
    {
    auto Sm = ITensor(s,prime(s));
    Sm.set(s=1,prime(s)=2,1);
    return Sm;
    }

ITensor
makeSz(Index const& s)
    {
    auto Sz = ITensor(s,prime(s));
    Sz.set(s=1,prime(s)=1,+0.5);
    Sz.set(s=2,prime(s)=2,-0.5);
    return Sz;
    }


int main()
    {
    //
    // Initial product state
    //
    auto s1 = Index(2,"s1");
    auto s2 = Index(2,"s2");
    auto psi = ITensor(s1,s2);
    psi.set(s1=1,s2=2,1.0);

    PrintData(psi);

    //
    // Single-site operators
    //
    auto Sz1 = makeSz(s1);
    auto Sz2 = makeSz(s2);
    auto Sp1 = makeSp(s1);
    auto Sp2 = makeSp(s2);
    auto Sm1 = makeSm(s1);
    auto Sm2 = makeSm(s2);

    //
    // Two-site Heisenberg Hamiltonian
    //
    auto H = Sz1*Sz2 + 0.5*Sp1*Sm2 + 0.5*Sm1*Sp2;

    // Initial energy expectation value
    auto initEn = elt(dag(prime(psi)) * H * psi);
    printfln("\nInitial energy = %.10f",initEn);


    //
    // Make exp(-beta*H)
    //
    // TODO
    //
    // 3. Adjust beta to get the ground state
    //
    Real beta = 0.;
    auto expH = expHermitian(H,-beta);

    // Here we apply exp(-beta*H), normalize
    // and unprime
    auto psibeta = expH*psi;
    psibeta.noPrime();
    psibeta /= norm(psibeta);

    PrintData(psibeta);

    auto En = elt(dag(prime(psibeta)) * H * psibeta);
    printfln("At beta=%.1f, energy = %.10f",beta,En);

    //
    // TODO
    //
    // 1. Adjust the following code to
    //    truncate to dimension 1.
    //    HINT: use the ITensor named argument
    //    system, e.g. {"MaxDim=",...}
    //

    auto [U,D,V] = svd(psibeta,{s1} /* Your code here */ );

    PrintData(D);

    //
    // TODO
    //
    // 2. Calculate the overlap of the new
    //    wavefunction with the old wavefunction.
    //    Print your results with PrintData(...).
    //    HINT: use U*D*V to calculate the new,
    //    truncated wavefunction
    //    

    /* Your code here */

    //
    // TODO
    //
    // 3. Increase beta (defined above) to get the
    //    ground state. How does the overlap 
    //    change?
    //    

    return 0;
    }
