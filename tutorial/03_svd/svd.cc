//
// 2019 Many Electron Collaboration Summer School
// ITensor Tutorial
//
#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

int main()
    {
    //
    // SVD of matrix M
    //

    int Nrow = 4;
    int Ncol = 3;
    auto maxdim = std::min(Nrow,Ncol);

    auto M = Matrix(Nrow,Ncol);
    M(0,0) = 0.435839; M(0,1) = 0.223707; M(0,2) = 0.10;
    M(1,0) = 0.435839; M(1,1) = 0.223707; M(1,2) = -0.10;
    M(2,0) = 0.223707; M(2,1) = 0.435839; M(2,2) = 0.10;
    M(3,0) = 0.223707; M(3,1) = 0.435839; M(3,2) = -0.10;
    Print(M);

    Matrix U,V;
    Vector d;
    SVD(M,U,d,V);

    Print(U);
    Print(d);
    Print(V);

    int nkeep = 2;
    auto Dtrunc = Matrix(maxdim,maxdim);
    for(auto j : range(nkeep))
        {
        Dtrunc(j,j) = d(j);
        }

    auto Mtrunc = U*Dtrunc*transpose(V);
    Print(Mtrunc);

    auto diff = norm(M-Mtrunc);
    auto diff2 = sqr(diff);

    printfln("|M-Mtrunc|^2 = %.2f",diff2);

    println();
    

    //
    // SVD of two-site wavefunction
    //
    
    auto s1 = Index(2,"s1");
    auto s2 = Index(2,"s2");

    auto sing = ITensor(s1,s2);
    auto prod = ITensor(s1,s2);

    //Make sing a singlet
    sing.set(s1=1,s2=2, 1./sqrt(2));
    sing.set(s1=2,s2=1,-1./sqrt(2));

    //Make prod a product state
    prod.set(s1=1,s2=2,1.);

    for(Real mix = 0; mix <= 1.; mix += 0.1)
        {
        //
        // TODO: ADD CODE here to create
        // a new wavefunction that is 
        // (1-mix) times a product state plus (mix)
        // times a singlet (i.e. maximally entangled state).
        //
        // SVD this wavefunction and analyze the results.
        // Try computing and plotting the entanglement entropy.
        //
        }


    return 0;
    }
