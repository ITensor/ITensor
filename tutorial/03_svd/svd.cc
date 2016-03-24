#include "itensor/util/print_macro.h"
#include "itensor/decomp.h"
#include "itensor/tensor/algs.h"

using namespace itensor;

int
main(int argc, char* argv[])
    {

    //
    // SVD of matrix M
    //

    int Nrow = 4;
    int Ncol = 3;
    auto maxm = std::min(Nrow,Ncol);

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

    auto Dtrunc = Matrix(maxm,maxm);
    stdx::fill(Dtrunc,0);

    int nkeep = 2;
    for(auto j : range1(nkeep))
        Dtrunc(j,j) = d(j);

    auto MM = U*Dtrunc*V;

    auto Diff = MM-M;
    auto D2 = transpose(Diff)*Diff;
    Real n2 = D2(1,1) + D2(2,2) + D2(3,3);

    Print(n2);

    println();
    

    //
    // SVD of two-site wavefunction
    //
    
    auto s1 = Index("s1",2,Site);
    auto s2 = Index("s2",2,Site);

    auto sing = ITensor(s1,s2);
    auto prod = ITensor(s1,s2);

    //Make sing a singlet
    sing.set(s1(1),s2(2), 1./sqrt(2));
    sing.set(s1(2),s2(1),-1./sqrt(2));

    //Make prod a product state
    prod.set(s1(1),s2(2),1.);

    for(Real mix = 0; mix <= 1.; mix += 0.1)
        {
        //
        // Create a new wavefunction that is 
        // (1-mix) times a product state plus (mix)
        // times a singlet (i.e. maximally entangled state).
        //
        // SVD this wavefunction and analyze the results.
        // Try computing and plotting the entanglement entropy.
        //
        }


    return 0;
    }
