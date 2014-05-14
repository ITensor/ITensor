#include "core.h"

using namespace std;
using boost::format;
using namespace itensor;

int
main(int argc, char* argv[])
    {

    //
    // SVD of matrix M
    //

    const int Nrow = 4;
    const int Ncol = 3;
    const int maxm = min(Nrow,Ncol);

    Matrix M(Nrow,Ncol);
    M(1,1) = 0.435839; M(1,2) = 0.223707; M(1,3) = 0.10;
    M(2,1) = 0.435839; M(2,2) = 0.223707; M(2,3) = -0.10;
    M(3,1) = 0.223707; M(3,2) = 0.435839; M(3,3) = 0.10;
    M(4,1) = 0.223707; M(4,2) = 0.435839; M(4,3) = -0.10;
    Print(M);

    Matrix U,V;
    Vector d;
    SVD(M,U,d,V);

    Print(U);
    Print(d);
    Print(V);

    Matrix Dtrunc(maxm,maxm);
    Dtrunc = 0;

    const int nkeep = 2;
    for(int j = 1; j <= nkeep; ++j)
        Dtrunc(j,j) = d(j);

    Matrix MM = U*Dtrunc*V;

    Matrix Diff = MM-M;
    Matrix D2 = Diff.t()*Diff;
    Real n2 = D2(1,1) + D2(2,2) + D2(3,3);

    Print(n2);
    

    //
    // SVD of two-site wavefunction
    //
    
    Index s1("s1",2,Site),
          s2("s2",2,Site);

    ITensor sing(s1,s2),
            prod(s1,s2);

    //Make sing a singlet
    sing(s1(1),s2(2)) =  1./sqrt(2);
    sing(s1(2),s2(1)) = -1./sqrt(2);

    //Make prod a product state
    prod(s1(1),s2(2)) =  1.;

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
