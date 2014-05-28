#include "test.h"
#include "iqtensor.h"

using namespace itensor;
using namespace std;

TEST_CASE("Homepage")
    {
    Index a("a",2), b("b",2), c("c",2);
    ITensor Z(a,b), X(b,c);

    commaInit(Z,a,b) << 1, 0, 
                        0, -1;

    commaInit(X,b,c) << 0, 1, 
                        1, 0;

    ITensor R = Z * X;

    CHECK_CLOSE(R(a(1),c(1)),0,1E-10);
    CHECK_CLOSE(R(a(1),c(2)),+1,1E-10);
    CHECK_CLOSE(R(a(2),c(1)),-1,1E-10);
    CHECK_CLOSE(R(a(2),c(2)),0,1E-10);
    }

TEST_CASE("TutorialIntro")
    {
    //cout << "Begin: TutorialIntro -----------------------" << endl;
    Index i("i",10), j("j",20), k("k",30);
    ITensor A(i,j), B(j,k), C;
    A(i(2),j(1)) = 0.5;
    A(i(1),j(2)) = 0.5;

    B(j(1),k(1)) = 1;
    B(j(1),k(1)) = -1;

    C = A * B;
    //PrintDat(C);
    //cout << "End: TutorialIntro -----------------------" << endl;
    }

TEST_CASE("TutorialSimpleMeasurement")
    {
    Index s("s",2,Site);

    ITensor ket(s);

    Real theta = Pi/4;
    ket(s(1)) = cos(theta/2);
    ket(s(2)) = sin(theta/2);


    ITensor Sz(s,prime(s)),
            Sx(s,prime(s));

    commaInit(Sz,s,prime(s)) << 0.5, 0, 
                                 0, -0.5;

    commaInit(Sx,s,prime(s)) << 0, 0.5, 
                                 0.5, 0;

    ITensor bra = dag(prime(ket));

    Real zz = (bra * Sz * ket).toReal();
    Real xx = (bra * Sx * ket).toReal();

    CHECK_CLOSE(zz,0.5*sin(Pi/4),1E-5);
    CHECK_CLOSE(xx,0.5*sin(Pi/4),1E-5);

    //println("<Sz> = ",zz);
    //println("<Sx> = ",xx);
    }
