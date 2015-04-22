#include "test.h"

#include "matrix/simplematrix.h"
#include "autovector.h"
#include "global.h"
#include "count.h"

using namespace itensor;
using namespace std;

autovector<Real>
randomData(long first, long last)
    {
    autovector<Real> data(first,last);
    for(auto& el : data) el = Global::random();
    return data;
    }

matrix
randomMatrix(long nr, long nc, bool trans = false)
    {
    matrix M(nr,nc,trans);
    for(auto& el : M) el = Global::random();
    return M;
    }

TEST_CASE("Test matrixref")
{

SECTION("Constructors")
    {
    SECTION("Regular Constructor")
        {
        auto Ar = 3,
             Ac = 4;
        auto data = randomData(1,Ar*Ac);
        CHECK(data.size() == Ar*Ac);
        auto A = matrixref(data.begin(),Ar,Ac);

        //println("data = ");
        //for(const auto& el : data) print(el,", ");
        //println();
        //Print(A);
        //exit(0);
        for(auto r : count1(Ar))
        for(auto c : count1(Ac))
            {
            CHECK_CLOSE(A(r,c),data[r+Ar*(c-1)]);
            }
        }
    SECTION("Transpose Constructor")
        {
        auto Ar = 3,
             Ac = 4;
        auto data = randomData(1,Ar*Ac);
        CHECK(data.size() == Ar*Ac);
        auto A = matrixref(data.begin(),Ar,Ac,true);

        for(auto r : count1(Ar))
        for(auto c : count1(Ac))
            {
            CHECK_CLOSE(A(c,r),data[r+Ar*(c-1)]);
            }
        }
    }

SECTION("Transpose Method")
    {
    auto Ar = 5,
         Ac = 3;
    auto data = randomData(1,Ar*Ac);
    auto A = matrixref(data.begin(),Ar,Ac);

    auto At = A.t();

    for(auto r : count1(Ar))
    for(auto c : count1(Ac))
        {
        CHECK_CLOSE(A(r,c),At(c,r));
        }
    }

SECTION("Read Only")
    {
    auto Ar = 5,
         Ac = 3;
    auto data = randomData(1,Ar*Ac);

    auto A = matrixref(data.begin(),Ar,Ac);
    CHECK(!A.readOnly());
    A(1,2) = 7;
    CHECK(A(1,2) == 7);

    const auto& cdata = data;
    auto cA = matrixref(cdata.begin(),Ar,Ac);
    CHECK(cA.readOnly());
    }

SECTION("Test matrixref mult")
    {
    SECTION("Case 1")
        {
        auto Ar = 3,
             K  = 4,
             Bc = 5;
        auto dataA = randomData(1,Ar*K);
        auto dataB = randomData(1,K*Bc);
        auto dataC = autovector<Real>(1,Ar*Bc);

        auto A = matrixref(dataA.begin(),Ar,K);
        auto B = matrixref(dataB.begin(),K,Bc);
        auto C = matrixref(dataC.begin(),Ar,Bc);

        mult(A,B,C);
        for(auto r : count1(C.Nrows()))
        for(auto c : count1(C.Ncols()))
            {
            Real val = 0;
            for(auto k : count1(K)) val += A(r,k)*B(k,c);
            CHECK_CLOSE(C(r,c),val);
            }
        }

    SECTION("Case 2")
        {
        auto Ac = 3,
             K  = 4,
             Bc = 5;
        auto dataA = randomData(1,K*Ac);
        auto dataB = randomData(1,K*Bc);
        auto dataC = autovector<Real>(1,Ac*Bc);

        auto A = matrixref(dataA.begin(),K,Ac);
        auto B = matrixref(dataB.begin(),K,Bc);
        auto C = matrixref(dataC.begin(),Ac,Bc);

        auto At = A.t();
        mult(At,B,C);
        for(auto r : count1(C.Nrows()))
        for(auto c : count1(C.Ncols()))
            {
            Real val = 0;
            for(auto k : count1(K)) val += At(r,k)*B(k,c);
            CHECK_CLOSE(C(r,c),val);
            }
        }

    SECTION("Case 3")
        {
        auto Ar = 3,
             K =  4,
             Br = 5;
        auto dataA = randomData(1,Ar*K);
        auto dataB = randomData(1,Br*K);
        auto dataC = autovector<Real>(1,Ar*Br);

        auto A = matrixref(dataA.begin(),Ar,K);
        auto B = matrixref(dataB.begin(),Br,K);
        auto C = matrixref(dataC.begin(),Ar,Br);

        auto Bt = B.t();
        mult(A,Bt,C);
        for(auto r : count1(C.Nrows()))
        for(auto c : count1(C.Ncols()))
            {
            Real val = 0;
            for(auto k : count1(K)) val += A(r,k)*Bt(k,c);
            CHECK_CLOSE(C(r,c),val);
            }
        }

    SECTION("Case 4")
        {
        auto Ac = 3,
             K =  4,
             Br = 5;
        auto dataA = randomData(1,K*Ac);
        auto dataB = randomData(1,Br*K);
        auto dataC = autovector<Real>(1,Ac*Br);

        auto A = matrixref(dataA.begin(),K,Ac);
        auto B = matrixref(dataB.begin(),Br,K);
        auto C = matrixref(dataC.begin(),Ac,Br);

        auto At = A.t();
        auto Bt = B.t();
        mult(At,Bt,C);
        for(auto r : count1(C.Nrows()))
        for(auto c : count1(C.Ncols()))
            {
            Real val = 0;
            for(auto k : count1(K)) val += At(r,k)*Bt(k,c);
            CHECK_CLOSE(C(r,c),val);
            }
        }

    SECTION("Case 5")
        {
        auto Ar = 3,
             K =  4,
             Bc = 5;
        auto dataA = randomData(1,Ar*K);
        auto dataB = randomData(1,K*Bc);
        auto dataC = autovector<Real>(1,Ar*Bc);

        auto A = matrixref(dataA.begin(),Ar,K);
        auto B = matrixref(dataB.begin(),K,Bc);
        auto C = matrixref(dataC.begin(),Bc,Ar);

        auto Ct = C.t();
        mult(A,B,Ct);
        for(auto r : count1(Ct.Nrows()))
        for(auto c : count1(Ct.Ncols()))
            {
            Real val = 0;
            for(auto k : count1(K)) val += A(r,k)*B(k,c);
            CHECK_CLOSE(Ct(r,c),val);
            }
        }

    SECTION("Case 6")
        {
        auto Ac = 3,
             K =  4,
             Bc = 5;
        auto dataA = randomData(1,K*Ac);
        auto dataB = randomData(1,K*Bc);
        auto dataC = autovector<Real>(1,Bc*Ac);

        auto A = matrixref(dataA.begin(),K,Ac);
        auto B = matrixref(dataB.begin(),K,Bc);
        auto C = matrixref(dataC.begin(),Bc,Ac);

        auto At = A.t();
        auto Ct = C.t();
        mult(At,B,Ct);
        for(auto r : count1(Ct.Nrows()))
        for(auto c : count1(Ct.Ncols()))
            {
            Real val = 0;
            for(auto k : count1(K)) val += At(r,k)*B(k,c);
            CHECK_CLOSE(Ct(r,c),val);
            }
        }

    SECTION("Case 7")
        {
        auto Ar = 3,
             K =  4,
             Br = 5;
        auto dataA = randomData(1,Ar*K);
        auto dataB = randomData(1,Br*K);
        auto dataC = autovector<Real>(1,Br*Ar);

        auto A = matrixref(dataA.begin(),Ar,K);
        auto B = matrixref(dataB.begin(),Br,K);
        auto C = matrixref(dataC.begin(),Br,Ar);

        auto Bt = B.t();
        auto Ct = C.t();
        mult(A,Bt,Ct);
        for(auto r : count1(Ct.Nrows()))
        for(auto c : count1(Ct.Ncols()))
            {
            Real val = 0;
            for(auto k : count1(K)) val += A(r,k)*Bt(k,c);
            CHECK_CLOSE(Ct(r,c),val);
            }
        }

    SECTION("Case 8")
        {
        auto Ac = 3,
             K =  4,
             Br = 5;
        auto dataA = randomData(1,K*Ac);
        auto dataB = randomData(1,Br*K);
        auto dataC = autovector<Real>(1,Br*Ac);

        auto A = matrixref(dataA.begin(),K,Ac);
        auto B = matrixref(dataB.begin(),Br,K);
        auto C = matrixref(dataC.begin(),Br,Ac);

        auto At = A.t();
        auto Bt = B.t();
        auto Ct = C.t();
        mult(At,Bt,Ct);
        for(auto r : count1(Ct.Nrows()))
        for(auto c : count1(Ct.Ncols()))
            {
            Real val = 0;
            for(auto k : count1(K)) val += At(r,k)*Bt(k,c);
            CHECK_CLOSE(Ct(r,c),val);
            }
        }
    }

SECTION("Test mult_add")
    {
    auto Ar = 3,
         K  = 4,
         Bc = 5;
    auto dataA = randomData(1,Ar*K);
    auto dataB = randomData(1,K*Bc);
    auto dataC = autovector<Real>(1,Ar*Bc);

    auto A = matrixref(dataA.begin(),Ar,K);
    auto B = matrixref(dataB.begin(),K,Bc);
    auto C = matrixref(dataC.begin(),Ar,Bc);

    //Save a copy of C's original data in order
    //to explicitly carry out mult_add alg. below
    auto orig_dataC = dataC;
    auto origC = matrixref(orig_dataC.begin(),Ar,Bc);

    mult_add(A,B,C);
    for(auto r : count1(C.Nrows()))
    for(auto c : count1(C.Ncols()))
        {
        Real val = 0;
        for(auto k : count1(K)) val += A(r,k)*B(k,c) + origC(r,c);
        CHECK_CLOSE(C(r,c),val);
        }
    }
}

TEST_CASE("Test matrix")
{
SECTION("Constructors")
    {
    SECTION("Constructor Basics")
        {
        auto Ar = 3,
             Ac = 4;
        auto A = matrix(Ar,Ac);
        CHECK(A.Nrows() == Ar);
        CHECK(A.Ncols() == Ac);
        CHECK(A);

        matrix B;
        CHECK(!B);
        }

    SECTION("Regular Constructor")
        {
        auto Ar = 3,
             Ac = 4;
        auto A = randomMatrix(Ar,Ac);

        const auto *data = A.store();
        for(auto r : count(Ar))
        for(auto c : count(Ac))
            {
            CHECK_CLOSE(A(1+r,1+c),data[r+Ar*c]);
            }
        }
    SECTION("Transpose Constructor")
        {
        auto Ar = 3,
             Ac = 4;
        auto A = randomMatrix(Ar,Ac,true);

        const auto *data = A.store();
        for(auto r : count(Ar))
        for(auto c : count(Ac))
            {
            CHECK_CLOSE(A(1+c,1+r),data[r+Ar*c]);
            }
        }
    }

SECTION("Test matrix mult")
    {
    auto Ar = 3,
         K  = 4,
         Bc = 5;

    //Multiply matrices A*B, store in matrix C
    auto A = matrix(Ar,K);
    auto B = matrix(K,Bc);
    auto C = matrix(Ar,Bc);
    mult(A,B,C);
    for(auto r : count1(C.Nrows()))
    for(auto c : count1(C.Ncols()))
        {
        Real val = 0;
        for(auto k : count1(K)) val += A(r,k)*B(k,c);
        CHECK_CLOSE(C(r,c),val);
        }


    //Store result in a matrixref instead
    auto dataC = autovector<Real>(1,Ar*Bc);
    auto Cref = matrixref(dataC.begin(),Ar,Bc);
    mult(A,B,Cref);
    for(auto r : count1(C.Nrows()))
    for(auto c : count1(C.Ncols()))
        {
        Real val = 0;
        for(auto k : count1(K)) val += A(r,k)*B(k,c);
        CHECK_CLOSE(Cref(r,c),val);
        }
    }

} // Test matrix
