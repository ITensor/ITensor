#include "test.h"

#include "matrix/matrix.h"
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

void
sliceFunc(const vecref& v1r,
           vecref& v2r)
    {
    vecref nr = v1r;
    v2r = nr;
    }

void
sliceFunc(const matrixref& M1r,
          matrixref& M2r)
    {
    matrixref nr = M1r;
    M2r = nr;
    }

TEST_CASE("Check no slicing")
    {
    SECTION("vector")
        {
        //
        // Without defining vecref::operator=(vecref) to
        // be virtual and defining proper overload,
        // sliceFunc1 and 2 could cause v2's store() and data_.data()
        // pointers to get out of sync
        //
        vec v1(10),
            v2(10);
        CHECK(v1.store() == v1.data());
        CHECK(v2.store() == v2.data());

        sliceFunc(v1,v2);
        CHECK(v1.store() == v1.data());
        CHECK(v2.store() == v2.data());
        }

    SECTION("matrix")
        {
        matrix M1(10,10),
               M2(10,10);
        CHECK(M1.store() == M1.data());
        CHECK(M2.store() == M2.data());

        sliceFunc(M1,M2);
        CHECK(M1.store() == M1.data());
        CHECK(M2.store() == M2.data());
        }
    }

TEST_CASE("Test vec and vecref")
{
SECTION("Constructors")
    {
    auto size = 10;
    auto data = randomData(1,size);
    auto vr = vecref(data.begin(),size);
    CHECK(vr);
    CHECK(vr.size() == size);

    for(auto i : count1(size))
        {
        CHECK_CLOSE(vr(i),data(i));
        }

    auto v = vec(size);
    CHECK(v);
    CHECK(v.size() == size);
    const auto* p = v.store();
    for(auto i : count1(size))
        {
        CHECK_CLOSE(v(i),p[i-1]);
        }
    }

SECTION("Construct and assign from vecref")
    {
    auto size = 10;
    auto data = randomData(1,size);
    auto vr = vecref(data.begin(),size);

    //print("data = "); for(auto& el : data) print(el," "); println();
    //print("vr = "); for(auto& el : vr) print(el," "); println();

    vec v1(vr);
    CHECK(v1.size() == size);
    for(auto i : count1(size))
        {
        CHECK_CLOSE(v1(i),data(i));
        }

    vec v2(size+2);
    v2 = vr;
    CHECK(v2.size() == size);
    for(auto i : count1(size))
        {
        CHECK_CLOSE(v2(i),data(i));
        }
    }

SECTION("Scalar multiply, divide")
    {
    auto size = 10;
    auto data = randomData(1,size);
    auto origdata = data;
    auto vr = vecref(data.begin(),size);
    auto fac = Global::random();

    vr *= fac;
    for(auto i : count1(size))
        {
        CHECK_CLOSE(data(i),fac*origdata(i));
        }

    vec v1(vr);
    v1 *= fac;
    for(auto i : count1(size))
        {
        CHECK_CLOSE(v1(i),fac*data(i));
        }
    }

SECTION("Test += -= operators")
    {
    auto size = 10;
    auto dataA = randomData(1,size);
    auto dataB = randomData(1,size);
    auto origdataA = dataA;
    auto A = vecref(dataA.begin(),size);
    auto B = vecref(dataB.begin(),size);

    A += B;
    for(auto i : count1(size))
        {
        CHECK_CLOSE(A(i),B(i)+origdataA(i));
        }

    dataA = origdataA;
    auto cstride = 3;
    auto dataC = randomData(1,cstride*size);
    //Only access every third element:
    auto C = vecref(dataC.begin(),size,cstride);
    A -= C;
    for(auto i : count1(size))
        {
        CHECK_CLOSE(A(i),origdataA(i)-C(i));
        }
    }

SECTION("vec + and -")
    {
    auto size = 20;
    auto v1 = randomVec(size),
         v2 = randomVec(size);
    auto origv1 = v1;

    v1 = v1+v2;
    for(auto i : count1(size))
        {
        CHECK_CLOSE(v1(i),origv1(i)+v2(i));
        }

    v1 = origv1;
    v1 = v1-v2;
    for(auto i : count1(size))
        {
        CHECK_CLOSE(v1(i),origv1(i)-v2(i));
        }

    v1 = origv1;
    auto origv2 = v2;
    v1 = v1+std::move(v2);
    CHECK(!v2);
    for(auto i : count1(size))
        {
        CHECK_CLOSE(v1(i),origv1(i)+origv2(i));
        }
    }
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
        CHECK(A.transposed());

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
    CHECK(At.transposed());

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

    SECTION("Multipy with self")
        {
        auto N = 8;
        auto dataA = randomData(1,N*N);
        auto dataC = randomData(1,N*N);

        auto A = matrixref(dataA.begin(),N,N);
        auto C = matrixref(dataC.begin(),N,N);

        mult(A,A,C);
        for(auto r : count1(C.Nrows()))
        for(auto c : count1(C.Ncols()))
            {
            Real val = 0;
            for(auto k : count1(N)) val += A(r,k)*A(k,c);
            CHECK_CLOSE(C(r,c),val);
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

SECTION("Test += -= operators")
    {
    auto N = 5;
    auto dataA = randomData(1,N*N);
    auto dataB = randomData(1,N*N);
    auto origdataA = dataA;

    auto A = matrixref(dataA.begin(),N,N);
    auto origA = matrixref(origdataA.begin(),N,N);
    auto B = matrixref(dataB.begin(),N,N);
    auto At = A.t();
    auto Bt = B.t();

    A += B;
    for(auto r : count1(N))
    for(auto c : count1(N))
        {
        CHECK_CLOSE(A(r,c),B(r,c)+origA(r,c));
        }

    dataA = origdataA;
    At += Bt;
    for(auto r : count1(N))
    for(auto c : count1(N))
        {
        CHECK_CLOSE(A(r,c),B(r,c)+origA(r,c));
        }

    dataA = origdataA;
    At += B;
    for(auto r : count1(N))
    for(auto c : count1(N))
        {
        CHECK_CLOSE(A(r,c),B(c,r)+origA(r,c));
        }

    dataA = origdataA;
    A -= Bt;
    for(auto r : count1(N))
    for(auto c : count1(N))
        {
        CHECK_CLOSE(A(r,c),-B(c,r)+origA(r,c));
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
        auto A = matrix(Ar,Ac,true);
        for(auto& el : A) el = Global::random();
        CHECK(A.transposed());

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

    //Use operator*
    auto R = A*B;
    for(auto r : count1(C.Nrows()))
    for(auto c : count1(C.Ncols()))
        {
        Real val = 0;
        for(auto k : count1(K)) val += A(r,k)*B(k,c);
        CHECK_CLOSE(R(r,c),val);
        }
    }

SECTION("Assign from ref")
    {
    auto N = 2;
    auto M1 = randomMatrix(N,N);
    auto M2 = randomMatrix(N,N);

    M2 = M1.t();

    auto M1t = M1.t();
    for(auto r : count1(N))
    for(auto c : count1(N))
        CHECK_CLOSE(M2(r,c),M1t(r,c));
    }

SECTION("Addition / Subtraction")
    {
    auto Nr = 4,
         Nc = 5;
    auto A = randomMatrix(Nr,Nc);
    auto B = randomMatrix(Nr,Nc);

    auto C = A;
    C += B;
    for(auto r : count1(Nr))
    for(auto c : count1(Nc))
        CHECK_CLOSE(C(r,c),A(r,c)+B(r,c));

    C = A;
    C -= B;
    for(auto r : count1(Nr))
    for(auto c : count1(Nc))
        CHECK_CLOSE(C(r,c),A(r,c)-B(r,c));

    auto D = B;
    C = A + std::move(D);
    CHECK(!D);
    CHECK(B);
    for(auto r : count1(Nr))
    for(auto c : count1(Nc))
        CHECK_CLOSE(C(r,c),A(r,c)+B(r,c));

    D = B;
    CHECK(D);
    C = A - std::move(D);
    CHECK(!D);
    CHECK(B);
    for(auto r : count1(Nr))
    for(auto c : count1(Nc))
        CHECK_CLOSE(C(r,c),A(r,c)-B(r,c));
    }

SECTION("Scalar multiply, divide")
    {
    auto N = 10;
    auto A = randomMatrix(N,N);
    auto origA = A;
    auto fac = Global::random();

    A *= fac;
    for(auto r : count1(N))
    for(auto c : count1(N))
        CHECK_CLOSE(A(r,c),origA(r,c)*fac);

    A = origA;
    A /= fac;
    for(auto r : count1(N))
    for(auto c : count1(N))
        CHECK_CLOSE(A(r,c),origA(r,c)/fac);

    A = origA;
    auto At = A.t();
    At *= fac;
    for(auto r : count1(N))
    for(auto c : count1(N))
        CHECK_CLOSE(A(r,c),origA(r,c)*fac);

    A = origA;
    auto S = subMatrix(A,1,N/2,1,N/2);
    S *= fac;
    for(auto r : count1(N/2))
    for(auto c : count1(N/2))
        CHECK_CLOSE(A(r,c),origA(r,c)*fac);
    }

} // Test matrix


TEST_CASE("Test slicing")
{


SECTION("Diagonal")
    {
    auto A = randomMatrix(4,4);
    auto d = diagonal(A);
    long i = 1;
    for(const auto& el : d)
        {
        CHECK_CLOSE(el,A(i,i));
        ++i;
        }

    //Set diag els of A to 100 through d
    for(auto& el : d) el = 100;

    for(auto j : count1(d.size()))
        {
        CHECK_CLOSE(100,A(j,j));
        }

    A = randomMatrix(5,10);
    d = diagonal(A);
    i = 1;
    for(const auto& el : d)
        {
        CHECK_CLOSE(el,A(i,i));
        ++i;
        }

    A = randomMatrix(10,5);
    d = diagonal(A);
    i = 1;
    for(const auto& el : d)
        {
        CHECK_CLOSE(el,A(i,i));
        ++i;
        }
    for(auto j : count1(d.size()))
        {
        CHECK_CLOSE(d(j),A(j,j));
        }
    }

SECTION("Row / Col Slicing")
    {
    auto N = 10;
    auto A = randomMatrix(N,N);

    for(auto r : count1(N))
        {
        auto R = row(A,r);
        for(auto c : count1(N))
            CHECK_CLOSE(A(r,c),R(c));
        }

    for(auto c : count1(N))
        {
        auto C = column(A,c);
        for(auto r : count1(N))
            CHECK_CLOSE(A(r,c),C(r));
        }
    }

SECTION("Transpose")
    {
    auto nr = 10,
         nc = 15;
    auto A = randomMatrix(nr,nc);
    auto At = matrixref(A.cstore(),transpose(A.ind()));
    CHECK(At.transposed());

    for(auto i : count1(nr))
    for(auto j : count1(nc))
        {
        CHECK_CLOSE(A(i,j),At.get(j,i));
        }
    }

SECTION("Sub Vector")
    {
    auto v = vec(20);
    for(auto& el : v) el = Global::random();

    long start = 1,
         stop = 10;
    auto s = subVector(v,start,stop);
    long count = 0;
    for(auto j : count1(start,stop))
        {
        CHECK_CLOSE(s(j),v(start-1+j));
        ++count;
        }
    CHECK(count == s.size());

    start = 3;
    stop = 12;
    s = subVector(v,start,stop);
    count = 0;
    for(auto j : count1(start,stop))
        {
        CHECK_CLOSE(s(j),v(start-1+j));
        ++count;
        }
    CHECK(count == s.size());

    //Set elements of v through s
    for(auto& el : s) el = -1;
    for(auto j : count1(start,stop))
        CHECK_CLOSE(-1,v(j));
    }

SECTION("Sub Matrix")
    {
    auto nr = 10,
         nc = 15;
    auto A = randomMatrix(nr,nc);

    auto rstart = 1,
         rstop = 3,
         cstart = 1,
         cstop = 4;
    auto S = subMatrix(A,rstart,rstop,cstart,cstop);
    for(auto r : count1(S.Nrows()))
    for(auto c : count1(S.Ncols()))
        {
        CHECK_CLOSE(S(r,c),A(rstart-1+r,cstart-1+c));
        }

    rstart = 2;
    cstart = 3;
    auto At = A.t();
    CHECK(At.transposed());
    S = subMatrix(At,rstart,rstop,cstart,cstop);
    for(auto r : count1(S.Nrows()))
    for(auto c : count1(S.Ncols()))
        {
        CHECK_CLOSE(S(r,c),At(rstart-1+r,cstart-1+c));
        }

    rstart = 3;
    rstop = 8;
    cstart = 5;
    cstop = 10;
    S = subMatrix(A,rstart,rstop,cstart,cstop);
    for(auto r : count1(S.Nrows()))
    for(auto c : count1(S.Ncols()))
        {
        CHECK_CLOSE(S(r,c),A(rstart-1+r,cstart-1+c));
        }

    }
} //Test slicing

TEST_CASE("Matrix Algorithms and Decompositions")
{

SECTION("diagSymmetric")
    {
    auto N = 100;
    auto M = randomMatrix(N,N);
    //Symmetrize:
    M = M+M.t();

    matrix U;
    vec d;
    diagSymmetric(M,U,d);

    auto D = matrix(N,N);
    diagonal(D) = d;
    auto R = U*D*U.t();

    for(auto r : count1(N))
    for(auto c : count1(N))
        {
        CHECK_CLOSE(R(r,c),M(r,c));
        }
    CHECK(norm(R-M) < 1E-12*norm(M));
    }

}
