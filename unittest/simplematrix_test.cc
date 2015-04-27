#include "test.h"

#include "autovector.h"
#include "global.h"
#include "count.h"
#include "matrix/algs.h"

using namespace itensor;
using namespace std;

autovector<Real>
randomData(long first, long last)
    {
    autovector<Real> data(first,last);
    for(auto& el : data) el = Global::random();
    return data;
    }

TEST_CASE("Test VecRef")
{
SECTION("Constructors")
    {
    auto size = 10;
    auto data = randomData(1,size);
    auto vr = VecRef(data.begin(),size);
    CHECK(vr);
    CHECK(vr.size() == size);

    for(auto i : count1(size))
        {
        CHECK_CLOSE(vr(i),data(i));
        }

    auto v = Vec(size);
    CHECK(v);
    CHECK(v.size() == size);
    const auto* p = v.data();
    for(auto i : count1(size))
        {
        CHECK_CLOSE(v(i),p[i-1]);
        }
    }

SECTION("Automatic Conversion")
    {
    auto size = 10;
    auto data1 = randomData(1,size);
    auto data2 = randomData(1,size);
    auto vr1 = VecRef(data1.begin(),size);
    auto cvr2 = CVecRef(data2.begin(),size);

    //Assigning to CVecRef from VecRef is ok
    cvr2 = vr1;
    CHECK(cvr2.data() == vr1.data());

    //Constructing CVecRef from VecRef is ok
    CVecRef cref(vr1);
    CHECK(cref.data() == vr1.data());

    //This is not allowed - no conversion
    //from CVecRef back to VecRef
    //vr1 = cvr2;
    }

SECTION("Construct and assign from VecRef")
    {
    auto size = 10;
    auto data = randomData(1,size);
    auto vr = VecRef(data.begin(),size);

    //print("data = "); for(auto& el : data) print(el," "); println();
    //print("vr = "); for(auto& el : vr) print(el," "); println();

    Vec v1(vr);
    CHECK(v1.size() == size);
    for(auto i : count1(size))
        {
        CHECK_CLOSE(v1(i),data(i));
        }

    Vec v2(size+2);
    v2 = vr;
    CHECK(v2.size() == size);
    for(auto i : count1(size))
        {
        CHECK_CLOSE(v2(i),data(i));
        }
    }

SECTION("Copy data")
    {
    auto size = 10;
    auto data1 = randomData(1,size);
    auto data2 = randomData(1,size);
    auto vr1 = VecRef(data1.begin(),size);
    auto vr2 = VecRef(data2.begin(),size);
    auto v1 = randomVec(size);
    auto v2 = randomVec(4*size);

    vr1 &= v1;
    for(auto i : count1(size))
        CHECK_CLOSE(data1(i),v1(i));

    vr2 &= makeVecRef(v2,size,4);
    for(auto i : count1(size))
        CHECK_CLOSE(data2(i),v2(1+4*(i-1)));
    }

SECTION("Scalar multiply, divide")
    {
    auto size = 10;
    auto data = randomData(1,size);
    auto origdata = data;
    auto vr = VecRef(data.begin(),size);
    auto fac = Global::random();

    vr *= fac;
    for(auto i : count1(size))
        {
        CHECK_CLOSE(data(i),fac*origdata(i));
        }

    Vec v1(vr);
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
    auto A = VecRef(dataA.begin(),size);
    auto B = VecRef(dataB.begin(),size);

    A += B;
    for(auto i : count1(size))
        {
        CHECK_CLOSE(A(i),B(i)+origdataA(i));
        }

    dataA = origdataA;
    auto cstride = 3;
    auto dataC = randomData(1,cstride*size);
    //Only access every third element:
    auto C = VecRef(dataC.begin(),size,cstride);
    A -= C;
    for(auto i : count1(size))
        {
        CHECK_CLOSE(A(i),origdataA(i)-C(i));
        }
    }

SECTION("Vec + and -")
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

SECTION("Dot product")
    {
    auto N = 10;

    auto v = randomVec(N);
    auto vnrm = norm(v);
    CHECK_CLOSE(v*v,vnrm*vnrm);

    ////Non-trivial stride case:
    //auto M1 = randomMat(N,N);
    //auto M2 = randomMat(N,N);
    //auto d1 = diagonal(M1);
    //auto d2 = diagonal(M2);

    //auto dot = d1*d2;
    //Real val = 0;
    //for(auto j : count1(N))
    //    {
    //    val += M1(j,j)*M2(j,j);
    //    }
    //CHECK_CLOSE(val,dot);
    }
}

TEST_CASE("Test MatRef")
{

SECTION("Constructors")
    {
    SECTION("Regular Constructor")
        {
        auto Ar = 3,
             Ac = 4;
        auto data = randomData(1,Ar*Ac);
        CHECK(data.size() == Ar*Ac);
        auto A = makeMatRef(data.begin(),Ar,Ac);

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
        auto A = makeMatRef(data.begin(),Ar,Ac,true);
        CHECK(A.transposed());

        for(auto r : count1(Ar))
        for(auto c : count1(Ac))
            {
            CHECK_CLOSE(A(c,r),data[r+Ar*(c-1)]);
            }
        }
    }

SECTION("Automatic Conversion")
    {
    auto N = 10;
    auto data1 = randomData(1,N*N);
    auto data2 = randomData(1,N*N);
    auto mr1 = MatRef(data1.begin(),N,N);
    auto cmr2 = CMatRef(data2.begin(),N,N);

    //Assigning to CMatRef from MatRef is ok
    cmr2 = mr1;
    CHECK(cmr2.data() == mr1.data());

    //Constructing CVecRef from VecRef is ok
    CMatRef cref(mr1);
    CHECK(cref.data() == mr1.data());

    //This is not allowed - no conversion
    //from CMatRef back to MatRef
    //mr1 = cmr2;
    }

SECTION("Transpose Method")
    {
    auto Ar = 5,
         Ac = 3;
    auto data = randomData(1,Ar*Ac);
    auto A = makeMatRef(data.begin(),Ar,Ac);

    auto At = A.t();
    CHECK(At.transposed());

    for(auto r : count1(Ar))
    for(auto c : count1(Ac))
        {
        CHECK_CLOSE(A(r,c),At(c,r));
        }
    }

SECTION("Test += -= operators")
    {
    auto N = 5;
    auto dataA = randomData(1,N*N);
    auto dataB = randomData(1,N*N);
    auto origdataA = dataA;

    auto A = makeMatRef(dataA.begin(),N,N);
    auto origA = makeMatRef(origdataA.cbegin(),N,N);
    auto B = makeMatRef(dataB.cbegin(),N,N);
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

//SECTION("Test addition of refs to same data")
//    {
//    auto N = 4;
//    auto M = randomMat(N,N);
//    auto origM = M;
//    matrixref& Mr1 = M;
//    const matrixref& Mr2 = M;
//
//    Mr1 += Mr2;
//    for(auto r : count1(N)) 
//    for(auto c : count1(N)) 
//        CHECK_CLOSE(M(r,c),2*origM(r,c));
//
//    Mr1 -= Mr2;
//    for(auto& el : M) CHECK(el < 1E-10);
//    }

SECTION("Test MatRef mult")
    {
    SECTION("Case 1")
        {
        auto Ar = 3,
             K  = 4,
             Bc = 5;
        auto dataA = randomData(1,Ar*K);
        auto dataB = randomData(1,K*Bc);
        auto dataC = autovector<Real>(1,Ar*Bc);

        auto A = makeMatRef(dataA.begin(),Ar,K);
        auto B = makeMatRef(dataB.begin(),K,Bc);
        auto C = makeMatRef(dataC.begin(),Ar,Bc);

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

        auto A = makeMatRef(dataA.begin(),K,Ac);
        auto B = makeMatRef(dataB.begin(),K,Bc);
        auto C = makeMatRef(dataC.begin(),Ac,Bc);

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

        auto A = makeMatRef(dataA.begin(),Ar,K);
        auto B = makeMatRef(dataB.begin(),Br,K);
        auto C = makeMatRef(dataC.begin(),Ar,Br);

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

        auto A = makeMatRef(dataA.begin(),K,Ac);
        auto B = makeMatRef(dataB.begin(),Br,K);
        auto C = makeMatRef(dataC.begin(),Ac,Br);

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

        auto A = makeMatRef(dataA.begin(),Ar,K);
        auto B = makeMatRef(dataB.begin(),K,Bc);
        auto C = makeMatRef(dataC.begin(),Bc,Ar);

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

        auto A = makeMatRef(dataA.begin(),K,Ac);
        auto B = makeMatRef(dataB.begin(),K,Bc);
        auto C = makeMatRef(dataC.begin(),Bc,Ac);

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

        auto A = makeMatRef(dataA.begin(),Ar,K);
        auto B = makeMatRef(dataB.begin(),Br,K);
        auto C = makeMatRef(dataC.begin(),Br,Ar);

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

        auto A = makeMatRef(dataA.begin(),K,Ac);
        auto B = makeMatRef(dataB.begin(),Br,K);
        auto C = makeMatRef(dataC.begin(),Br,Ac);

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

        auto A = makeMatRef(dataA.begin(),N,N);
        auto C = makeMatRef(dataC.begin(),N,N);

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

SECTION("Test multAdd")
    {
    auto Ar = 3,
         K  = 4,
         Bc = 5;
    auto dataA = randomData(1,Ar*K);
    auto dataB = randomData(1,K*Bc);
    auto dataC = autovector<Real>(1,Ar*Bc);

    auto A = makeMatRef(dataA.begin(),Ar,K);
    auto B = makeMatRef(dataB.begin(),K,Bc);
    auto C = makeMatRef(dataC.begin(),Ar,Bc);

    //Save a copy of C's original data in order
    //to explicitly carry out multAdd alg. below
    auto orig_dataC = dataC;
    auto origC = makeMatRef(orig_dataC.begin(),Ar,Bc);

    multAdd(A,B,C);
    for(auto r : count1(C.Nrows()))
    for(auto c : count1(C.Ncols()))
        {
        Real val = 0;
        for(auto k : count1(K)) val += A(r,k)*B(k,c) + origC(r,c);
        CHECK_CLOSE(C(r,c),val);
        }
    }


} //Test MatRef


TEST_CASE("Test Mat")
{
SECTION("Constructors")
    {
    SECTION("Constructor Basics")
        {
        auto Ar = 3,
             Ac = 4;
        auto A = Mat(Ar,Ac);
        CHECK(A.Nrows() == Ar);
        CHECK(A.Ncols() == Ac);
        CHECK(A);

        Mat B;
        CHECK(!B);
        }

    SECTION("Regular Constructor")
        {
        auto Ar = 3,
             Ac = 4;
        auto A = randomMat(Ar,Ac);

        const auto *data = A.data();
        for(auto r : count(Ar))
        for(auto c : count(Ac))
            {
            CHECK_CLOSE(A(1+r,1+c),data[r+Ar*c]);
            }
        }
    }

SECTION("Assign from ref")
    {
    auto N = 4;
    auto M1 = randomMat(N,N);
    auto M2 = randomMat(N,N);

    auto M1t = transpose(M1);
    M2 = M1t;

    for(auto r : count1(N))
    for(auto c : count1(N))
        CHECK_CLOSE(M2(r,c),M1t(r,c));

    auto M1tt = transpose(transpose(M1));
    M2 = M1tt;
    for(auto r : count1(N))
    for(auto c : count1(N))
        CHECK_CLOSE(M2(r,c),M1tt(r,c));
    }

SECTION("Test Mat multiplication")
    {
    auto Ar = 3,
         K  = 4,
         Bc = 5;

    //Multiply matrices A*B, store in matrix C
    auto A = Mat(Ar,K);
    auto B = Mat(K,Bc);
    auto C = Mat(Ar,Bc);
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
    auto Cref = MatRef(dataC.begin(),Ar,Bc);
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


SECTION("Addition / Subtraction")
    {
    auto Nr = 4,
         Nc = 5;
    auto A = randomMat(Nr,Nc);
    auto B = randomMat(Nr,Nc);

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
    auto A = randomMat(N,N);
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
    auto At = transpose(A);
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

SECTION("Matrix-vector product")
    {
    SECTION("Square case")
        {
        auto N = 10;
        auto M = randomMat(N,N);
        auto x = randomVec(N);

        auto y = M*x;
        for(auto r : count1(N))
            {
            Real val = 0;
            for(auto c : count1(N)) val += M(r,c)*x(c);
            CHECK_CLOSE(y(r),val);
            }

        y = transpose(M)*x;
        for(auto r : count1(N))
            {
            Real val = 0;
            for(auto c : count1(N)) val += M(c,r)*x(c);
            CHECK_CLOSE(y(r),val);
            }

        y = x*M;
        for(auto c : count1(N))
            {
            Real val = 0;
            for(auto r : count1(N)) val += x(r)*M(r,c);
            CHECK_CLOSE(y(c),val);
            }

        y = x*transpose(M);
        for(auto c : count1(N))
            {
            Real val = 0;
            for(auto r : count1(N)) val += x(r)*M(c,r);
            CHECK_CLOSE(y(c),val);
            }

        //Check a case where vector is strided
        auto d = diagonal(M);
        y = M*d;
        for(auto r : count1(N))
            {
            Real val = 0;
            for(auto c : count1(N)) val += M(r,c)*M(c,c);
            CHECK_CLOSE(y(r),val);
            }
        } //Square case

    SECTION("Rectangular case")
        {
        auto Ar = 5,
             Ac = 10;
        auto A = randomMat(Ar,Ac);
        auto vR = randomVec(Ac);
        auto vL = randomVec(Ar);

        auto res = A*vR;
        for(auto r : count1(Ar))
            {
            Real val = 0;
            for(auto c : count1(Ac)) val += A(r,c)*vR(c);
            CHECK_CLOSE(res(r),val);
            }

        res = vL*A;
        for(auto c : count1(Ac))
            {
            Real val = 0;
            for(auto r : count1(Ar)) val += vL(r)*A(r,c);
            CHECK_CLOSE(res(c),val);
            }

        res = transpose(A)*vL;
        for(auto c : count1(Ac))
            {
            Real val = 0;
            for(auto r : count1(Ar)) val += vL(r)*A(r,c);
            CHECK_CLOSE(res(c),val);
            }

        res = vR*transpose(A);
        for(auto r : count1(Ar))
            {
            Real val = 0;
            for(auto c : count1(Ac)) val += A(r,c)*vR(c);
            CHECK_CLOSE(res(r),val);
            }

        } //Rectangular case
    }

} // Test matrix


TEST_CASE("Test slicing")
{


SECTION("Diagonal")
    {
    auto A = randomMat(4,4);
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

    A = randomMat(5,10);
    d = diagonal(A);
    i = 1;
    for(const auto& el : d)
        {
        CHECK_CLOSE(el,A(i,i));
        ++i;
        }

    A = randomMat(10,5);
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
    auto A = randomMat(N,N);

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

    auto S = transpose(subMatrix(A,1,N/2,1,N/2));
    for(auto c : count1(N/2))
        {
        auto C = column(S,c);
        for(auto r : count1(N/2))
            CHECK_CLOSE(S(r,c),C(r));
        }
    }

SECTION("Transpose")
    {
    auto nr = 10,
         nc = 15;
    auto A = randomMat(nr,nc);
    auto At = transpose(MatRef(A.data(),A.ind()));
    CHECK(At.transposed());

    for(auto i : count1(nr))
    for(auto j : count1(nc))
        {
        CHECK_CLOSE(A(i,j),At(j,i));
        }
    }

SECTION("Sub Vector")
    {
    auto v = Vec(20);
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
    auto A = randomMat(nr,nc);

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
    auto At = transpose(A);
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
    auto N = 10;
    auto M = randomMat(N,N);
    //Symmetrize:
    M = M+transpose(M);

    Mat U;
    Vec d;
    diagSymmetric(M,U,d);

    auto D = Mat(N,N);
    diagonal(D) &= d;
    auto R = U*D*transpose(U);

    for(auto r : count1(N))
    for(auto c : count1(N))
        {
        CHECK_CLOSE(R(r,c),M(r,c));
        }
    CHECK(norm(R-M) < 1E-12*norm(M));
    }

SECTION("Orthogonalize")
    {
    auto N = 10;
    auto M = randomMat(N,N);

    orthog(M);

    auto R = transpose(M)*M;
    for(auto r : count1(N))
    for(auto c : count1(N))
        {
        if(r == c) CHECK_CLOSE(R(r,c),1);
        else       CHECK(R(r,c) < 1E-12);
        }
    }

SECTION("Singular Value Decomp")
    {
    //auto N = 10;
    //auto M = randomMat(N,N);

    //Mat U,V;
    //Vec d;
    //SVD(M,U,d,V);

    }
}
