#include "test.h"

#include "itensor/util/autovector.h"
#include "itensor/util/count.h"
#include "itensor/tensor/algs.h"
#include "itensor/global.h"

using namespace itensor;
using namespace std;

autovector<Real>
randomData(long first, long last)
    {
    autovector<Real> data(first,last);
    for(auto& el : data) el = Global::random();
    return data;
    }

autovector<Real>
randomData(long size)
    {
    autovector<Real> data(0,size-1);
    for(auto& el : data) el = Global::random();
    return data;
    }

TEST_CASE("Test VectorRef and Vector")
{

SECTION("Test makeRef")
    {
    auto size = 10;
    auto data = randomData(size);

    auto vr = makeVecRef(data.begin(),size);
    auto cvr = makeVecRefc(data.begin(),size);
    auto v = Vector(size);
    const Vector& cv = v;

    auto ref1 = makeRef(vr);
    CHECK(ref1.data() == data.begin());

    auto ref2 = makeRef(cvr);
    CHECK(ref2.data() == data.begin());
    CHECK((std::is_same<decltype(ref2),VectorRefc>::value));

    auto ref3 = makeRef(v);
    CHECK(ref3.data() == v.data());
    CHECK((std::is_same<decltype(ref3),VectorRef>::value));

    auto ref4 = makeRef(cv);
    CHECK(ref4.data() == cv.data());
    CHECK((std::is_same<decltype(ref4),VectorRefc>::value));

    //Not allowed to make ref of temporary:
    //auto ref5 = makeRef(Vector(size)); //not allowed
    //auto ref6 = makeRef(std::move(v)); //not allowed
    }


SECTION("Constructors")
    {
    auto size = 10;
    auto data = randomData(size);
    auto vr = makeVecRef(data.begin(),size);
    CHECK(vr);
    CHECK(vr.size() == size);
    CHECK(isContiguous(vr));

    for(auto i : count(size))
        {
        CHECK_CLOSE(vr(i),data(i));
        }

    auto v = Vector(size);
    CHECK(v);
    CHECK(v.size() == size);
    const auto* p = v.data();
    for(auto i : count(size))
        {
        CHECK_CLOSE(v(i),p[i]);
        }
    CHECK(isContiguous(v));
    }

SECTION("VectorRef Conversion")
    {
    auto size = 10;
    auto data1 = randomData(size);
    auto data2 = randomData(size);
    auto vr1 = makeVecRef(data1.begin(),size);
    auto cvr2 = makeVecRef(data2.begin(),size);

    //Assigning to VectorRefc from VectorRef is ok
    cvr2 = vr1;
    CHECK(cvr2.data() == vr1.data());

    //Constructing VectorRefc from VectorRef is ok
    VectorRefc cref(vr1);
    CHECK(cref.data() == vr1.data());

    //This is not allowed - no conversion
    //from VectorRefc back to VectorRef
    //vr1 = cvr2;
    }

SECTION("Vector to Ref Conversion")
    {
    auto size = 5;
    auto f1 = [](VectorRefc ref) { return ref.data(); };
    auto f2 = [](const VectorRefc& ref) { return ref.data(); };
    auto f3 = [](VectorRef ref) { ref *= 2; };

    Vector v = randomVec(size);
    auto origv = v;
    const Vector cv = randomVec(size);

    CHECK(f1(v) == v.data());
    CHECK(f1(cv) == cv.data());
    CHECK(f2(v) == v.data());
    CHECK(f2(cv) == cv.data());

    //This case not allowed:
    //CHECK(f3(cv) == cv.data());

    f3(v);
    for(auto j : count(size))
        CHECK_CLOSE(v(j),2*origv(j));

    //Not allowed to convert/assign to either
    //type of Ref from a temporary
    //f1(randomVec(size)); //not allowed
    //f3(randomVec(size)); //not allowed

    
    //Not allowed to construct a VectorRef from const Vec:
    //VectorRef vref1 = cv; //not allowed
    //VectorRef vref2(cv); //not allowed

    //Not allowed to assign/convert to VectorRef from const Vec:
    //VectorRef vref3;
    //vref3 = cv; //not allowed

    VectorRefc ref2(cv);
    CHECK(ref2.data() == cv.data());
    ref2.clear();
    ref2 = cv;
    CHECK(ref2.data() == cv.data());

    VectorRef ref3(v);
    CHECK(ref3.data() == v.data());
    ref3.clear();
    ref3 = v;
    CHECK(ref3.data() == v.data());

    VectorRefc ref4 = v;
    CHECK(ref4.data() == v.data());
    ref4.clear();
    ref4 = v;
    CHECK(ref4.data() == v.data());
    }

SECTION("Construct and assign from VectorRef")
    {
    auto size = 10;
    auto data = randomData(size);
    auto vr = makeVecRef(data.begin(),size);

    //print("data = "); for(auto& el : data) print(el," "); println();
    //print("vr = "); for(auto& el : vr) print(el," "); println();

    Vector v1(vr);
    CHECK(v1.size() == size);
    for(auto i : count(size))
        {
        CHECK_CLOSE(v1(i),data(i));
        }

    Vector v2(size+2);
    v2 = vr;
    CHECK(v2.size() == size);
    for(auto i : count(size))
        {
        CHECK_CLOSE(v2(i),data(i));
        }
    }

SECTION("Copy data")
    {
    auto size = 10;
    auto data1 = randomData(size);
    auto data2 = randomData(size);
    auto vr1 = makeVecRef(data1.begin(),size);
    auto vr2 = makeVecRef(data2.begin(),size);
    auto v1 = randomVec(size);
    auto v2 = randomVec(4*size);

    vr1 &= v1;
    for(auto i : count(size))
        CHECK_CLOSE(data1(i),v1(i));

    vr2 &= makeVecRef(v2.data(),size,4);
    for(auto i : count(size))
        CHECK_CLOSE(data2(i),v2(4*i));
    }

SECTION("Scalar multiply, divide")
    {
    auto size = 10;
    auto data = randomData(size);
    auto origdata = data;
    auto vr = makeVecRef(data.begin(),size);
    auto fac = Global::random();

    vr *= fac;
    for(auto i : count(size))
        {
        CHECK_CLOSE(data(i),fac*origdata(i));
        }

    Vector v1(vr);
    v1 *= fac;
    for(auto i : count(size))
        {
        CHECK_CLOSE(v1(i),fac*data(i));
        }
    }

SECTION("Test += -= operators")
    {
    auto size = 10;
    auto dataA = randomData(size);
    auto dataB = randomData(size);
    auto origdataA = dataA;
    auto A = makeVecRef(dataA.begin(),size);
    auto B = makeVecRef(dataB.begin(),size);

    A += B;
    for(auto i : count(size))
        {
        CHECK_CLOSE(A(i),B(i)+origdataA(i));
        }

    dataA = origdataA;
    auto cstride = 3;
    auto dataC = randomData(cstride*size);
    //Only access every third element:
    auto C = makeVecRef(dataC.begin(),size,cstride);
    A -= C;
    for(auto i : count(size))
        {
        CHECK_CLOSE(A(i),origdataA(i)-C(i));
        }
    }

SECTION("Vector + and -")
    {
    auto size = 20;
    auto v1 = randomVec(size),
         v2 = randomVec(size);
    auto origv1 = v1;

    v1 = v1+v2;
    for(auto i : count(size))
        {
        CHECK_CLOSE(v1(i),origv1(i)+v2(i));
        }

    v1 = origv1;
    v1 = v1-v2;
    for(auto i : count(size))
        {
        CHECK_CLOSE(v1(i),origv1(i)-v2(i));
        }

    v1 = origv1;
    auto origv2 = v2;
    v1 = v1+std::move(v2);
    CHECK(!v2);
    for(auto i : count(size))
        {
        CHECK_CLOSE(v1(i),origv1(i)+origv2(i));
        }

    v1 = origv1;
    CHECK(v1.data() != origv1.data());
    auto ref1 = makeRef(origv1);
    CHECK(ref1.data() == origv1.data());
    CHECK((std::is_same<decltype(ref1),VectorRef>::value));
    v1 += ref1;
    for(auto i : count(size))
        {
        CHECK_CLOSE(v1(i),2*origv1(i));
        }

    v1 = origv1;
    v1 = v1+ref1;
    for(auto i : count(size))
        {
        CHECK_CLOSE(v1(i),2*origv1(i));
        }

    v1 = origv1;
    auto cref1 = makeRefc(origv1);
    v1 = cref1+ref1;
    for(auto i : count(size))
        {
        CHECK_CLOSE(v1(i),2*origv1(i));
        }
    }

SECTION("Dot product")
    {
    auto N = 10;

    auto v = randomVec(N);
    auto vnrm = norm(v);
    CHECK_CLOSE(v*v,vnrm*vnrm);

    auto ref = makeRef(v);
    CHECK_CLOSE(ref*ref,vnrm*vnrm);
    CHECK_CLOSE(ref*v,vnrm*vnrm);
    CHECK_CLOSE(v*ref,vnrm*vnrm);

    const auto& cv = v;
    auto cref = makeRef(cv);
    CHECK_CLOSE(cref*cref,vnrm*vnrm);
    CHECK_CLOSE(ref*cref,vnrm*vnrm);

    ////Non-trivial stride case:
    //auto M1 = randomMat(N,N);
    //auto M2 = randomMat(N,N);
    //auto d1 = diagonal(M1);
    //auto d2 = diagonal(M2);

    //auto dot = d1*d2;
    //Real val = 0;
    //for(auto j : count(N))
    //    {
    //    val += M1(j,j)*M2(j,j);
    //    }
    //CHECK_CLOSE(val,dot);
    }

SECTION("Misc Vector Functions")
    {
    auto size = 20;
    auto v = Vector(size);

    //Test that randomize works
    randomize(v);

    //Test norm function
    Real calcnrm = 0;
    for(auto& el : v) calcnrm += el*el;
    calcnrm = std::sqrt(calcnrm);
    CHECK_CLOSE(norm(v),calcnrm);
    }


SECTION("Sub Vector")
    {
    auto v = randomVec(20);

    long start = 0,
         stop = 10;
    auto s = subVector(v,start,stop);
    size_t n = 0;
    for(auto j : count(s.size()))
        {
        CHECK_CLOSE(s(j),v(start+j));
        ++n;
        }
    CHECK(n == s.size());

    start = 2;
    stop = 12;
    s = subVector(v,start,stop);
    n = 0;
    for(auto j : count(s.size()))
        {
        CHECK_CLOSE(s(j),v(start+j));
        ++n;
        }
    CHECK(n == s.size());

    //Set elements of v through s
    for(auto& el : s) el = -1;
    for(auto j : count(start,stop))
        CHECK_CLOSE(-1,v(j));

    //Not allowed: would return ref to temporary
    //auto ref = subVector(Vector(20),1,10);
    }


SECTION("Test Resize")
    {
    auto size = 20;
    auto v = Vector(size);
    for(auto i : count(size)) v(i) = i*i;
    auto origv = v;

    //Check that downsizing doesn't change any elements
    resize(v,size/2);
    for(auto i : count(v.size())) 
        {
        CHECK(v(i) == i*i);
        }

    //Check that upsizing pads with zeros
    v = origv;
    resize(v,2*size);
    for(auto i : count(size)) 
        {
        CHECK(v(i) == i*i);
        }
    for(auto i : count(size+1,2*size)) 
        {
        CHECK_CLOSE(v(i),0);
        }
    }
}


TEST_CASE("Test MatrixRef")
{

SECTION("Constructors")
    {
    SECTION("Regular Constructor")
        {
        auto Ar = 3,
             Ac = 4;
        auto data = randomData(Ar*Ac);
        CHECK(data.size() == Ar*Ac);
        auto A = makeMatRef(data.begin(),data.size(),Ar,Ac);

        for(auto c : count(Ac))
        for(auto r : count(Ar))
            {
            CHECK_CLOSE(A(r,c),data[r+Ar*c]);
            }
        }

//    SECTION("Transpose Constructor")
//        {
//        auto Ar = 3,
//             Ac = 4;
//        auto data = randomData(1,Ar*Ac);
//        CHECK(data.size() == Ar*Ac);
//
//        //auto MC = MatrixRef(data.begin(),Ar,Ac,true);
//        //printfln("MC: %s (%d,%d,%d,%d)",MC.data(),MC.ind().rn,MC.ind().rs,MC.ind().cn,MC.ind().cs);
//
//        //auto MA = MatrixRef(data.begin(),Ar,Ac);
//        //MA.applyTrans();
//        //printfln("MA: %s (%d,%d,%d,%d)",MA.data(),MA.ind().rn,MA.ind().rs,MA.ind().cn,MA.ind().cs);
//
//        //auto MR = MatrixRefc(data.begin(),Ar,Ac);
//        //MR = MatrixRefc(MR.data(),MR.Nrows(),MR.Ncols(),true);
//        //printfln("MR: %s (%d,%d,%d,%d)",MR.data(),MR.ind().rn,MR.ind().rs,MR.ind().cn,MR.ind().cs);
//
//        auto A = MatrixRef(data.begin(),Ar,Ac,true);
//        CHECK(A.transposed());
//        for(auto r : count(Ar))
//        for(auto c : count(Ac))
//            {
//            CHECK_CLOSE(A(c,r),data[r+Ar*(c-1)]);
//            }
//
//        auto A2 = MatrixRef(data.begin(),Ar,Ac);
//        A2.applyTrans();
//        CHECK(A2.transposed());
//
//        for(auto r : count(Ar))
//        for(auto c : count(Ac))
//            {
//            CHECK_CLOSE(A2(c,r),data[r+Ar*(c-1)]);
//            }
//        }
    }

SECTION("Automatic Conversion")
    {
    auto N = 10;
    auto data1 = randomData(N*N);
    auto data2 = randomData(N*N);
    auto mr1 = makeMatRef(data1.begin(),data1.size(),N,N);
    auto cmr2 = makeMatRefc(data2.begin(),data2.size(),N,N);

    //Assigning to MatrixRefc from MatrixRef is ok
    cmr2 = mr1;
    CHECK(cmr2.data() == mr1.data());

    //Constructing VectorRefc from VectorRef is ok
    MatrixRefc cref(mr1);
    CHECK(cref.data() == mr1.data());

    //This is not allowed - no conversion
    //from MatrixRefc back to MatrixRef
    //mr1 = cmr2;
    }

SECTION("Test makeRef")
    {
    auto N = 10;
    auto data = randomData(N*N);

    auto Mr = makeMatRef(data.begin(),data.size(),N,N);
    auto cMr = makeMatRefc(data.begin(),data.size(),N,N);
    Matrix M(N,N);
    const auto& cM = M;

    auto ref1 = makeRef(Mr);
    CHECK(ref1.data() == data.begin());

    auto ref2 = makeRef(cMr);
    CHECK(ref2.data() == data.begin());
    CHECK((std::is_same<decltype(ref2),MatrixRefc>::value));

    auto ref3 = makeRef(M);
    CHECK(ref3.data() == M.data());
    CHECK((std::is_same<decltype(ref3),MatrixRef>::value));

    auto ref4 = makeRef(cM);
    CHECK(ref4.data() == cM.data());
    CHECK((std::is_same<decltype(ref4),MatrixRefc>::value));

    //Not allowed to make ref of temporary:
    //auto ref5 = makeRef(Matrix(N,N)); //not allowed
    //auto ref6 = makeRef(std::move(M)); //not allowed
    }



SECTION("Test += -= operators")
    {
    auto N = 5;
    auto dataA = randomData(N*N);
    auto dataB = randomData(N*N);
    auto origdataA = dataA;

    auto A = makeMatRef(dataA.begin(),dataA.size(),N,N);
    auto origA = makeMatRefc(origdataA.cbegin(),origdataA.size(),N,N);
    auto B = makeMatRefc(dataB.cbegin(),dataB.size(),N,N);
    auto At = transpose(makeMatRef(dataA.begin(),dataA.size(),N,N));
    auto Bt = transpose(makeMatRefc(dataB.cbegin(),dataB.size(),N,N));

    A += B;
    for(auto r : count(N))
    for(auto c : count(N))
        {
        CHECK_CLOSE(A(r,c),B(r,c)+origA(r,c));
        }

    dataA = origdataA;
    At += Bt;
    for(auto r : count(N))
    for(auto c : count(N))
        {
        CHECK_CLOSE(A(r,c),B(r,c)+origA(r,c));
        }

    dataA = origdataA;
    At += B;
    for(auto r : count(N))
    for(auto c : count(N))
        {
        CHECK_CLOSE(A(r,c),B(c,r)+origA(r,c));
        }

    dataA = origdataA;
    A -= Bt;
    for(auto r : count(N))
    for(auto c : count(N))
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
//    for(auto r : count(N)) 
//    for(auto c : count(N)) 
//        CHECK_CLOSE(M(r,c),2*origM(r,c));
//
//    Mr1 -= Mr2;
//    for(auto& el : M) CHECK(el < 1E-10);
//    }

SECTION("Test MatrixRef mult")
    {
    SECTION("Case 1")
        {
        auto Ar = 3,
             K  = 4,
             Bc = 5;
        auto dataA = randomData(Ar*K);
        auto dataB = randomData(K*Bc);
        auto dataC = autovector<Real>(1,Ar*Bc);

        auto A = makeMatRef(dataA.begin(),dataA.size(),Ar,K);
        auto B = makeMatRef(dataB.begin(),dataB.size(),K,Bc);
        auto C = makeMatRef(dataC.begin(),dataC.size(),Ar,Bc);

        mult(A,B,C);
        for(auto r : count(nrows(C)))
        for(auto c : count(ncols(C)))
            {
            Real val = 0;
            for(auto k : count(K)) val += A(r,k)*B(k,c);
            CHECK_CLOSE(C(r,c),val);
            }
        }

    SECTION("Case 2")
        {
        auto Ac = 3,
             K  = 4,
             Bc = 5;
        auto dataA = randomData(K*Ac);
        auto dataB = randomData(K*Bc);
        auto dataC = autovector<Real>(1,Ac*Bc);

        auto A = makeMatRef(dataA.begin(),dataA.size(),K,Ac);
        auto B = makeMatRef(dataB.begin(),dataB.size(),K,Bc);
        auto C = makeMatRef(dataC.begin(),dataC.size(),Ac,Bc);

        auto At = transpose(A);
        mult(At,B,C);
        for(auto r : count(nrows(C)))
        for(auto c : count(ncols(C)))
            {
            Real val = 0;
            for(auto k : count(K)) val += At(r,k)*B(k,c);
            CHECK_CLOSE(C(r,c),val);
            }
        }

    SECTION("Case 3")
        {
        auto Ar = 3,
             K =  4,
             Br = 5;
        auto dataA = randomData(Ar*K);
        auto dataB = randomData(Br*K);
        auto dataC = autovector<Real>(1,Ar*Br);

        auto A = makeMatRef(dataA.begin(),dataA.size(),Ar,K);
        auto B = makeMatRef(dataB.begin(),dataB.size(),Br,K);
        auto C = makeMatRef(dataC.begin(),dataC.size(),Ar,Br);

        auto Bt = transpose(B);
        mult(A,Bt,C);
        for(auto r : count(nrows(C)))
        for(auto c : count(ncols(C)))
            {
            Real val = 0;
            for(auto k : count(K)) val += A(r,k)*Bt(k,c);
            CHECK_CLOSE(C(r,c),val);
            }
        }

    SECTION("Case 4")
        {
        auto Ac = 3,
             K =  4,
             Br = 5;
        auto dataA = randomData(K*Ac);
        auto dataB = randomData(Br*K);
        auto dataC = autovector<Real>(1,Ac*Br);

        auto A = makeMatRef(dataA.begin(),dataA.size(),K,Ac);
        auto B = makeMatRef(dataB.begin(),dataB.size(),Br,K);
        auto C = makeMatRef(dataC.begin(),dataC.size(),Ac,Br);

        auto At = transpose(A);
        auto Bt = transpose(B);
        mult(At,Bt,C);
        for(auto r : count(nrows(C)))
        for(auto c : count(ncols(C)))
            {
            Real val = 0;
            for(auto k : count(K)) val += At(r,k)*Bt(k,c);
            CHECK_CLOSE(C(r,c),val);
            }
        }

    SECTION("Case 5")
        {
        auto Ar = 3,
             K =  4,
             Bc = 5;
        auto dataA = randomData(Ar*K);
        auto dataB = randomData(K*Bc);
        auto dataC = autovector<Real>(1,Ar*Bc);

        auto A = makeMatRef(dataA.begin(),dataA.size(),Ar,K);
        auto B = makeMatRef(dataB.begin(),dataB.size(),K,Bc);
        auto C = makeMatRef(dataC.begin(),dataC.size(),Bc,Ar);

        auto Ct = transpose(C);
        mult(A,B,Ct);
        for(auto r : count(nrows(Ct)))
        for(auto c : count(ncols(Ct)))
            {
            Real val = 0;
            for(auto k : count(K)) val += A(r,k)*B(k,c);
            CHECK_CLOSE(Ct(r,c),val);
            }
        }

    SECTION("Case 6")
        {
        auto Ac = 3,
             K =  4,
             Bc = 5;
        auto dataA = randomData(K*Ac);
        auto dataB = randomData(K*Bc);
        auto dataC = autovector<Real>(1,Bc*Ac);

        auto A = makeMatRef(dataA.begin(),dataA.size(),K,Ac);
        auto B = makeMatRef(dataB.begin(),dataB.size(),K,Bc);
        auto C = makeMatRef(dataC.begin(),dataC.size(),Bc,Ac);

        auto At = transpose(A);
        auto Ct = transpose(C);
        mult(At,B,Ct);
        for(auto r : count(nrows(Ct)))
        for(auto c : count(ncols(Ct)))
            {
            Real val = 0;
            for(auto k : count(K)) val += At(r,k)*B(k,c);
            CHECK_CLOSE(Ct(r,c),val);
            }
        }

    SECTION("Case 7")
        {
        auto Ar = 3,
             K =  4,
             Br = 5;
        auto dataA = randomData(Ar*K);
        auto dataB = randomData(Br*K);
        auto dataC = autovector<Real>(1,Br*Ar);

        auto A = makeMatRef(dataA.begin(),dataA.size(),Ar,K);
        auto B = makeMatRef(dataB.begin(),dataB.size(),Br,K);
        auto C = makeMatRef(dataC.begin(),dataC.size(),Br,Ar);

        auto Bt = transpose(B);
        auto Ct = transpose(C);
        mult(A,Bt,Ct);
        for(auto r : count(nrows(Ct)))
        for(auto c : count(ncols(Ct)))
            {
            Real val = 0;
            for(auto k : count(K)) val += A(r,k)*Bt(k,c);
            CHECK_CLOSE(Ct(r,c),val);
            }
        }

    SECTION("Case 8")
        {
        auto Ac = 3,
             K =  4,
             Br = 5;
        auto dataA = randomData(K*Ac);
        auto dataB = randomData(Br*K);
        auto dataC = autovector<Real>(1,Br*Ac);

        auto A = makeMatRef(dataA.begin(),dataA.size(),K,Ac);
        auto B = makeMatRef(dataB.begin(),dataB.size(),Br,K);
        auto C = makeMatRef(dataC.begin(),dataC.size(),Br,Ac);

        auto At = transpose(A);
        auto Bt = transpose(B);
        auto Ct = transpose(C);
        mult(At,Bt,Ct);
        for(auto r : count(nrows(Ct)))
        for(auto c : count(ncols(Ct)))
            {
            Real val = 0;
            for(auto k : count(K)) val += At(r,k)*Bt(k,c);
            CHECK_CLOSE(Ct(r,c),val);
            }
        }

    SECTION("Multipy with self")
        {
        auto N = 8;
        auto dataA = randomData(N*N);
        auto dataC = randomData(N*N);

        auto A = makeMatRef(dataA.begin(),dataA.size(),N,N);
        auto C = makeMatRef(dataC.begin(),dataC.size(),N,N);

        mult(A,A,C);
        for(auto r : count(nrows(C)))
        for(auto c : count(ncols(C)))
            {
            Real val = 0;
            for(auto k : count(N)) val += A(r,k)*A(k,c);
            CHECK_CLOSE(C(r,c),val);
            }
        }
    }

SECTION("Test multAdd")
    {
    auto Ar = 3,
         K  = 4,
         Bc = 5;
    auto dataA = randomData(Ar*K);
    auto dataB = randomData(K*Bc);
    auto dataC = autovector<Real>(1,Ar*Bc);

    auto A = makeMatRef(dataA.begin(),dataA.size(),Ar,K);
    auto B = makeMatRef(dataB.begin(),dataB.size(),K,Bc);
    auto C = makeMatRef(dataC.begin(),dataC.size(),Ar,Bc);

    //Save a copy of C's original data in order
    //to explicitly carry out multAdd alg. below
    auto orig_dataC = dataC;
    auto origC = makeMatRef(orig_dataC.begin(),orig_dataC.size(),Ar,Bc);

    multAdd(A,B,C);
    for(auto r : count(nrows(C)))
    for(auto c : count(ncols(C)))
        {
        Real val = 0;
        for(auto k : count(K)) val += A(r,k)*B(k,c) + origC(r,c);
        CHECK_CLOSE(C(r,c),val);
        }
    }


} //Test MatrixRef


TEST_CASE("Test Matrix")
{
SECTION("Constructors")
    {
    SECTION("Constructor Basics")
        {
        auto Ar = 3,
             Ac = 4;
        auto A = Matrix(Ar,Ac);
        CHECK(nrows(A) == Ar);
        CHECK(ncols(A) == Ac);
        CHECK(A);

        Matrix B;
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
            CHECK_CLOSE(A(r,c),data[r+Ar*c]);
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

    for(auto r : count(N))
    for(auto c : count(N))
        CHECK_CLOSE(M2(r,c),M1t(r,c));

    auto M1tt = transpose(transpose(M1));
    M2 = M1tt;
    for(auto r : count(N))
    for(auto c : count(N))
        CHECK_CLOSE(M2(r,c),M1tt(r,c));
    }

SECTION("Test Matrix multiplication")
    {
    auto Ar = 3,
         K  = 4,
         Bc = 5;

    //Multiply matrices A*B, store in matrix C
    auto A = Matrix(Ar,K);
    auto B = Matrix(K,Bc);
    auto C = Matrix(Ar,Bc);
    mult(A,B,C);
    for(auto r : count(nrows(C)))
    for(auto c : count(ncols(C)))
        {
        Real val = 0;
        for(auto k : count(K)) val += A(r,k)*B(k,c);
        CHECK_CLOSE(C(r,c),val);
        }


    //Store result in a matrixref instead
    auto dataC = autovector<Real>(1,Ar*Bc);
    auto Cref = makeMatRef(dataC.begin(),dataC.size(),Ar,Bc);
    mult(A,B,Cref);
    for(auto r : count(nrows(C)))
    for(auto c : count(ncols(C)))
        {
        Real val = 0;
        for(auto k : count(K)) val += A(r,k)*B(k,c);
        CHECK_CLOSE(Cref(r,c),val);
        }

    //Use operator*
    auto R = A*B;
    for(auto r : count(nrows(C)))
    for(auto c : count(ncols(C)))
        {
        Real val = 0;
        for(auto k : count(K)) val += A(r,k)*B(k,c);
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
    for(auto r : count(Nr))
    for(auto c : count(Nc))
        CHECK_CLOSE(C(r,c),A(r,c)+B(r,c));

    C = A;
    C -= B;
    for(auto r : count(Nr))
    for(auto c : count(Nc))
        CHECK_CLOSE(C(r,c),A(r,c)-B(r,c));

    auto D = B;
    C = A + std::move(D);
    CHECK(!D);
    CHECK(B);
    for(auto r : count(Nr))
    for(auto c : count(Nc))
        CHECK_CLOSE(C(r,c),A(r,c)+B(r,c));

    D = B;
    CHECK(D);
    C = A - std::move(D);
    CHECK(!D);
    CHECK(B);
    for(auto r : count(Nr))
    for(auto c : count(Nc))
        CHECK_CLOSE(C(r,c),A(r,c)-B(r,c));
    }

SECTION("Scalar multiply, divide")
    {
    auto N = 10;
    auto A = randomMat(N,N);
    auto origA = A;
    auto fac = Global::random();

    A *= fac;
    for(auto r : count(N))
    for(auto c : count(N))
        CHECK_CLOSE(A(r,c),origA(r,c)*fac);

    A = origA;
    A /= fac;
    for(auto r : count(N))
    for(auto c : count(N))
        CHECK_CLOSE(A(r,c),origA(r,c)/fac);

    A = origA;
    auto At = transpose(A);
    At *= fac;
    for(auto r : count(N))
    for(auto c : count(N))
        CHECK_CLOSE(A(r,c),origA(r,c)*fac);

    A = origA;
    auto S = subMatrix(A,0,N/2,0,N/2);
    S *= fac;
    for(auto r : count(N/2))
    for(auto c : count(N/2))
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
        for(auto r : count(N))
            {
            Real val = 0;
            for(auto c : count(N)) val += M(r,c)*x(c);
            CHECK_CLOSE(y(r),val);
            }

        y = transpose(M)*x;
        for(auto r : count(N))
            {
            Real val = 0;
            for(auto c : count(N)) val += M(c,r)*x(c);
            CHECK_CLOSE(y(r),val);
            }

        y = x*M;
        for(auto c : count(N))
            {
            Real val = 0;
            for(auto r : count(N)) val += x(r)*M(r,c);
            CHECK_CLOSE(y(c),val);
            }

        y = x*transpose(M);
        for(auto c : count(N))
            {
            Real val = 0;
            for(auto r : count(N)) val += x(r)*M(c,r);
            CHECK_CLOSE(y(c),val);
            }

        //Check a case where vector is strided
        auto d = diagonal(M);
        y = M*d;
        for(auto r : count(N))
            {
            Real val = 0;
            for(auto c : count(N)) val += M(r,c)*M(c,c);
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
        for(auto r : count(Ar))
            {
            Real val = 0;
            for(auto c : count(Ac)) val += A(r,c)*vR(c);
            CHECK_CLOSE(res(r),val);
            }

        res = vL*A;
        for(auto c : count(Ac))
            {
            Real val = 0;
            for(auto r : count(Ar)) val += vL(r)*A(r,c);
            CHECK_CLOSE(res(c),val);
            }

        res = transpose(A)*vL;
        for(auto c : count(Ac))
            {
            Real val = 0;
            for(auto r : count(Ar)) val += vL(r)*A(r,c);
            CHECK_CLOSE(res(c),val);
            }

        res = vR*transpose(A);
        for(auto r : count(Ar))
            {
            Real val = 0;
            for(auto c : count(Ac)) val += A(r,c)*vR(c);
            CHECK_CLOSE(res(r),val);
            }

        } //Rectangular case
    }

SECTION("Test reduceColsTo")
    {
    auto Nr = 10,
         Nc = 10;
    auto M = randomMat(Nr,Nc);
    auto origM = M;

    reduceCols(M,Nc/2);
    for(auto r : count(Nr))
    for(auto c : count(Nc/2))
        {
        CHECK_CLOSE(M(r,c),origM(r,c));
        }
    }

} // Test matrix


TEST_CASE("Test slicing")
{


SECTION("Diagonal")
    {
    auto A = randomMat(4,4);
    auto d = diagonal(A);
    long i = 0;
    for(const auto& el : d)
        {
        CHECK_CLOSE(el,A(i,i));
        ++i;
        }

    //Set diag els of A to 100 through d
    for(auto& el : d) el = 100;

    for(auto j : count(d.size()))
        {
        CHECK_CLOSE(100,A(j,j));
        }

    A = randomMat(5,10);
    d = diagonal(A);
    i = 0;
    for(const auto& el : d)
        {
        CHECK_CLOSE(el,A(i,i));
        ++i;
        }

    A = randomMat(10,5);
    d = diagonal(A);
    i = 0;
    for(const auto& el : d)
        {
        CHECK_CLOSE(el,A(i,i));
        ++i;
        }
    for(auto j : count(d.size()))
        {
        CHECK_CLOSE(d(j),A(j,j));
        }
    }

SECTION("Row / Col Slicing")
    {
    auto N = 10;
    auto A = randomMat(N,N);

    for(auto r : count(N))
        {
        auto R = row(A,r);
        for(auto c : count(N))
            CHECK_CLOSE(A(r,c),R(c));
        }

    for(auto c : count(N))
        {
        auto C = column(A,c);
        for(auto r : count(N))
            CHECK_CLOSE(A(r,c),C(r));
        }

    auto S = transpose(subMatrix(A,0,N/2,0,N/2));
    for(auto c : count(N/2))
        {
        auto C = column(S,c);
        for(auto r : count(N/2))
            CHECK_CLOSE(S(r,c),C(r));
        }
    }

SECTION("Transpose")
    {
    auto nr = 10,
         nc = 15;
    auto A = randomMat(nr,nc);
    auto At = transpose(A);
    CHECK(!isTransposed(A));
    CHECK(isTransposed(At));

    for(auto i : count(nr))
    for(auto j : count(nc))
        {
        CHECK_CLOSE(A(i,j),At(j,i));
        }

    auto B = randomMat(nc,nr);
    transpose(B) &= A;
    for(auto i : count(nc))
    for(auto j : count(nr))
        {
        CHECK_CLOSE(B(i,j),A(j,i));
        }
    }


SECTION("Sub Matrix")
    {
    auto nr = 10,
         nc = 15;
    auto A = randomMat(nr,nc);

    auto rstart = 0,
         rstop = 3,
         cstart = 0,
         cstop = 4;
    auto S = subMatrix(A,rstart,rstop,cstart,cstop);
    for(auto r : count(nrows(S)))
    for(auto c : count(ncols(S)))
        {
        CHECK_CLOSE(S(r,c),A(rstart+r,cstart+c));
        }

    rstart = 1;
    cstart = 2;
    auto At = transpose(A);
    CHECK(isTransposed(At));
    S = subMatrix(At,rstart,rstop,cstart,cstop);
    for(auto r : count(nrows(S)))
    for(auto c : count(ncols(S)))
        {
        CHECK_CLOSE(S(r,c),At(rstart+r,cstart+c));
        }

    rstart = 2;
    rstop = 8;
    cstart = 4;
    cstop = 10;
    S = subMatrix(A,rstart,rstop,cstart,cstop);
    for(auto r : count(nrows(S)))
    for(auto c : count(ncols(S)))
        {
        CHECK_CLOSE(S(r,c),A(rstart+r,cstart+c));
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

    Matrix U;
    Vector d;
    diagSymmetric(M,U,d);

    auto D = Matrix(N,N);
    diagonal(D) &= d;
    auto R = U*D*transpose(U);

    for(auto r : count(N))
    for(auto c : count(N))
        {
        CHECK_CLOSE(R(r,c),M(r,c));
        }
    CHECK(norm(R-M) < 1E-12*norm(M));
    }

SECTION("Orthogonalize")
    {
    auto N = 4;
    auto M = randomMat(N,N);

    orthog(M);

    auto R = transpose(M)*M;
    for(auto r : count(N))
    for(auto c : count(N))
        {
        if(r == c) CHECK_CLOSE(R(r,c),1);
        else       CHECK(R(r,c) < 1E-12);
        }
    }

SECTION("Singular Value Decomp")
    {
    Matrix U,V,D;
    Vector d;

    auto dMat = 
        [](Vector const& d) 
        { 
        Matrix D(d.size(),d.size()); 
        diagonal(D) &= d; 
        return D; 
        };

    auto Nr = 10,
         Nc = 8;
    auto M = randomMat(Nr,Nc);
    SVD(M,U,d,V);
    auto R = U*dMat(d)*transpose(V);
    CHECK((norm(R-M)/norm(M)) < 1E-13);

    Nr = 10;
    Nc = 10;
    M = randomMat(Nr,Nc);
    SVD(M,U,d,V);
    auto R2 = U*dMat(d)*transpose(V);
    CHECK((norm(R2-M)/norm(M)) < 1E-13);

    Nr = 10;
    Nc = 20;
    M = randomMat(Nr,Nc);
    SVD(M,U,d,V);
    auto R3 = U*dMat(d)*transpose(V);
    CHECK((norm(R3-M)/norm(M)) < 1E-13);
    }
} //Test matrix algorithms
