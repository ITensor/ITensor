#include "test.h"
#include "itensor.h"
#include "contract.h"

using namespace itensor;

struct GetFront
    {
    const Real* front;

    GetFront() {}

    NewData
    operator()(const ITDense<Real>& d)
        {
        front = d.data.data();
        return NewData();
        }

    NewData
    operator()(const ITDiag<Real>& d)
        {
        Error("Encountered Diag ITensor");
        return NewData();
        }

    NewData
    operator()(const ITDense<Complex>& d)
        {
        Error("Encountered Complex ITensor");
        return NewData();
        }

    NewData
    operator()(const ITDiag<Complex>& d)
        {
        Error("Encountered Complex Diag ITensor");
        return NewData();
        }

    //template<typename T>
    //NewData
    //operator()(T& d) 
    //    { 
    //    Error("GetFront not supported"); 
    //    return NewData(); 
    //    }
    };

void 
simple_from_ITensor(const ITensor& it, RTensor& st)
    {
    if(isComplex(it)) error("Complex case not implemented yet");
    auto r = it.r();
    std::vector<long> n(r);
    const auto& is = it.indices();
    for(int i = 0; i < r; ++i)
        n[i] = is[i].m();
    //PRI(n)
    auto front = applyFunc<GetFront>(it.data()).front;
    st = RTensor(front,std::move(n));
    }

// C += A * B
void 
mult_add(ITensor& A, 
         ITensor& B, 
         ITensor& C,
         const Args& args)
    {
    //Convert ITensors to RTensors
    RTensor a,b,c;
    println("Calling simple_from_ITensor A");
    simple_from_ITensor(A,a);
    println("Calling simple_from_ITensor B");
    simple_from_ITensor(B,b);
    println("Calling simple_from_ITensor C");
    simple_from_ITensor(C,c);

    //Make a small_map which maps each unique Index
    //to a unique integer
    int i = 0;
    small_map<Index,int> imap;
    for(auto ind : A.inds())
        if(imap[ind] == 0) imap[ind] = ++i;
    for(auto ind : B.inds())
        if(imap[ind] == 0) imap[ind] = ++i;

    //Create vectors of integers instructing
    //contractloop how to contract A with B
    Label ai(A.r()),
          bi(B.r()),
          ci(C.r());
    for(int j = 0; j < A.r(); ++j)
        ai.at(j) = imap[A.inds()[j]];
    for(int j = 0; j < B.r(); ++j)
        bi.at(j) = imap[B.inds()[j]];
    for(int j = 0; j < C.r(); ++j)
        ci.at(j) = imap[C.inds()[j]];

    contractloop(a,ai,b,bi,c,ci,args);
    }

TEST_CASE("Contract Test")
    {

    //SECTION("Contract Reshape Basic")
    //    {
    //    RTensor A(2,2),
    //            B(2,2),
    //            C(2,2);
    //    A(0,0) = 1; A(0,1) = 2;
    //    A(1,0) = 3; A(1,1) = 4;

    //    B(0,0) = 5; B(0,1) = 6;
    //    B(1,0) = 7; B(1,1) = 8;
    //    //
    //    // 1 2  5 6   19 22
    //    // 3 4  7 8   43 50
    //    //

    //    contract_reshape(A,{1,2},B,{2,3},C,{1,3});
    //    for(int r = 0; r < 2; ++r)
    //    for(int c = 0; c < 2; ++c)
    //        {
    //        Real val = 0;
    //        for(int k = 0; k < 2; ++k)
    //            {
    //            //printfln("A(%d,%d)*B(%d,%d)=%d*%d=%d",r,k,k,c,A(r,k),B(k,c),A(r,k)*B(k,c));
    //            val += A(r,k)*B(k,c);
    //            }
    //        REQUIRE(C(r,c) == val);
    //        }

    //    contract_reshape(A,{1,2},B,{3,2},C,{1,3});
    //    for(int r = 0; r < 2; ++r)
    //    for(int c = 0; c < 2; ++c)
    //        {
    //        Real val = 0;
    //        for(int k = 0; k < 2; ++k)
    //            {
    //            val += A(r,k)*B(c,k);
    //            }
    //        REQUIRE(C(r,c) == val);
    //        }

    //    contract_reshape(A,{1,2},B,{2,3},C,{3,1});
    //    for(int r = 0; r < 2; ++r)
    //    for(int c = 0; c < 2; ++c)
    //        {
    //        Real val = 0;
    //        for(int k = 0; k < 2; ++k)
    //            {
    //            val += A(r,k)*B(k,c);
    //            }
    //        REQUIRE(C(c,r) == val);
    //        }

    //    contract_reshape(A,{1,2},B,{3,2},C,{3,1});
    //    for(int r = 0; r < 2; ++r)
    //    for(int c = 0; c < 2; ++c)
    //        {
    //        Real val = 0;
    //        for(int k = 0; k < 2; ++k)
    //            {
    //            val += A(r,k)*B(c,k);
    //            }
    //        REQUIRE(C(c,r) == val);
    //        }

    //    contract_reshape(A,{2,1},B,{2,3},C,{1,3});
    //    for(int r = 0; r < 2; ++r)
    //    for(int c = 0; c < 2; ++c)
    //        {
    //        Real val = 0;
    //        for(int k = 0; k < 2; ++k)
    //            {
    //            val += A(k,r)*B(k,c);
    //            }
    //        REQUIRE(C(r,c) == val);
    //        }

    //    contract_reshape(A,{2,1},B,{3,2},C,{1,3});
    //    for(int r = 0; r < 2; ++r)
    //    for(int c = 0; c < 2; ++c)
    //        {
    //        Real val = 0;
    //        for(int k = 0; k < 2; ++k)
    //            {
    //            val += A(k,r)*B(c,k);
    //            }
    //        REQUIRE(C(r,c) == val);
    //        }

    //    contract_reshape(A,{2,1},B,{2,3},C,{3,1});
    //    for(int r = 0; r < 2; ++r)
    //    for(int c = 0; c < 2; ++c)
    //        {
    //        Real val = 0;
    //        for(int k = 0; k < 2; ++k)
    //            {
    //            val += A(k,r)*B(k,c);
    //            }
    //        REQUIRE(C(c,r) == val);
    //        }

    //    contract_reshape(A,{2,1},B,{3,2},C,{3,1});
    //    for(int r = 0; r < 2; ++r)
    //    for(int c = 0; c < 2; ++c)
    //        {
    //        Real val = 0;
    //        for(int k = 0; k < 2; ++k)
    //            {
    //            val += A(k,r)*B(c,k);
    //            }
    //        REQUIRE(C(c,r) == val);
    //        }
    //    }

    SECTION("Contract Reshape")
        {
        auto randomize = [](auto& t)
            {
            for(auto& elt : t) elt = Global::random();
            };

        SECTION("Case 1")
            {
            RTensor A(2,3,4),
                    B(3,7,2),
                    C;
            randomize(A);
            randomize(B);
            contract_reshape(A,{2,3,4},B,{3,7,2},C,{4,7});
            for(int i4 = 0; i4 < 4; ++i4)
            for(int i7 = 0; i7 < 7; ++i7)
                {
                Real val = 0;
                for(int i2 = 0; i2 < 2; ++i2)
                for(int i3 = 0; i3 < 3; ++i3)
                    {
                    val += A(i2,i3,i4)*B(i3,i7,i2);
                    }
                CHECK_REQUAL(C(i4,i7),val);
                }
            }

        SECTION("Case 2")
            {
            RTensor A(2,3,4),
                    B(3,7,2),
                    C;
            randomize(A);
            randomize(B);
            contract_reshape(A,{2,3,4},B,{3,7,2},C,{7,4});
            for(int i4 = 0; i4 < 4; ++i4)
            for(int i7 = 0; i7 < 7; ++i7)
                {
                Real val = 0;
                for(int i2 = 0; i2 < 2; ++i2)
                for(int i3 = 0; i3 < 3; ++i3)
                    {
                    val += A(i2,i3,i4)*B(i3,i7,i2);
                    }
                CHECK_REQUAL(C(i7,i4),val);
                }
            }

        SECTION("Case 3")
            {
            RTensor A(2,4,3),
                    B(3,7,2),
                    C;
            randomize(A);
            randomize(B);
            contract_reshape(A,{2,4,3},B,{3,7,2},C,{7,4});
            for(int i4 = 0; i4 < 4; ++i4)
            for(int i7 = 0; i7 < 7; ++i7)
                {
                Real val = 0;
                for(int i2 = 0; i2 < 2; ++i2)
                for(int i3 = 0; i3 < 3; ++i3)
                    {
                    val += A(i2,i4,i3)*B(i3,i7,i2);
                    }
                CHECK_REQUAL(C(i7,i4),val);
                }
            }

        SECTION("Case 4")
            {
            RTensor A(2,4,3),
                    B(3,7,2),
                    C;
            randomize(A);
            randomize(B);
            contract_reshape(A,{2,4,3},B,{3,7,2},C,{4,7});
            for(int i4 = 0; i4 < 4; ++i4)
            for(int i7 = 0; i7 < 7; ++i7)
                {
                Real val = 0;
                for(int i2 = 0; i2 < 2; ++i2)
                for(int i3 = 0; i3 < 3; ++i3)
                    {
                    val += A(i2,i4,i3)*B(i3,i7,i2);
                    }
                CHECK_REQUAL(C(i4,i7),val);
                }
            }
        }

    //SECTION("Contract Loop")
    //    {
    //    int fac = 1;
    //    Index i1("i1",200*fac), 
    //          i2("i2",200*fac), 
    //          i3("i3",200*fac), 
    //          i4("i4",30), 
    //          i5("i5",10),
    //          i6("i6",10);
    //    auto itA = randIT(i1,i2,i4,i5);
    //    Print(itA);
    //    auto itB = randIT(i1,i3,i4,i6); 
    //    Print(itB);
    //    auto itC = randIT(i2,i3,i5,i6);
    //    Print(itC);
    //    mult_add(itA, itB, itC,{"NThread",4});
    //    }
    }
