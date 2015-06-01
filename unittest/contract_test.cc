#include "test.h"
#include "itensor/util/cputime.h"
#include "itensor/util/count.h"
#include "itensor/tensor/contract.h"
#include "itensor/util/set_scoped.h"

using namespace itensor;

TEST_CASE("Contract Test")
    {
    auto randomize = [](auto& t)
        {
        for(auto& elt : t) elt = Global::random();
        };

    SECTION("Contract Reshape Basic")
        {
        Ten A(2,2),
            B(2,2),
            C(2,2);
        A(0,0) = 1; A(0,1) = 2;
        A(1,0) = 3; A(1,1) = 4;

        B(0,0) = 5; B(0,1) = 6;
        B(1,0) = 7; B(1,1) = 8;

        SECTION("Case 1")
            {
            //
            // 1 2  5 6   19 22
            // 3 4  7 8   43 50
            //
            contract(A,{1,2},B,{2,3},C,{1,3});
            for(int r = 0; r < 2; ++r)
            for(int c = 0; c < 2; ++c)
                {
                Real val = 0;
                for(int k = 0; k < 2; ++k)
                    {
                    //printfln("A(%d,%d)*B(%d,%d)=%d*%d=%d",r,k,k,c,A(r,k),B(k,c),A(r,k)*B(k,c));
                    val += A(r,k)*B(k,c);
                    }
                REQUIRE(C(r,c) == val);
                }
            }

        SECTION("Case 2")
            {
            //
            // 1 2  5 7   17 23
            // 3 4  6 8   39 53
            //
            //println("Case 2:");
            contract(A,{1,2},B,{3,2},C,{1,3});
            for(int r = 0; r < 2; ++r)
            for(int c = 0; c < 2; ++c)
                {
                Real val = 0;
                for(int k = 0; k < 2; ++k)
                    {
                    val += A(r,k)*B(c,k);
                    }
                REQUIRE(C(r,c) == val);
                }
            }

        SECTION("Case 3")
            {
            contract(A,{1,2},B,{2,3},C,{3,1});
            for(int r = 0; r < 2; ++r)
            for(int c = 0; c < 2; ++c)
                {
                Real val = 0;
                for(int k = 0; k < 2; ++k)
                    {
                    val += A(r,k)*B(k,c);
                    }
                REQUIRE(C(c,r) == val);
                }
            }

        SECTION("Case 4")
            {
            contract(A,{1,2},B,{3,2},C,{3,1});
            for(int r = 0; r < 2; ++r)
            for(int c = 0; c < 2; ++c)
                {
                Real val = 0;
                for(int k = 0; k < 2; ++k)
                    {
                    val += A(r,k)*B(c,k);
                    }
                REQUIRE(C(c,r) == val);
                }
            }

        SECTION("Case 5")
            {
            contract(A,{2,1},B,{2,3},C,{1,3});
            for(int r = 0; r < 2; ++r)
            for(int c = 0; c < 2; ++c)
                {
                Real val = 0;
                for(int k = 0; k < 2; ++k)
                    {
                    val += A(k,r)*B(k,c);
                    }
                REQUIRE(C(r,c) == val);
                }
            }

        SECTION("Case 6")
            {
            contract(A,{2,1},B,{3,2},C,{1,3});
            for(int r = 0; r < 2; ++r)
            for(int c = 0; c < 2; ++c)
                {
                Real val = 0;
                for(int k = 0; k < 2; ++k)
                    {
                    val += A(k,r)*B(c,k);
                    }
                REQUIRE(C(r,c) == val);
                }
            }

        SECTION("Case 7")
            {
            contract(A,{2,1},B,{2,3},C,{3,1});
            for(int r = 0; r < 2; ++r)
            for(int c = 0; c < 2; ++c)
                {
                Real val = 0;
                for(int k = 0; k < 2; ++k)
                    {
                    val += A(k,r)*B(k,c);
                    }
                REQUIRE(C(c,r) == val);
                }
            }

        SECTION("Case 8")
            {
            contract(A,{2,1},B,{3,2},C,{3,1});
            for(int r = 0; r < 2; ++r)
            for(int c = 0; c < 2; ++c)
                {
                Real val = 0;
                for(int k = 0; k < 2; ++k)
                    {
                    val += A(k,r)*B(c,k);
                    }
                REQUIRE(C(c,r) == val);
                }
            }
        }

    SECTION("Contract Reshape Non-Matrix")
        {
        SECTION("Case 1")
            {
            Ten A(2,3,4),
                B(3,7,2),
                C(4,7);
            randomize(A);
            randomize(B);
            contract(A,{2,3,4},B,{3,7,2},C,{4,7});
            for(int i4 = 0; i4 < 4; ++i4)
            for(int i7 = 0; i7 < 7; ++i7)
                {
                Real val = 0;
                for(int i2 = 0; i2 < 2; ++i2)
                for(int i3 = 0; i3 < 3; ++i3)
                    {
                    val += A(i2,i3,i4)*B(i3,i7,i2);
                    }
                CHECK_CLOSE(C(i4,i7),val);
                }
            }

        SECTION("Case 2")
            {
            Ten A(2,3,4),
                B(3,7,2),
                C(7,4);
            randomize(A);
            randomize(B);
            Global::debug3() = true;
            contract(A,{2,3,4},B,{3,7,2},C,{7,4});
            Global::debug3() = false;
            for(int i4 = 0; i4 < 4; ++i4)
            for(int i7 = 0; i7 < 7; ++i7)
                {
                Real val = 0;
                for(int i2 = 0; i2 < 2; ++i2)
                for(int i3 = 0; i3 < 3; ++i3)
                    {
                    val += A(i2,i3,i4)*B(i3,i7,i2);
                    }
                CHECK_CLOSE(C(i7,i4),val);
                }
            }

        SECTION("Case 3")
            {
            Ten A(2,4,3),
                B(3,7,2),
                C(7,4);
            randomize(A);
            randomize(B);
            contract(A,{2,4,3},B,{3,7,2},C,{7,4});
            for(int i4 = 0; i4 < 4; ++i4)
            for(int i7 = 0; i7 < 7; ++i7)
                {
                Real val = 0;
                for(int i2 = 0; i2 < 2; ++i2)
                for(int i3 = 0; i3 < 3; ++i3)
                    {
                    val += A(i2,i4,i3)*B(i3,i7,i2);
                    }
                CHECK_CLOSE(C(i7,i4),val);
                }
            }

        SECTION("Case 4")
            {
            Ten A(2,4,3),
                B(3,7,2),
                C(4,7);
            randomize(A);
            randomize(B);
            contract(A,{2,4,3},B,{3,7,2},C,{4,7});
            for(int i4 = 0; i4 < 4; ++i4)
            for(int i7 = 0; i7 < 7; ++i7)
                {
                Real val = 0;
                for(int i2 = 0; i2 < 2; ++i2)
                for(int i3 = 0; i3 < 3; ++i3)
                    {
                    val += A(i2,i4,i3)*B(i3,i7,i2);
                    }
                CHECK_CLOSE(C(i4,i7),val);
                }
            }

        SECTION("Case NM1")
            {
            Ten A(2,3,4),
                B(7,3,2),
                C(4,7);
            randomize(A);
            randomize(B);
            contract(A,{2,3,4},B,{7,3,2},C,{4,7});
            for(int i4 = 0; i4 < 4; ++i4)
            for(int i7 = 0; i7 < 7; ++i7)
                {
                Real val = 0;
                for(int i2 = 0; i2 < 2; ++i2)
                for(int i3 = 0; i3 < 3; ++i3)
                    {
                    val += A(i2,i3,i4)*B(i7,i3,i2);
                    }
                CHECK_CLOSE(C(i4,i7),val);
                }
            }

        SECTION("Case NM2")
            {
            Ten A(2,3,4,5),
                B(7,6,3,2),
                C(5,4,6,7);
            randomize(A);
            randomize(B);
            contract(A,{2,3,4,5},B,{7,6,3,2},C,{5,4,6,7});
            REQUIRE(C.extent(0) == 5);
            REQUIRE(C.extent(1) == 4);
            REQUIRE(C.extent(2) == 6);
            REQUIRE(C.extent(3) == 7);
            for(int i4 = 0; i4 < 4; ++i4)
            for(int i5 = 0; i5 < 5; ++i5)
            for(int i6 = 0; i6 < 6; ++i6)
            for(int i7 = 0; i7 < 7; ++i7)
                {
                Real val = 0;
                for(int i2 = 0; i2 < 2; ++i2)
                for(int i3 = 0; i3 < 3; ++i3)
                    {
                    val += A(i2,i3,i4,i5)*B(i7,i6,i3,i2);
                    }
                CHECK_CLOSE(C(i5,i4,i6,i7),val);
                }
            }

        SECTION("Case NM3")
            {
            Ten A(2,3,4,5),
                B(7,6,3,2),
                C(5,4,6,7);
            randomize(A);
            randomize(B);
            contract(B,{7,6,3,2},A,{2,3,4,5},C,{5,4,6,7});
            REQUIRE(C.extent(0) == 5);
            REQUIRE(C.extent(1) == 4);
            REQUIRE(C.extent(2) == 6);
            REQUIRE(C.extent(3) == 7);
            for(int i4 = 0; i4 < 4; ++i4)
            for(int i5 = 0; i5 < 5; ++i5)
            for(int i6 = 0; i6 < 6; ++i6)
            for(int i7 = 0; i7 < 7; ++i7)
                {
                Real val = 0;
                for(int i2 = 0; i2 < 2; ++i2)
                for(int i3 = 0; i3 < 3; ++i3)
                    {
                    val += A(i2,i3,i4,i5)*B(i7,i6,i3,i2);
                    }
                CHECK_CLOSE(C(i5,i4,i6,i7),val);
                }
            }

        SECTION("Case NM4")
            {
            SET_SCOPED(Global::debug1()) = true;
            Ten A(2,3,4,5,6,7),
                B(8,7,5,6,9),
                C(2,8,4,3,9);
            randomize(A);
            randomize(B);
            contract(B,{8,7,5,6,9},A,{2,3,4,5,6,7},C,{2,8,4,3,9});
            for(int i2 = 0; i2 < 2; ++i2)
            for(int i3 = 0; i3 < 3; ++i3)
            for(int i4 = 0; i4 < 4; ++i4)
            for(int i8 = 0; i8 < 8; ++i8)
            for(int i9 = 0; i9 < 9; ++i9)
                {
                Real val = 0;
                for(int i5 = 0; i5 < 5; ++i5)
                for(int i6 = 0; i6 < 6; ++i6)
                for(int i7 = 0; i7 < 7; ++i7)
                    {
                    val += A(i2,i3,i4,i5,i6,i7)*B(i8,i7,i5,i6,i9);
                    }
                CHECK_CLOSE(C(i2,i8,i4,i3,i9),val);
                }
            }


        } // Contract Reshape Non-Matrix

    SECTION("Contract Reshape Matrix")
        {
        SECTION("Case M1")
            {
            Ten A(2,3,4),
                B(7,2,3),
                C(4,7);
            randomize(A);
            randomize(B);
            contract(A,{2,3,4},B,{7,2,3},C,{4,7});
            for(int i4 = 0; i4 < 4; ++i4)
            for(int i7 = 0; i7 < 7; ++i7)
                {
                Real val = 0;
                for(int i2 = 0; i2 < 2; ++i2)
                for(int i3 = 0; i3 < 3; ++i3)
                    {
                    val += A(i2,i3,i4)*B(i7,i2,i3);
                    }
                CHECK_CLOSE(C(i4,i7),val);
                }
            }

        SECTION("Case M2")
            {
            Ten A(4,2,3),
                B(7,2,3),
                C(4,7);
            randomize(A);
            randomize(B);
            contract(A,{4,2,3},B,{7,2,3},C,{4,7});
            for(int i4 = 0; i4 < 4; ++i4)
            for(int i7 = 0; i7 < 7; ++i7)
                {
                Real val = 0;
                for(int i2 = 0; i2 < 2; ++i2)
                for(int i3 = 0; i3 < 3; ++i3)
                    {
                    val += A(i4,i2,i3)*B(i7,i2,i3);
                    }
                CHECK_CLOSE(C(i4,i7),val);
                }
            }

        SECTION("Case M3")
            {
            Ten A(4,2,3),
                B(2,3,7),
                C(4,7);
            randomize(A);
            randomize(B);
            contract(A,{4,2,3},B,{2,3,7},C,{4,7});
            for(int i4 = 0; i4 < 4; ++i4)
            for(int i7 = 0; i7 < 7; ++i7)
                {
                Real val = 0;
                for(int i2 = 0; i2 < 2; ++i2)
                for(int i3 = 0; i3 < 3; ++i3)
                    {
                    val += A(i4,i2,i3)*B(i2,i3,i7);
                    }
                CHECK_CLOSE(C(i4,i7),val);
                }
            }

        SECTION("Case M4")
            {
            Ten A(2,3,4),
                B(2,3,7),
                C(4,7);
            randomize(A);
            randomize(B);
            contract(A,{2,3,4},B,{2,3,7},C,{4,7});
            for(int i4 = 0; i4 < 4; ++i4)
            for(int i7 = 0; i7 < 7; ++i7)
                {
                Real val = 0;
                for(int i2 = 0; i2 < 2; ++i2)
                for(int i3 = 0; i3 < 3; ++i3)
                    {
                    val += A(i2,i3,i4)*B(i2,i3,i7);
                    }
                CHECK_CLOSE(C(i4,i7),val);
                }
            }

        SECTION("Regression Test 1")
            {
            Ten A(4,3,2),
                B(5,4,3,2),
                C(5);
            randomize(A);
            randomize(B);
            contract(A,{4,3,2},B,{5,4,3,2},C,{5});
            REQUIRE(C.extent(0) == 5);
            for(int i5 = 0; i5 < 5; ++i5)
                {
                Real val = 0;
                for(int i2 = 0; i2 < 2; ++i2)
                for(int i3 = 0; i3 < 3; ++i3)
                for(int i4 = 0; i4 < 4; ++i4)
                    {
                    val += A(i4,i3,i2)*B(i5,i4,i3,i2);
                    }
                CHECK_CLOSE(C(i5),val);
                }
            }

        } // Contract Reshape Matrix

    SECTION("Zero Rank Cases")
        {
        SECTION("Case 1")
            {
            Ten A(4,3,2),
                C(2,4,3);
            randomize(A);
            auto Bval = 2.;
            std::vector<Real> Bdat(1,Bval);
            Range Br;
            auto B = makeTensorRef(Bdat.data(),Br);

            contract(makeRefc(A),{4,3,2},makeRefc(B),{},makeRef(C),{2,4,3});
            for(auto i2 : count(2))
            for(auto i3 : count(3))
            for(auto i4 : count(4))
                {
                CHECK_CLOSE(C(i2,i4,i3),Bval * A(i4,i3,i2));
                }
            }

        SECTION("Case 2")
            {
            Ten A(2,3,4),
                C(2,4,3);
            randomize(A);
            auto Bval = Global::random();
            std::vector<Real> Bdat(1,Bval);
            Range Br;
            auto B = makeTensorRef(Bdat.data(),Br);

            contract<Range>(B,{},A,{2,3,4},C,{2,4,3});
            for(auto i2 : count(2))
            for(auto i3 : count(3))
            for(auto i4 : count(4))
                {
                CHECK_CLOSE(C(i2,i4,i3),Bval * A(i2,i3,i4));
                }
            }
        }

    SECTION("Contract Loop")
        {
        SECTION("Case 1: Bik Akj = Cij")
            {
            int m1 = 10,
                m2 = 20,
                m3 = 30;
            Ten A(m1,m2,4,5),
                    B(m3,m1,4,6),
                    C(m3,m2,5,6);
            randomize(A);
            randomize(B);
            contractloop<Range>(A,{1,2,4,5},B,{3,1,4,6},C,{3,2,5,6});
            for(int i2 = 0; i2 < m2; ++i2)
            for(int i3 = 0; i3 < m3; ++i3)
            for(int i5 = 0; i5 < 5; ++i5)
            for(int i6 = 0; i6 < 6; ++i6)
                {
                Real val = 0;
                for(int i1 = 0; i1 < m1; ++i1)
                for(int i4 = 0; i4 < 4; ++i4)
                    {
                    val += A(i1,i2,i4,i5)*B(i3,i1,i4,i6);
                    }
                CHECK_CLOSE(C(i3,i2,i5,i6),val);
                }
            }

        SECTION("Case 2: Bki Akj = C_ij")
            {
            int m1 = 10,
                m2 = 20,
                m3 = 30;
            Ten A(m1,m2,4,5),
                B(m1,m3,4,6),
                C(m3,m2,5,6);
            randomize(A);
            randomize(B);
            contractloop(A,{1,2,4,5},B,{1,3,4,6},C,{3,2,5,6});
            for(int i2 = 0; i2 < m2; ++i2)
            for(int i3 = 0; i3 < m3; ++i3)
            for(int i5 = 0; i5 < 5; ++i5)
            for(int i6 = 0; i6 < 6; ++i6)
                {
                Real val = 0;
                for(int i1 = 0; i1 < m1; ++i1)
                for(int i4 = 0; i4 < 4; ++i4)
                    {
                    val += A(i1,i2,i4,i5)*B(i1,i3,i4,i6);
                    }
                CHECK_CLOSE(C(i3,i2,i5,i6),val);
                }
            }

        SECTION("Case 3: Aki Bjk = Cij")
            {
            int m1 = 10,
                m2 = 20,
                m3 = 30;
            Ten A(m1,m2,4,5),
                B(m3,m1,4,6),
                C(m2,m3,5,6);
            randomize(A);
            randomize(B);
            contractloop(A,{1,2,4,5},B,{3,1,4,6},C,{2,3,5,6});
            for(int i2 = 0; i2 < m2; ++i2)
            for(int i3 = 0; i3 < m3; ++i3)
            for(int i5 = 0; i5 < 5; ++i5)
            for(int i6 = 0; i6 < 6; ++i6)
                {
                Real val = 0;
                for(int i1 = 0; i1 < m1; ++i1)
                for(int i4 = 0; i4 < 4; ++i4)
                    {
                    val += A(i1,i2,i4,i5)*B(i3,i1,i4,i6);
                    }
                CHECK_CLOSE(C(i2,i3,i5,i6),val);
                }
            }

        SECTION("Case 4: Aki Bkj = Cij")
            {
            int m1 = 10,
                m2 = 20,
                m3 = 30;
            Ten A(m1,m2,4,5),
                B(m1,m3,4,6),
                C(m2,m3,5,6);
            randomize(A);
            randomize(B);
            contractloop(A,{1,2,4,5},B,{1,3,4,6},C,{2,3,5,6});
            for(int i2 = 0; i2 < m2; ++i2)
            for(int i3 = 0; i3 < m3; ++i3)
            for(int i5 = 0; i5 < 5; ++i5)
            for(int i6 = 0; i6 < 6; ++i6)
                {
                Real val = 0;
                for(int i1 = 0; i1 < m1; ++i1)
                for(int i4 = 0; i4 < 4; ++i4)
                    {
                    val += A(i1,i2,i4,i5)*B(i1,i3,i4,i6);
                    }
                CHECK_CLOSE(C(i2,i3,i5,i6),val);
                }
            }

        SECTION("Case 5: Bik Ajk = Cij")
            {
            int m1 = 10,
                m2 = 20,
                m3 = 30;
            Ten A(m3,m1,4,5),
                B(m2,m1,4,6),
                C(m2,m3,5,6);
            randomize(A);
            randomize(B);
            contractloop(A,{3,1,4,5},B,{2,1,4,6},C,{2,3,5,6});
            for(int i2 = 0; i2 < m2; ++i2)
            for(int i3 = 0; i3 < m3; ++i3)
            for(int i5 = 0; i5 < 5; ++i5)
            for(int i6 = 0; i6 < 6; ++i6)
                {
                Real val = 0;
                for(int i1 = 0; i1 < m1; ++i1)
                for(int i4 = 0; i4 < 4; ++i4)
                    {
                    val += A(i3,i1,i4,i5)*B(i2,i1,i4,i6);
                    }
                CHECK_CLOSE(C(i2,i3,i5,i6),val);
                }
            }

        SECTION("Case 6: Aik Bjk = Cij")
            {
            int m1 = 10,
                m2 = 20,
                m3 = 30;
            Ten A(m2,m1,4,5),
                B(m3,m1,4,6),
                C(m2,m3,5,6);
            randomize(A);
            randomize(B);
            contractloop(A,{2,1,4,5},B,{3,1,4,6},C,{2,3,5,6});
            for(int i2 = 0; i2 < m2; ++i2)
            for(int i3 = 0; i3 < m3; ++i3)
            for(int i5 = 0; i5 < 5; ++i5)
            for(int i6 = 0; i6 < 6; ++i6)
                {
                Real val = 0;
                for(int i1 = 0; i1 < m1; ++i1)
                for(int i4 = 0; i4 < 4; ++i4)
                    {
                    val += A(i2,i1,i4,i5)*B(i3,i1,i4,i6);
                    }
                CHECK_CLOSE(C(i2,i3,i5,i6),val);
                }
            }

        SECTION("Case 7: Bki Ajk = Cij")
            {
            int m1 = 10,
                m2 = 20,
                m3 = 30;
            Ten A(m3,m1,4,5),
                    B(m1,m2,4,6),
                    C(m2,m3,5,6);
            randomize(A);
            randomize(B);
            contractloop(A,{3,1,4,5},B,{1,2,4,6},C,{2,3,5,6});
            for(int i2 = 0; i2 < m2; ++i2)
            for(int i3 = 0; i3 < m3; ++i3)
            for(int i5 = 0; i5 < 5; ++i5)
            for(int i6 = 0; i6 < 6; ++i6)
                {
                Real val = 0;
                for(int i1 = 0; i1 < m1; ++i1)
                for(int i4 = 0; i4 < 4; ++i4)
                    {
                    val += A(i3,i1,i4,i5)*B(i1,i2,i4,i6);
                    }
                CHECK_CLOSE(C(i2,i3,i5,i6),val);
                }

            Ten Ap(m2,m1,4,5),
                Bp(m1,m3,4,6),
                Cp(m3,m2,5,6);
            randomize(A);
            randomize(B);
            contractloop(Ap,{2,1,4,5},Bp,{1,3,4,6},Cp,{3,2,5,6});
            for(int i2 = 0; i2 < m2; ++i2)
            for(int i3 = 0; i3 < m3; ++i3)
            for(int i5 = 0; i5 < 5; ++i5)
            for(int i6 = 0; i6 < 6; ++i6)
                {
                Real val = 0;
                for(int i1 = 0; i1 < m1; ++i1)
                for(int i4 = 0; i4 < 4; ++i4)
                    {
                    val += Ap(i2,i1,i4,i5)*Bp(i1,i3,i4,i6);
                    }
                CHECK_CLOSE(Cp(i3,i2,i5,i6),val);
                }
            }

        SECTION("Case 8: Aik Bkj = Cij")
            {
            int m1 = 10,
                m2 = 20,
                m3 = 30;
            Ten A(m2,m1,4,5),
                B(m1,m3,4,6),
                C(m2,m3,5,6);
            randomize(A);
            randomize(B);
            contractloop(A,{2,1,4,5},B,{1,3,4,6},C,{2,3,5,6});
            for(int i2 = 0; i2 < m2; ++i2)
            for(int i3 = 0; i3 < m3; ++i3)
            for(int i5 = 0; i5 < 5; ++i5)
            for(int i6 = 0; i6 < 6; ++i6)
                {
                Real val = 0;
                for(int i1 = 0; i1 < m1; ++i1)
                for(int i4 = 0; i4 < 4; ++i4)
                    {
                    val += A(i2,i1,i4,i5)*B(i1,i3,i4,i6);
                    }
                CHECK_CLOSE(C(i2,i3,i5,i6),val);
                }
            }

//#define DO_TIMING

#ifdef DO_TIMING
        SECTION("Timing")
            {
            int m1 = 100,
                m2 = 200,
                m3 = 30,
                m4 = 4,
                m5 = 5,
                m6 = 6;
            Ten A(m1,m2,m4,m5),
                B(m1,m3,m4,m6),
                C(m2,m3,m5,m6);
            randomize(A);
            randomize(B);
            cpu_time cpu;
            contractloop(A,{1,2,4,5},B,{1,3,4,6},C,{2,3,5,6});
            println("Time = ",cpu.sincemark());
            for(int i2 = 0; i2 < m2; ++i2)
            for(int i3 = 0; i3 < m3; ++i3)
            for(int i5 = 0; i5 < 5; ++i5)
            for(int i6 = 0; i6 < 6; ++i6)
                {
                Real val = 0;
                for(int i1 = 0; i1 < m1; ++i1)
                for(int i4 = 0; i4 < 4; ++i4)
                    {
                    val += A(i1,i2,i4,i5)*B(i1,i3,i4,i6);
                    }
                CHECK_CLOSE(C(i2,i3,i5,i6),val);
                }
            }
#endif
    
        } // Contract Loop
    }
