#include "test.h"
#include "itensor/tensor/contract.h"
#include "itensor/util/iterate.h"

using namespace itensor;

TEST_CASE("Non-Contracting Product Test")
    {

    SECTION("Rank 2")
        {
        Tensor A(2,2),
               B(2,2),
               C(2,2,2);
        randomize(A);
        randomize(B);
        SECTION("Case 1")
            {
            ncprod(A,{0,1},B,{2,1},C,{0,1,2});
            for(auto i0 : range(A.extent(0)))
            for(auto i1 : range(A.extent(1)))
            for(auto i2 : range(B.extent(0)))
                {
                CHECK_CLOSE(C(i0,i1,i2), A(i0,i1)*B(i2,i1));
                }
            }
        SECTION("Case 2")
            {
            ncprod(A,{0,1},B,{1,2},C,{0,1,2});
            for(auto i0 : range(A.extent(0)))
            for(auto i1 : range(A.extent(1)))
            for(auto i2 : range(B.extent(0)))
                {
                CHECK_CLOSE(C(i0,i1,i2), A(i0,i1)*B(i1,i2));
                }
            }
        SECTION("Case 3")
            {
            ncprod(A,{1,0},B,{1,2},C,{0,1,2});
            for(auto i0 : range(A.extent(0)))
            for(auto i1 : range(A.extent(1)))
            for(auto i2 : range(B.extent(0)))
                {
                CHECK_CLOSE(C(i0,i1,i2), A(i1,i0)*B(i1,i2));
                }
            }
        SECTION("Case 4")
            {
            ncprod(A,{1,0},B,{1,2},C,{2,1,0});
            for(auto i0 : range(A.extent(0)))
            for(auto i1 : range(A.extent(1)))
            for(auto i2 : range(B.extent(0)))
                {
                CHECK_CLOSE(C(i2,i1,i0), A(i1,i0)*B(i1,i2));
                }
            }
        }

    SECTION("Ranks 3*2=4")
        {
        auto m0 = 2,
             m1 = 3,
             m2 = 4,
             m3 = 5;

        SECTION("Case 1")
            {
            Tensor A(m0,m1,m2),
                   B(m3,m1),
                   C(m0,m1,m2,m3);
            randomize(A);
            randomize(B);
            ncprod(A,{0,1,2},B,{3,1},C,{0,1,2,3});
            for(auto i0 : range(A.extent(0)))
            for(auto i1 : range(A.extent(1)))
            for(auto i2 : range(A.extent(2)))
            for(auto i3 : range(B.extent(0)))
                {
                CHECK_CLOSE(C(i0,i1,i2,i3), A(i0,i1,i2)*B(i3,i1));
                }
            }

        SECTION("Case 2")
            {
            Tensor A(m2,m1,m0),
                   B(m3,m1),
                   C(m0,m1,m2,m3);
            randomize(A);
            randomize(B);
            ncprod(A,{2,1,0},B,{3,1},C,{0,1,2,3});
            for(auto i2 : range(A.extent(0)))
            for(auto i1 : range(A.extent(1)))
            for(auto i0 : range(A.extent(2)))
            for(auto i3 : range(B.extent(0)))
                {
                CHECK_CLOSE(C(i0,i1,i2,i3), A(i2,i1,i0)*B(i3,i1));
                }
            }
        }

    }
