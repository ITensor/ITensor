#include "test.h"
#include "itensor/util/print.h"
#include "itensor/util/count.h"
#include "itensor/detail/algs.h"
#include "itensor/matrix/permutation.h"
#include "itensor/matrix/sliceten.h"

using namespace itensor;


TEST_CASE("Tensor and TensorRef")
{

SECTION("Range")
    {
    SECTION("RangeBuilder")
        {
        auto r = 3;
        auto B = RangeBuilder(r);
        B.nextIndex(4);
        B.nextIndex(3);
        B.nextIndex(2);

        SECTION("Basics")
            {
            CHECK(B);
            auto R = B.build();
            CHECK((not B));
            CHECK(R.r() == r);
            }
        SECTION("Auto Strides")
            {
            auto R = B.build();
            CHECK(R.stride(0) == 1);
            CHECK(R.stride(1) == 4);
            CHECK(R.stride(2) == 12);
            CHECK(isContiguous(R));
            }

        SECTION("Manual Strides")
            {
            B.setStride(0,6);
            B.setStride(1,2);
            B.setStride(2,1);
            auto R = B.build();
            CHECK(R.stride(0) == 6);
            CHECK(R.stride(1) == 2);
            CHECK(R.stride(2) == 1);
            CHECK(isContiguous(R));
            }
        }

    } // Range

SECTION("TensorRef")
    {
    auto v = std::vector<Real>{11,21,31,41,51,
                               12,22,32,42,52};
    auto ind = Range{5,2};

    SECTION("Owns range")
        {
        auto t = makeTenRef(v.data(),v.size(),std::move(ind));
        CHECK(t.ownRange() == true);
        }

    SECTION("Doesn't own range")
        {
        auto t = makeTenRef(v.data(),v.size(),&ind);
        CHECK(t.ownRange() == false);
        }

    SECTION("Constructor Basics")
        {
        auto t = makeTenRef(v.data(),v.size(),&ind);

        CHECK(t);
        CHECK(t.r() == 2);
        CHECK(t.extent(0) == 5);
        CHECK(t.extent(1) == 2);
        CHECK(t.size() == 10);
        }

    SECTION("Non-const Element Access")
        {
        auto t = makeTenRef(v.data(),v.size(),&ind);

        CHECK_CLOSE(t(0,0),11);
        CHECK_CLOSE(t(0,1),12);
        CHECK_CLOSE(t(2,0),31);
        CHECK_CLOSE(t(2,1),32);
        CHECK_CLOSE(t(3,0),41);
        CHECK_CLOSE(t(3,1),42);
        CHECK_CLOSE(t(4,0),51);
        CHECK_CLOSE(t(4,1),52);

        t(4,0) = 5.1;
        CHECK_CLOSE(t(4,0),5.1);
        }

    SECTION("Const Element Access")
        {
        const auto* cv = v.data();
        auto ct = makeTenRef(cv,v.size(),&ind);
        CHECK_CLOSE(ct(4,0),51);
        static_assert(std::is_same<decltype(ct(5,1)),const Real&>::value,
                      "Type of ct(5,1) is not const Real&");
        }
    }

SECTION("Tensor")
    {
    SECTION("Default Constructor")
        {
        Tensor s;
        CHECK(s.r() == 0);
        CHECK(s.size() == 0);
        CHECK(!s);
        }

    SECTION("Tensor Constructors")
        {
        SECTION("Case 1")
            {
            long m0 = 12;
            auto s = Tensor(m0);
            CHECK(s);
            CHECK(s.r() == 1);
            CHECK(s.extent(0) == m0);
            CHECK(s.size() == m0);
            CHECK_CLOSE(s(0),0);
            CHECK_CLOSE(s(1),0);
            //...
            CHECK_CLOSE(s(m0),0);
            }

        SECTION("Case 2")
            {
            long m0 = 3,
                 m1 = 7;
            auto s = Tensor(m0,m1);
            CHECK(s.r() == 2);
            CHECK(s.extent(0) == m0);
            CHECK(s.extent(1) == m1);
            CHECK(s.size() == m0*m1);
            }

        SECTION("Case 3")
            {
            long m0 = 3,
                 m1 = 7,
                 m2 = 2;
            auto s = Tensor(m0,m1,m2);
            CHECK(s.r() == 3);
            CHECK(s.extent(0) == m0);
            CHECK(s.extent(1) == m1);
            CHECK(s.extent(2) == m2);
            CHECK(s.size() == m0*m1*m2);

            Real total = 0;
            for(auto x : s) total += x;
            CHECK(total == 0);
            }
        }


    SECTION("Construct from Data")
        {
        auto ind = Range{5,2};
        auto v = Tensor::storage_type({11,21,31,41,51,
                                       12,22,32,42,52});
        SECTION("Move in Data")
            {
            auto t = Tensor(std::move(v),std::move(ind));

            //Check that ind and v got moved
            CHECK(ind.empty());
            CHECK(v.empty());

            CHECK(t.extent(0) == 5);
            CHECK(t.extent(1) == 2);
            CHECK(t.size() == 5*2);

            CHECK_CLOSE(t(0,0),11);
            CHECK_CLOSE(t(0,1),12);
            CHECK_CLOSE(t(2,0),31);
            CHECK_CLOSE(t(2,1),32);
            CHECK_CLOSE(t(3,0),41);
            CHECK_CLOSE(t(3,1),42);
            CHECK_CLOSE(t(4,0),51);
            CHECK_CLOSE(t(4,1),52);
            }
        }
    } // Tensor

SECTION("Slicing")
    {
    SECTION("Permute")
        {
        auto T = Tensor(4,2,3);
        for(auto& el : T) el = detail::quickran();

        SECTION("Case 1")
            {
            auto PT = permute(T,Label{0,2,1});
            for(auto i0 : count(T.extent(0)))
            for(auto i1 : count(T.extent(1)))
            for(auto i2 : count(T.extent(2)))
                {
                CHECK_CLOSE(PT(i0,i1,i2), T(i0,i2,i1));
                }
            }

        SECTION("Case 1.5")
            {
            auto P = Permutation(3);
            P.setFromTo(0,0);
            P.setFromTo(1,2);
            P.setFromTo(2,1);
            //println("P=",P);
            //print("P="); for(auto el : P) print(el," "); println();
            auto PT = permute(T,P);
            for(auto i0 : count(T.extent(0)))
            for(auto i1 : count(T.extent(1)))
            for(auto i2 : count(T.extent(2)))
                {
                CHECK_CLOSE(PT(i0,i1,i2), T(i0,i2,i1));
                }
            }

        SECTION("Case 2")
            {
            auto PT = permute(T,Label{2,0,1});
            for(auto& i : PT.range())
                {
                CHECK_CLOSE(PT(i), T(i[2],i[0],i[1]));
                }
            }

        SECTION("Case 3")
            {
            auto PT = permute(T,Label{2,1,0});
            for(auto& i : PT.range())
                {
                CHECK_CLOSE(PT(i), T(i[2],i[1],i[0]));
                }
            }

        SECTION("Case 4")
            {
            auto PT = permute(T,Label{1,2,0});
            for(auto& i : PT.range())
                {
                CHECK_CLOSE(PT(i), T(i[1],i[2],i[0]));
                }
            }

        }

    SECTION("Sub Tensor")
        {
        auto T = Tensor(7,3,8,6);
        for(auto& el : T) el = detail::quickran();

        SECTION("Case 1")
            {
            Label start = {0,0,0,0},
                  stop  = {7,3,8,6};
            auto S = subTensor(T,start,stop);
            for(auto& i : S.range())
                {
                CHECK_CLOSE(S(i), T(start[0]+i[0],
                                    start[1]+i[1],
                                    start[2]+i[2],
                                    start[3]+i[3]));
                CHECK_CLOSE(S(i),T(i));
                }
            }

        SECTION("Case 2")
            {
            Label start = {1,0,1,0},
                  stop  = {7,3,8,6};
            auto S = subTensor(T,start,stop);
            for(auto& i : S.range())
                {
                CHECK_CLOSE(S(i), T(start[0]+i[0],
                                    start[1]+i[1],
                                    start[2]+i[2],
                                    start[3]+i[3]));
                }
            }

        SECTION("Case 3")
            {
            Label start = {1,0,1,0},
                  stop  = {7,3,8,6};
            auto S = subTensor(T,start,stop);
            for(auto& i : S.range())
                {
                CHECK_CLOSE(S(i), T(start[0]+i[0],
                                    start[1]+i[1],
                                    start[2]+i[2],
                                    start[3]+i[3]));
                }
            }

        SECTION("Case 4")
            {
            Label start = {2,1,1,2},
                  stop  = {4,3,4,5};
            auto S = subTensor(T,start,stop);
            for(auto& i : S.range())
                {
                CHECK_CLOSE(S(i), T(start[0]+i[0],
                                    start[1]+i[1],
                                    start[2]+i[2],
                                    start[3]+i[3]));
                }
            }

        SECTION("Case 5")
            {
            Label start = {5,2,7,5},
                  stop  = {7,3,8,6};
            auto S = subTensor(T,start,stop);
            for(auto& i : S.range())
                {
                CHECK_CLOSE(S(i), T(start[0]+i[0],
                                    start[1]+i[1],
                                    start[2]+i[2],
                                    start[3]+i[3]));
                }
            }

        SECTION("Case 6")
            {
            Label start = {5,1,5,4},
                  stop  = {6,2,7,5};
            auto S = subTensor(T,start,stop);
            for(auto& i : S.range())
                {
                CHECK_CLOSE(S(i), T(start[0]+i[0],
                                    start[1]+i[1],
                                    start[2]+i[2],
                                    start[3]+i[3]));
                }
            }

        }

    } // Slicing
}
