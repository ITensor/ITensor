#include "test.h"
#include "print.h"
#include "tensor/ten.h"
#include "types.h"

using namespace itensor;

TEST_CASE("Tensor and TensorRef")
{

SECTION("TensorRef")
    {
    auto v = std::vector<Real>{11,21,31,41,51,
                               12,22,32,42,52};
    Range ind{5,2};
    auto t = makeTenRef(v.data(),ind);

    CHECK(t);
    CHECK(t.r() == 2);
    CHECK(t.dim(0) == 5);
    CHECK(t.dim(1) == 2);
    CHECK(t.size() == 10);

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

    const auto* cv = v.data();
    auto ct = makeTenRef(cv,ind);
    CHECK_CLOSE(ct(4,0),5.1);
    }

SECTION("Tensor")
    {
    SECTION("Default Constructor")
        {
        Ten s;
        CHECK(s.r() == 0);
        CHECK(s.size() == 0);
        CHECK(!s);
        }

    SECTION("Tensor Constructors")
        {
        SECTION("Case 1")
            {
            long m0 = 12;
            Ten s(m0);
            CHECK(s);
            CHECK(s.r() == 1);
            CHECK(s.dim(0) == m0);
            CHECK(s.size() == m0);
            CHECK_CLOSE(s(0),0);
            CHECK_CLOSE(s(1),0);
            //...
            CHECK_CLOSE(s(m0-1),0);
            }

        SECTION("Case 2")
            {
            long m0 = 3,
                 m1 = 7;
            Ten s(m0,m1);
            CHECK(s.r() == 2);
            CHECK(s.dim(0) == m0);
            CHECK(s.dim(1) == m1);
            CHECK(s.size() == m0*m1);
            }

        SECTION("Case 3")
            {
            long m0 = 3,
                 m1 = 7,
                 m2 = 2;
            Ten s(m0,m1,m2);
            CHECK(s.r() == 3);
            CHECK(s.dim(0) == m0);
            CHECK(s.dim(1) == m1);
            CHECK(s.dim(2) == m2);
            CHECK(s.size() == m0*m1*m2);

            Real total = 0;
            for(auto x : s) total += x;
            CHECK(total == 0);
            }
        }


    SECTION("Construct from Data")
        {
        auto ind = Range{5,2};
        auto v = Ten::storage_type({11,21,31,41,51,
                                    12,22,32,42,52});
        SECTION("Move in Data")
            {
            Ten t(std::move(ind),std::move(v));

            //Check that ind and v got moved
            CHECK(ind.empty());
            CHECK(v.empty());

            CHECK(t.dim(0) == 5);
            CHECK(t.dim(1) == 2);
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
    }
}
