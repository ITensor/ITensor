#include "test.h"
#include "print.h"
#include "simpletensor.h"
#include "types.h"

using namespace itensor;

using RTensor = tensor<Real>;
using RTRef = tensorref<Real,Range>;

TEST_CASE("tensor")
    {
    SECTION("Default Constructor")
        {
        RTensor s;
        CHECK(s.r() == 0);
        CHECK(s.size() == 0);
        CHECK(!s);
        }

    SECTION("Variadic Constructor")
        {
        SECTION("Case 1")
            {
            long m0 = 12;
            RTensor s(m0);
            CHECK(s);
            CHECK(s.r() == 1);
            CHECK(s.n(0) == m0);
            CHECK(s.size() == m0);
            CHECK_REQUAL(s(0),0);
            CHECK_REQUAL(s(1),0);
            //...
            CHECK_REQUAL(s(m0-1),0);
            CHECK(s.ownstorage());
            }

        SECTION("Case 2")
            {
            long m0 = 3,
                 m1 = 7;
            RTensor s(m0,m1);
            CHECK(s.r() == 2);
            CHECK(s.n(0) == m0);
            CHECK(s.n(1) == m1);
            CHECK(s.size() == m0*m1);
            CHECK(s.ownstorage());
            }

        SECTION("Case 3")
            {
            long m0 = 3,
                 m1 = 7,
                 m2 = 2;
            RTensor s(m0,m1,m2);
            CHECK(s.r() == 3);
            CHECK(s.n(0) == m0);
            CHECK(s.n(1) == m1);
            CHECK(s.n(2) == m2);
            CHECK(s.size() == m0*m1*m2);
            CHECK(s.ownstorage());

            Real total = 0;
            for(auto x : s.store()) total += x;
            CHECK(total == 0);
            }
        }

    SECTION("Non-owning")
        {
        std::vector<Real> v({11,21,31,41,51,
                             12,22,32,42,52});
        Range ind{5,2};
        RTRef t(v.data(),ind);

        CHECK(t);
        CHECK(!t.ownstorage());
        CHECK(t.n(0) == 5);
        CHECK(t.n(1) == 2);
        CHECK(t.size() == 10);

        CHECK_REQUAL(t(0,0),11);
        CHECK_REQUAL(t(0,1),12);
        CHECK_REQUAL(t(2,0),31);
        CHECK_REQUAL(t(2,1),32);
        CHECK_REQUAL(t(3,0),41);
        CHECK_REQUAL(t(3,1),42);
        CHECK_REQUAL(t(4,0),51);
        CHECK_REQUAL(t(4,1),52);
        }

    SECTION("Construct from Data")
        {
        auto ind = Range{5,2};
        auto v = RTensor::storage_type({11,21,31,41,51,
                                        12,22,32,42,52});

        SECTION("Copy in Data")
            {
            RTensor t(ind,v.begin(),v.end());

            CHECK(t.ownstorage());
            CHECK(t.n(0) == 5);
            CHECK(t.n(1) == 2);
            CHECK(t.size() == 5*2);

            CHECK_REQUAL(t(0,0),11);
            CHECK_REQUAL(t(0,1),12);
            CHECK_REQUAL(t(2,0),31);
            CHECK_REQUAL(t(2,1),32);
            CHECK_REQUAL(t(3,0),41);
            CHECK_REQUAL(t(3,1),42);
            CHECK_REQUAL(t(4,0),51);
            CHECK_REQUAL(t(4,1),52);
            }

        SECTION("Move in Data")
            {
            RTensor t(std::move(ind),std::move(v));

            //Check that ind and v got moved
            CHECK(ind.empty());
            CHECK(v.empty());

            CHECK(t.ownstorage());
            CHECK(t.n(0) == 5);
            CHECK(t.n(1) == 2);
            CHECK(t.size() == 5*2);

            CHECK_REQUAL(t(0,0),11);
            CHECK_REQUAL(t(0,1),12);
            CHECK_REQUAL(t(2,0),31);
            CHECK_REQUAL(t(2,1),32);
            CHECK_REQUAL(t(3,0),41);
            CHECK_REQUAL(t(3,1),42);
            CHECK_REQUAL(t(4,0),51);
            CHECK_REQUAL(t(4,1),52);
            }
        }

    }
