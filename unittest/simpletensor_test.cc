#include "test.h"
#include "print.h"
#include "simpletensor.h"

using namespace itensor;

using STensor = simpletensor<Real>;

TEST_CASE("Simpletensor")
    {
    SECTION("Default Constructor")
        {
        STensor s;
        CHECK(s.r() == 0);
        CHECK(s.size() == 0);
        CHECK(!s);
        }

    SECTION("Variadic Constructor")
        {
            {
            long m0 = 12;
            STensor s(m0);
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

            {
            long m0 = 3,
                 m1 = 7;
            STensor s(m0,m1);
            CHECK(s.r() == 2);
            CHECK(s.n(0) == m0);
            CHECK(s.n(1) == m1);
            CHECK(s.size() == m0*m1);
            CHECK(s.ownstorage());
            }

            {
            long m0 = 3,
                 m1 = 7,
                 m2 = 2;
            STensor s(m0,m1,m2);
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
        STensor::index ind = {5,2};
        STensor t(&(v.front()),ind);

        CHECK(t);
        CHECK(!t.ownstorage());
        CHECK(t.n(0) == 5);
        CHECK(t.n(1) == 2);
        CHECK(t.size() == 0);

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
        STensor::index ind = {5,2};
        STensor::storage v= {11,21,31,41,51,
                             12,22,32,42,52};

        STensor t(std::move(ind),std::move(v));

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
