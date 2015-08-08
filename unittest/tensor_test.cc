#include "test.h"
#include "itensor/util/print.h"
#include "itensor/tensor/ten.h"

using namespace itensor;

//template<size_t n>
//struct choice : choice<n+1>
//    {
//    constexpr choice(){}
//    };
//
//template<>
//struct choice<10>
//    {
//    constexpr choice(){}
//    };
//
//struct select_overload : choice<1> { };
//
//auto
//checkCompilesImpl(choice<1>,TensorRefc ct)
//    -> std::conditional_t<std::is_void<decltype(ct(4,0) = 5.1)>::value,bool,bool>
//    {
//    return true;
//    }
//
//bool
//checkCompilesImpl(choice<2>,TensorRefc ct)
//    {
//    return false;
//    }
//
//template<typename... VArgs>
//bool
//checkCompiles(VArgs&&... vargs) { return checkCompilesImpl(select_overload{},std::forward<VArgs>(vargs)...); }

TEST_CASE("Tensor and TensorRef")
{

SECTION("TensorRef")
    {
    auto v = std::vector<Real>{11,21,31,41,51,
                               12,22,32,42,52};
    auto ind = Range{5,2};

    SECTION("Owns range")
        {
        auto t = makeTenRef(v.data(),std::move(ind));
        CHECK(t.ownRange() == true);
        }

    SECTION("Doesn't own range")
        {
        auto t = makeTenRef(v.data(),ind);
        CHECK(t.ownRange() == false);
        }

    SECTION("Constructor Basics")
        {
        auto t = makeTenRef(v.data(),ind);

        CHECK(t);
        CHECK(t.r() == 2);
        CHECK(t.extent(0) == 5);
        CHECK(t.extent(1) == 2);
        CHECK(t.size() == 10);
        }

    SECTION("Non-const Element Access")
        {
        auto t = makeTenRef(v.data(),ind);

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
        auto ct = makeTenRef(cv,ind);
        CHECK_CLOSE(ct(4,0),51);
        static_assert(std::is_same<decltype(ct(4,0)),const Real&>::value,
                      "Type of ct(4,0) is not const Real&");
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
            CHECK_CLOSE(s(m0-1),0);
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
            auto t = Tensor(std::move(ind),std::move(v));

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
    }
}
