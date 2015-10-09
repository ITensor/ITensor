#include "test.h"
#include "itensor/detail/algs.h"

using namespace itensor;
using namespace std;

TEST_CASE("binaryFind test")
    {

    std::vector<int> ints = {{ 1,3,6,7,9,10,12,14 }};

    SECTION("const container")
        {
        const auto& cints = ints;
        auto res = detail::binaryFind(cints,3);
        CHECK(res);
        auto assignable = std::is_assignable<decltype(*res),int>::value;
        CHECK(!assignable);
        }

    SECTION("non-const container")
        {
        auto res = detail::binaryFind(ints,3);
        CHECK(res);
        auto assignable = std::is_assignable<decltype(*res),int>::value;
        CHECK(assignable);
        }

    SECTION("Basic Functionality")
        {
        for(int i = 0; i <= 30; ++i)
            {
            auto res = detail::binaryFind(ints,i);
            //Do linear search for i to check
            bool found = false;
            for(const auto& el : ints)
                if(el == i)
                    {
                    found = true;
                    break;
                    }
            if(found) 
                {
                CHECK(res);
                CHECK(i == *res);
                }
            else CHECK(!res);
            }
        }
    }
