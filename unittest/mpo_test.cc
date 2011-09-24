#include "test.h"
#include "hams.h"
#include <boost/test/unit_test.hpp>

using std::cout;
using std::cerr;
using std::endl;

struct MPODefaults
{
    const int N;
    SpinOne::Model s1model;

    MPODefaults() : 
    N(10),
    s1model(N)
    {
    }

};

BOOST_FIXTURE_TEST_SUITE(MPOTest,MPODefaults)

BOOST_AUTO_TEST_CASE(Constructors)
{
}

BOOST_AUTO_TEST_CASE(Position)
{
    MPO H = SpinOne::Heisenberg(s1model)();
    H.position(1);
    CHECK(H.is_ortho());
    CHECK_EQUAL(H.ortho_center(),1);
}

BOOST_AUTO_TEST_SUITE_END()
