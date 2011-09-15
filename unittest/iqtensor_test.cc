#include "test.h"
#include "iqtensor.h"
#include <boost/test/unit_test.hpp>

struct IQTensorDefaults
{
    IQTensorDefaults()
    {
    }

    ~IQTensorDefaults() { }

};

BOOST_FIXTURE_TEST_SUITE(IQTensorTest,IQTensorDefaults)

BOOST_AUTO_TEST_CASE(Null)
{
    IQTensor t1;

    CHECK(t1.is_null());
}

BOOST_AUTO_TEST_CASE(Constructors)
{
}


BOOST_AUTO_TEST_SUITE_END()
