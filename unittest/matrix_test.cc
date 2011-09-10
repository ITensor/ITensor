#include "test.h"
#include "matrix.h"
#include <boost/test/unit_test.hpp>

struct MatrixDefaults
{
    MatrixDefaults() {} 
};

BOOST_FIXTURE_TEST_SUITE(MatrixTest,MatrixDefaults)

BOOST_AUTO_TEST_CASE(EigenValues)
{
    Matrix M(2,2);
    M(1,1) = 5;
    CHECK_EQUAL(M(1,1),5);

}


BOOST_AUTO_TEST_SUITE_END()

