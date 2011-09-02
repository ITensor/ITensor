#define THIS_IS_MAIN
#include "dmrg.h"

#define BOOST_TEST_MODULE Hello
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( my_test )
{
    Index i1;
    BOOST_CHECK(i1.is_null());
    i1 = Index("i1");
    BOOST_CHECK(i1.is_not_null());

    cerr << format("i1 %s null\n")%(i1.is_null() ? "is" : "is not");

}
