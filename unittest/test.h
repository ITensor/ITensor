#ifndef __ITENSOR_UNITTEST_TEST_H
#define __ITENSOR_UNITTEST_TEST_H

#include "catch.hpp"

#define CHECK_EQUAL(X,Y) REQUIRE((X) == (Y))
#define CHECK_CLOSE(X,Y,T) REQUIRE(fabs((X)-(Y)) < T)

#endif
