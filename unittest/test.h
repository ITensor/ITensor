#ifndef __ITENSOR_UNITTEST_TEST_H
#define __ITENSOR_UNITTEST_TEST_H

#include "catch.hpp"

#define CHECK_EQUAL(X,Y) REQUIRE((X) == (Y))
#define CHECK_CLOSE(X,Y) REQUIRE(std::norm((X)-(Y)) < 1E-11)
#define CHECK_REQUAL(X,Y) CHECK_CLOSE(X,Y)
#define CHECK_DIFF(X,Y,T) REQUIRE(std::norm((X)-(Y)) < T)
#define DISABLE if(false)

#endif
