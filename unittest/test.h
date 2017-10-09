#ifndef __ITENSOR_UNITTEST_TEST_H
#define __ITENSOR_UNITTEST_TEST_H

#include "catch.hpp"

#define CHECK_EQUAL(X,Y) CHECK((X) == (Y))
#define CHECK_CLOSE(X,Y) CHECK(std::norm((X)-(Y)) < 1E-10)
//#define CHECK_CLOSE(X,Y)\
//    if(std::norm(X) < std::norm(Y)) CHECK(std::norm((X)-(Y)) < 1E-10*std::norm(Y));\
//    else CHECK(std::norm((X)-(Y)) < 1E-10*std::norm(X))

#define CHECK_REQUAL(X,Y) CHECK_CLOSE(X,Y)
#define CHECK_DIFF(X,Y,T) CHECK(std::norm((X)-(Y)) < T)
#define DISABLE if(false)

#endif
