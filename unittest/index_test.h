#ifndef __INDEX_TEST_H
#define __INDEX_TEST_H
#include "test.h"
#include "index.h"

BOOST_AUTO_TEST_SUITE(IndexTest)

BOOST_AUTO_TEST_CASE(Null)
{
    Index i1;
    CHECK(i1.is_null());
    CHECK_EQUAL(1,i1.m());

    Index i2("i2");
    CHECK(i2.is_not_null());
}

BOOST_AUTO_TEST_SUITE_END()

#endif
