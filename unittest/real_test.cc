#include "test.h"
#include "real.h"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(LogNumberTest)

BOOST_AUTO_TEST_CASE(Constructors)
{
    LogNumber l1;

    CHECK_EQUAL(l1.logNum(),0);
    CHECK_EQUAL(l1.sign(),1);
    CHECK_EQUAL(Real(l1),1);

    LogNumber l2(1);

    CHECK_CLOSE(l2.logNum(),0,LogNumber_Accuracy);
    CHECK_EQUAL(l2.sign(),1);
    CHECK_CLOSE(Real(l2),1,LogNumber_Accuracy);

    LogNumber l3(-1);

    CHECK_CLOSE(l3.logNum(),0,LogNumber_Accuracy);
    CHECK_EQUAL(l3.sign(),-1);
    CHECK_CLOSE(Real(l3),-1,LogNumber_Accuracy);

    LogNumber l4(0);

    CHECK_CLOSE(l4.logNum(),0,LogNumber_Accuracy);
    CHECK_EQUAL(l4.sign(),0);
    CHECK_CLOSE(Real(l4),0,LogNumber_Accuracy);
    CHECK(l4.isRealZero());

    const Real Big= 1E50;
    const int BigExp = 50;

    LogNumber l5(Big);

    CHECK_CLOSE(l5.logNum(),BigExp*log(10),LogNumber_Accuracy);
    CHECK_EQUAL(l5.sign(),1);
    CHECK_CLOSE(Real(l5),Big,LogNumber_Accuracy);

    LogNumber l6(-Big);

    CHECK_CLOSE(l6.logNum(),BigExp*log(10),LogNumber_Accuracy);
    CHECK_EQUAL(l6.sign(),-1);
    CHECK_CLOSE(Real(l6),-Big,LogNumber_Accuracy);

    Real r = ran1();
    LogNumber l7(r);

    CHECK_CLOSE(l7.logNum(),log(fabs(r)),LogNumber_Accuracy);
    CHECK_EQUAL(l7.sign(),(r > 0 ? 1 : -1));
    CHECK_CLOSE(Real(l7),r,LogNumber_Accuracy);
}

BOOST_AUTO_TEST_CASE(Operators)
{
    Real a = ran1(), b = ran1();

    const LogNumber la(a), lb(b);

    CHECK_CLOSE(Real(la*lb),a*b,LogNumber_Accuracy);
    CHECK_CLOSE(Real(la/lb),a/b,LogNumber_Accuracy);
    CHECK_CLOSE(Real(lb/la),b/a,LogNumber_Accuracy);

    LogNumber l1(a);
    l1 *= -1;

    CHECK_CLOSE(Real(l1)*-1,a,LogNumber_Accuracy);

    LogNumber l2(a);
    l2 *= b;

    CHECK_CLOSE(Real(l2),a*b,LogNumber_Accuracy);

    LogNumber l3(a), l4(b);
    l3 *= l4;

    CHECK_CLOSE(Real(l3),a*b,LogNumber_Accuracy);

    LogNumber l5(a),l6(b);
    l5 *= -l6;

    CHECK_CLOSE(Real(l5),-a*b,LogNumber_Accuracy);

    LogNumber l7(a);
    l7 /= b;

    CHECK_CLOSE(Real(l7),a/b,LogNumber_Accuracy);

    LogNumber l8(a),l9(b);
    l8 /= l9;

    CHECK_CLOSE(Real(l8),a/b,LogNumber_Accuracy);
}

BOOST_AUTO_TEST_CASE(Comparison)
{
    Real a = ran1(), b = ran1();

    const LogNumber la(a),lb(b);

    CHECK_EQUAL(la <  lb,a <  b);
    CHECK_EQUAL(la >  lb,a >  b);
    CHECK_EQUAL(la <= lb,a <= b);
    CHECK_EQUAL(la >= lb,a >= b);
    CHECK_EQUAL(la == lb,a == b);

    const LogNumber dla(a+0.1*LogNumber_Accuracy);
    CHECK(la.approxEquals(dla));
    CHECK(la != dla);

    const LogNumber mdla(a-0.1*LogNumber_Accuracy);
    CHECK(la.approxEquals(mdla));
    CHECK(la != mdla);

    CHECK(la.magnitudeLessThan(lb) || lb.magnitudeLessThan(la));

    const LogNumber p(1E-10),q(1E-12),r(-1E-12);
    CHECK(q.magnitudeLessThan(p));
    CHECK(r.magnitudeLessThan(p));
    CHECK(!p.magnitudeLessThan(p));

    LogNumber zero(0);

    CHECK(zero.approxEquals(0));

    zero *= -1;

    CHECK(zero.approxEquals(0));
}

BOOST_AUTO_TEST_CASE(ReadWrite)
{
    Real a = ran1();
    LogNumber la(a); 

    writeToFile(".read_write/LogNumber",la);

    LogNumber laR;
    readFromFile(".read_write/LogNumber",laR);

    CHECK(laR == la);
}

BOOST_AUTO_TEST_SUITE_END()

