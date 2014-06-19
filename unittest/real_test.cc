#include "test.h"

#include "global.h"
#include "real.h"

using namespace itensor;

TEST_CASE("Constructors")
    {
    LogNumber l1;

    REQUIRE( std::isnan(l1.logNum()) );
    CHECK_EQUAL(l1.sign(),1);
    REQUIRE( std::isnan(l1.real()) );

    LogNumber l2(1);

    CHECK_CLOSE(l2.logNum(),0,LogNumber_Accuracy);
    CHECK_EQUAL(l2.sign(),1);
    CHECK_CLOSE(l2.real(),1,LogNumber_Accuracy);

    LogNumber l3(-1);

    CHECK_CLOSE(l3.logNum(),0,LogNumber_Accuracy);
    CHECK_EQUAL(l3.sign(),-1);
    CHECK_CLOSE(l3.real(),-1,LogNumber_Accuracy);

    LogNumber l4(0);

    CHECK_CLOSE(l4.logNum(),0,LogNumber_Accuracy);
    CHECK_EQUAL(l4.sign(),0);
    CHECK_CLOSE(l4.real(),0,LogNumber_Accuracy);
    REQUIRE(l4.isRealZero());

    const Real Big= 1E50;
    const int BigExp = 50;

    LogNumber l5(Big);

    CHECK_CLOSE(l5.logNum(),BigExp*log(10),LogNumber_Accuracy);
    CHECK_EQUAL(l5.sign(),1);
    //CHECK_CLOSE(l5.real(),Big,LogNumber_Accuracy);

    LogNumber l6(-Big);

    CHECK_CLOSE(l6.logNum(),BigExp*log(10),LogNumber_Accuracy);
    CHECK_EQUAL(l6.sign(),-1);
    //CHECK_CLOSE(l6.real(),-Big,LogNumber_Accuracy);

    Real r = Global::random();
    LogNumber l7(r);

    CHECK_CLOSE(l7.logNum(),log(fabs(r)),LogNumber_Accuracy);
    CHECK_EQUAL(l7.sign(),(r > 0 ? 1 : -1));
    //CHECK_CLOSE(l7.real(),r,LogNumber_Accuracy);
    }

TEST_CASE("Operators")
    {
    Real a = Global::random(), b = Global::random();

    const LogNumber la(a), lb(b);

    CHECK_CLOSE((la*lb).real(),a*b,LogNumber_Accuracy);
    CHECK_CLOSE((la/lb).real(),a/b,LogNumber_Accuracy);
    CHECK_CLOSE((lb/la).real(),b/a,LogNumber_Accuracy);

    LogNumber l1(a);
    l1 *= -1;

    CHECK_CLOSE((-l1).real(),a,LogNumber_Accuracy);

    LogNumber l2(a);
    l2 *= b;

    CHECK_CLOSE(l2.real(),a*b,LogNumber_Accuracy);

    LogNumber l3(a), l4(b);
    l3 *= l4;

    CHECK_CLOSE(l3.real(),a*b,LogNumber_Accuracy);

    LogNumber l5(a),l6(b);
    l5 *= -l6;

    CHECK_CLOSE(l5.real(),-a*b,LogNumber_Accuracy);

    LogNumber l7(a);
    l7 /= b;

    CHECK_CLOSE(l7.real(),a/b,LogNumber_Accuracy);

    LogNumber l8(a),l9(b);
    l8 /= l9;

    CHECK_CLOSE(l8.real(),a/b,LogNumber_Accuracy);
    }

TEST_CASE("Comparison")
    {
    Real a = Global::random(), 
         b = Global::random();

    const LogNumber la(a),lb(b);

    CHECK_EQUAL(la <  lb,a <  b);
    CHECK_EQUAL(la >  lb,a >  b);
    CHECK_EQUAL(la <= lb,a <= b);
    CHECK_EQUAL(la >= lb,a >= b);
    CHECK_EQUAL(la == lb,a == b);

    const LogNumber dla(a+a*0.1*LogNumber_Accuracy);
    REQUIRE(la.approxEquals(dla));
    REQUIRE(la != dla);

    const LogNumber mdla(a-a*0.1*LogNumber_Accuracy);
    REQUIRE(la.approxEquals(mdla));
    REQUIRE(la != mdla);

    REQUIRE((la.magnitudeLessThan(lb) || lb.magnitudeLessThan(la)));

    const LogNumber p(1E-10),q(1E-12),r(-1E-12);
    REQUIRE(q.magnitudeLessThan(p));
    REQUIRE(r.magnitudeLessThan(p));
    REQUIRE(!p.magnitudeLessThan(p));

    LogNumber zero(0);

    REQUIRE(zero.approxEquals(zero));

    zero *= -1;

    REQUIRE(zero.approxEquals(zero));

    LogNumber one(1);

    REQUIRE(zero < one);
    }

