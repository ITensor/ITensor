#include "test.h"

#include "itensor/global.h"
#include "itensor/real.h"

using namespace itensor;

TEST_CASE("Constructors")
    {
    LogNum l1;

    REQUIRE( std::isnan(l1.logNum()) );
    CHECK_EQUAL(l1.sign(),1);
    REQUIRE( std::isnan(l1.real()) );

    LogNum l2(1);

    CHECK_DIFF(l2.logNum(),0,LogNum_Accuracy);
    CHECK_EQUAL(l2.sign(),1);
    CHECK_DIFF(l2.real(),1,LogNum_Accuracy);

    LogNum l3(-1);

    CHECK_DIFF(l3.logNum(),0,LogNum_Accuracy);
    CHECK_EQUAL(l3.sign(),-1);
    CHECK_DIFF(l3.real(),-1,LogNum_Accuracy);

    LogNum l4(0);

    CHECK_DIFF(l4.logNum(),0,LogNum_Accuracy);
    CHECK_EQUAL(l4.sign(),0);
    CHECK_DIFF(l4.real(),0,LogNum_Accuracy);
    REQUIRE(l4.isRealZero());

    const Real Big= 1E50;
    const int BigExp = 50;

    LogNum l5(Big);

    CHECK_DIFF(l5.logNum(),BigExp*log(10),LogNum_Accuracy);
    CHECK_EQUAL(l5.sign(),1);
    //CHECK_DIFF(l5.real(),Big,LogNum_Accuracy);

    LogNum l6(-Big);

    CHECK_DIFF(l6.logNum(),BigExp*log(10),LogNum_Accuracy);
    CHECK_EQUAL(l6.sign(),-1);
    //CHECK_DIFF(l6.real(),-Big,LogNum_Accuracy);

    Real r = Global::random();
    LogNum l7(r);

    CHECK_DIFF(l7.logNum(),log(fabs(r)),LogNum_Accuracy);
    CHECK_EQUAL(l7.sign(),(r > 0 ? 1 : -1));
    //CHECK_DIFF(l7.real(),r,LogNum_Accuracy);
    }

TEST_CASE("Operators")
    {
    Real a = Global::random(), b = Global::random();

    const LogNum la(a), lb(b);

    CHECK_DIFF((la*lb).real(),a*b,LogNum_Accuracy);
    CHECK_DIFF((la/lb).real(),a/b,LogNum_Accuracy);
    CHECK_DIFF((lb/la).real(),b/a,LogNum_Accuracy);

    LogNum l1(a);
    l1 *= -1;

    CHECK_DIFF((-l1).real(),a,LogNum_Accuracy);

    LogNum l2(a);
    l2 *= b;

    CHECK_DIFF(l2.real(),a*b,LogNum_Accuracy);

    LogNum l3(a), l4(b);
    l3 *= l4;

    CHECK_DIFF(l3.real(),a*b,LogNum_Accuracy);

    LogNum l5(a),l6(b);
    l5 *= -l6;

    CHECK_DIFF(l5.real(),-a*b,LogNum_Accuracy);

    LogNum l7(a);
    l7 /= b;

    CHECK_DIFF(l7.real(),a/b,LogNum_Accuracy);

    LogNum l8(a),l9(b);
    l8 /= l9;

    CHECK_DIFF(l8.real(),a/b,LogNum_Accuracy);
    }

TEST_CASE("Comparison")
    {
    Real a = Global::random(), 
         b = Global::random();

    const LogNum la(a),lb(b);

    CHECK_EQUAL(la <  lb,a <  b);
    CHECK_EQUAL(la >  lb,a >  b);
    CHECK_EQUAL(la <= lb,a <= b);
    CHECK_EQUAL(la >= lb,a >= b);
    CHECK_EQUAL(la == lb,a == b);

    const LogNum dla(a+a*0.1*LogNum_Accuracy);
    REQUIRE(la.approxEquals(dla));
    REQUIRE(la != dla);

    const LogNum mdla(a-a*0.1*LogNum_Accuracy);
    REQUIRE(la.approxEquals(mdla));
    REQUIRE(la != mdla);

    REQUIRE((la.magnitudeLessThan(lb) || lb.magnitudeLessThan(la)));

    const LogNum p(1E-10),q(1E-12),r(-1E-12);
    REQUIRE(q.magnitudeLessThan(p));
    REQUIRE(r.magnitudeLessThan(p));
    REQUIRE(!p.magnitudeLessThan(p));

    LogNum zero(0);

    REQUIRE(zero.approxEquals(zero));

    zero *= -1;

    REQUIRE(zero.approxEquals(zero));

    LogNum one(1);

    REQUIRE(zero < one);
    }

