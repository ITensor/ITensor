#ifndef __ITENSOR_TEST_H
#define __ITENSOR_TEST_H
#include "test.h"

struct ITensorDefaults
{
    const Index s1,s2,s3,s4,
          s1P,s2P,s3P,s4P,
          l1,l2,l3,l4,l5,l6,l7,l8,
          a1,a2,a3,a4;

    ITensorDefaults() :
    s1(Index("s1",2,Site)),
    s2(Index("s2",2,Site)),
    s3(Index("s3",2,Site)),
    s4(Index("s4",2,Site)),
    s1P(s1.primed()),
    s2P(s2.primed()),
    s3P(s3.primed()),
    s4P(s4.primed()),
    l1(Index("l1",2)),
    l2(Index("l2",2)),
    l3(Index("l3",2)),
    l4(Index("l4",2)),
    l5(Index("l5",2)),
    l6(Index("l6",2)),
    l7(Index("l7",2)),
    l8(Index("l8",2)),
    a1(Index("a1")),
    a2(Index("a2")),
    a3(Index("a3")),
    a4(Index("a4"))
    {
    }

    ~ITensorDefaults() { }

};

BOOST_FIXTURE_TEST_SUITE(ITensorTest,ITensorDefaults)

BOOST_AUTO_TEST_CASE(Null)
{
    ITensor t1;

    BOOST_CHECK(t1.is_null());

    ITensor t2(s1);

    BOOST_CHECK(t2.is_not_null());
}

BOOST_AUTO_TEST_CASE(Constructors)
{
    ITensor t1(l1);

    BOOST_CHECK_EQUAL(t1.r(),1);
    BOOST_CHECK_EQUAL(t1.r_n(),1);
    BOOST_CHECK_EQUAL(t1.r_1(),0);
    BOOST_CHECK(t1.hasindex(l1));
    BOOST_CHECK_CLOSE(t1.norm(),0,1E-10);

    ITensor t2(l1,l2);

    BOOST_CHECK_EQUAL(t2.r(),2);
    BOOST_CHECK_EQUAL(t2.r_n(),2);
    BOOST_CHECK_EQUAL(t2.r_1(),0);
    BOOST_CHECK(t2.hasindex(l1));
    BOOST_CHECK(t2.hasindex(l2));
    BOOST_CHECK_CLOSE(t2.norm(),0,1E-10);

    ITensor t3(l1,l2,l3);

    BOOST_CHECK_EQUAL(t3.r(),3);
    BOOST_CHECK_EQUAL(t3.r_n(),3);
    BOOST_CHECK_EQUAL(t3.r_1(),0);
    BOOST_CHECK(t3.hasindex(l1));
    BOOST_CHECK(t3.hasindex(l2));
    BOOST_CHECK(t3.hasindex(l3));
    BOOST_CHECK_CLOSE(t3.norm(),0,1E-10);

    ITensor t4(a1,l1);

    BOOST_CHECK_EQUAL(t4.r(),2);
    BOOST_CHECK_EQUAL(t4.r_n(),1);
    BOOST_CHECK_EQUAL(t4.r_1(),1);
    BOOST_CHECK(t4.hasindex(a1));
    BOOST_CHECK(t4.hasindex(l1));
    BOOST_CHECK_CLOSE(t4.norm(),0,1E-10);

    ITensor t5(l1,a1,l2);

    BOOST_CHECK_EQUAL(t5.r(),3);
    BOOST_CHECK_EQUAL(t5.r_n(),2);
    BOOST_CHECK_EQUAL(t5.r_1(),1);
    BOOST_CHECK(t5.hasindex(a1));
    BOOST_CHECK(t5.hasindex(l1));
    BOOST_CHECK(t5.hasindex(l2));
    BOOST_CHECK_CLOSE(t5.norm(),0,1E-10);

    ITensor t6(l1,a1,l2,a2);

    BOOST_CHECK_EQUAL(t6.r(),4);
    BOOST_CHECK_EQUAL(t6.r_n(),2);
    BOOST_CHECK_EQUAL(t6.r_1(),2);
    BOOST_CHECK(t6.hasindex(l1));
    BOOST_CHECK(t6.hasindex(a1));
    BOOST_CHECK(t6.hasindex(l2));
    BOOST_CHECK(t6.hasindex(a2));
    BOOST_CHECK_CLOSE(t6.norm(),0,1E-10);

    Real a = ran1();
    ITensor t7(l1,l2,a);

    BOOST_CHECK_EQUAL(t7.r(),2);
    BOOST_CHECK_EQUAL(t7.r_n(),2);
    BOOST_CHECK_EQUAL(t7.r_1(),0);
    BOOST_CHECK(t7.hasindex(l1));
    BOOST_CHECK(t7.hasindex(l2));
    BOOST_CHECK_CLOSE(t7.val2(1,1),a,1E-10);
    BOOST_CHECK_CLOSE(t7.val2(1,2),0,1E-10);
    BOOST_CHECK_CLOSE(t7.val2(2,1),0,1E-10);
    BOOST_CHECK_CLOSE(t7.val2(2,2),a,1E-10);
    BOOST_CHECK_CLOSE(t7.norm(),sqrt(min(l1.m(),l2.m()))*fabs(a),1E-10);

    Matrix M(l1.m(),l2.m()); 
    M(1,1) = ran1(); M(1,2) = ran1();
    M(2,1) = ran1(); M(2,2) = ran1();
    ITensor t8(l1,l2,M);

    BOOST_CHECK_EQUAL(t8.r(),2);
    BOOST_CHECK_EQUAL(t8.r_n(),2);
    BOOST_CHECK_EQUAL(t8.r_1(),0);
    BOOST_CHECK(t8.hasindex(l1));
    BOOST_CHECK(t8.hasindex(l2));
    BOOST_CHECK_CLOSE(t8.val2(1,1),M(1,1),1E-10);
    BOOST_CHECK_CLOSE(t8.val2(1,2),M(1,2),1E-10);
    BOOST_CHECK_CLOSE(t8.val2(2,1),M(2,1),1E-10);
    BOOST_CHECK_CLOSE(t8.val2(2,2),M(2,2),1E-10);
    BOOST_CHECK_CLOSE(t8.sumels(),M.TreatAsVector().sumels(),1E-10);
    BOOST_CHECK_CLOSE(t8.norm(),Norm(M.TreatAsVector()),1E-10);

    Real b = ran1();
    ITensor t9(b);

    BOOST_CHECK_CLOSE(t9.sumels(),b,1E-10);
    BOOST_CHECK_CLOSE(t9.norm(),fabs(b),1E-10);

    Index linkind("linkind",10);
    Vector V(linkind.m()); V.Randomize();
    ITensor t10(linkind,V);

    BOOST_CHECK_EQUAL(t10.r(),1);
    BOOST_CHECK_EQUAL(t10.r_n(),1);
    BOOST_CHECK_EQUAL(t10.r_1(),0);
    BOOST_CHECK(t10.hasindex(linkind));
    BOOST_CHECK_CLOSE(t10.sumels(),V.sumels(),1E-10);
    BOOST_CHECK_CLOSE(t10.norm(),Norm(V),1E-10);
}

BOOST_AUTO_TEST_CASE(IndexValConstructors)
{
    ITensor t1(l1(2));

    BOOST_CHECK_EQUAL(t1.r(),1);
    BOOST_CHECK_EQUAL(t1.r_n(),1);
    BOOST_CHECK_EQUAL(t1.r_1(),0);
    BOOST_CHECK(t1.hasindex(l1));
    BOOST_CHECK_CLOSE(t1.val1(1),0,1E-10);
    BOOST_CHECK_CLOSE(t1.val1(2),1,1E-10);
    BOOST_CHECK_CLOSE(t1.sumels(),1,1E-10);
    BOOST_CHECK_CLOSE(t1.norm(),1,1E-10);

    ITensor t2(l1(2),l2(1));

    BOOST_CHECK_EQUAL(t2.r(),2);
    BOOST_CHECK_EQUAL(t2.r_n(),2);
    BOOST_CHECK_EQUAL(t2.r_1(),0);
    BOOST_CHECK(t2.hasindex(l1));
    BOOST_CHECK(t2.hasindex(l2));
    BOOST_CHECK_CLOSE(t2.val2(1,1),0,1E-10);
    BOOST_CHECK_CLOSE(t2.val2(2,1),1,1E-10);
    BOOST_CHECK_CLOSE(t2.val2(1,2),0,1E-10);
    BOOST_CHECK_CLOSE(t2.val2(2,2),0,1E-10);
    BOOST_CHECK_CLOSE(t2.sumels(),1,1E-10);
    BOOST_CHECK_CLOSE(t2.norm(),1,1E-10);
}

BOOST_AUTO_TEST_CASE(IndexVectorConstructors)
{
    std::vector<Index> indices(4);
    indices[0] = a2;
    indices[1] = l3;
    indices[2] = l1;
    indices[3] = a4;

    ITensor t1(indices);

    BOOST_CHECK_EQUAL(t1.r(),4);
    BOOST_CHECK_EQUAL(t1.r_n(),2);
    BOOST_CHECK_EQUAL(t1.r_1(),2);
    BOOST_CHECK(t1.hasindex(a2));
    BOOST_CHECK(t1.hasindex(l3));
    BOOST_CHECK(t1.hasindex(l1));
    BOOST_CHECK(t1.hasindex(a4));
    BOOST_CHECK_CLOSE(t1.norm(),0,1E-10);

    Vector V(l1.m()*l3.m());
    V.Randomize();

    ITensor t2(indices,V);

    BOOST_CHECK_EQUAL(t2.r(),4);
    BOOST_CHECK_EQUAL(t2.r_n(),2);
    BOOST_CHECK_EQUAL(t2.r_1(),2);
    BOOST_CHECK(t2.hasindex(a2));
    BOOST_CHECK(t2.hasindex(l3));
    BOOST_CHECK(t2.hasindex(l1));
    BOOST_CHECK(t2.hasindex(a4));
    BOOST_CHECK_CLOSE(t2.norm(),Norm(V),1E-10);
    BOOST_CHECK_CLOSE(t2.sumels(),V.sumels(),1E-10);
}

BOOST_AUTO_TEST_CASE(Copy)
{
    std::vector<Index> indices(4);
    indices[0] = a2;
    indices[1] = l3;
    indices[2] = l1;
    indices[3] = a4;

    Vector V(l1.m()*l3.m());
    V.Randomize();

    ITensor t1(indices,V);

    BOOST_CHECK_EQUAL(t1.r(),4);
    BOOST_CHECK_EQUAL(t1.r_n(),2);
    BOOST_CHECK_EQUAL(t1.r_1(),2);
    BOOST_CHECK(t1.hasindex(a2));
    BOOST_CHECK(t1.hasindex(l3));
    BOOST_CHECK(t1.hasindex(l1));
    BOOST_CHECK(t1.hasindex(a4));
    BOOST_CHECK_CLOSE(t1.norm(),Norm(V),1E-10);
    BOOST_CHECK_CLOSE(t1.sumels(),V.sumels(),1E-10);

    //Use copy constructor
    ITensor t2(t1);
    t1 = ITensor(); //destroy t1

    BOOST_CHECK_EQUAL(t2.r(),4);
    BOOST_CHECK_EQUAL(t2.r_n(),2);
    BOOST_CHECK_EQUAL(t2.r_1(),2);
    BOOST_CHECK(t2.hasindex(a2));
    BOOST_CHECK(t2.hasindex(l3));
    BOOST_CHECK(t2.hasindex(l1));
    BOOST_CHECK(t2.hasindex(a4));
    BOOST_CHECK_CLOSE(t2.norm(),Norm(V),1E-10);
    BOOST_CHECK_CLOSE(t2.sumels(),V.sumels(),1E-10);

    //Use operator=
    ITensor t3 = t2;
    t2 = ITensor(); //destroy t2

    BOOST_CHECK_EQUAL(t3.r(),4);
    BOOST_CHECK_EQUAL(t3.r_n(),2);
    BOOST_CHECK_EQUAL(t3.r_1(),2);
    BOOST_CHECK(t3.hasindex(a2));
    BOOST_CHECK(t3.hasindex(l3));
    BOOST_CHECK(t3.hasindex(l1));
    BOOST_CHECK(t3.hasindex(a4));
    BOOST_CHECK_CLOSE(t3.norm(),Norm(V),1E-10);
    BOOST_CHECK_CLOSE(t3.sumels(),V.sumels(),1E-10);
}

BOOST_AUTO_TEST_SUITE_END()

#endif
