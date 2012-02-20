#include "test.h"
#include "mps.h"
#include "model/spinhalf.h"
#include <boost/test/unit_test.hpp>

struct MPSDefaults
    {
    static const int N = 10;
    SpinHalf shmodel;

    InitState shNeel, shFerro;

    MPSDefaults() :
    shmodel(N),
    shNeel(N),
    shFerro(N)
        {
        for(int j = 1; j <= N; ++j)
            {
            shNeel(j) = (j%2==1 ? shmodel.Up(j) : shmodel.Dn(j));
            shFerro(j) = shmodel.Up(j);
            }
        }

    ~MPSDefaults() { }

    };

BOOST_FIXTURE_TEST_SUITE(MPSTest,MPSDefaults)

BOOST_AUTO_TEST_CASE(Constructors)
    {
    }

BOOST_AUTO_TEST_CASE(QNCheck)
    {
    IQMPS psiNeel(shmodel,shNeel);
    CHECK(checkQNs(psiNeel));

    CHECK_EQUAL(totalQN(psiNeel),QN(0));

    IQMPS psiFerro(shmodel,shFerro);
    CHECK(checkQNs(psiFerro));

    CHECK_EQUAL(totalQN(psiFerro),QN(10));
    }


BOOST_AUTO_TEST_SUITE_END()
