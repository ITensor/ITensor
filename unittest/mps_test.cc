#include "test.h"
#include "mps.h"
#include "model/spinhalf.h"
#include "model/spinless.h"
#include <boost/test/unit_test.hpp>

using namespace itensor;

struct MPSDefaults
    {
    static const int N = 10;
    SpinHalf shmodel;

    InitState shNeel, 
              shFerro;

    MPSDefaults() :
    shmodel(N),
    shNeel(shmodel),
    shFerro(shmodel,"Up")
        {
        for(int j = 1; j <= N; ++j)
            {
            shNeel.set(j,j%2==1 ? "Up" : "Dn");
            }
        }

    ~MPSDefaults() { }

    };

BOOST_FIXTURE_TEST_SUITE(MPSTest,MPSDefaults)

TEST(Constructors)
    {
    }

TEST(QNCheck)
    {
    IQMPS psiNeel(shNeel);
    CHECK(checkQNs(psiNeel));

    CHECK_EQUAL(totalQN(psiNeel),QN(0));

    IQMPS psiFerro(shFerro);
    CHECK(checkQNs(psiFerro));

    CHECK_EQUAL(totalQN(psiFerro),QN(10));
    }

TEST(MPSAddition)
    {
    Spinless model(10);

    InitState i1(model,"Emp"),
              i2(model,"Emp");

    i1.set(1,"Occ");
    i2.set(2,"Occ");

    //"Valence bond" between sites 1 and 2
    MPS psi = ISqrt2*(MPS(i1) + MPS(i2));

    CHECK_CLOSE(psi.norm(),1,1E-5);

    IQMPS iqpsi = ISqrt2*(IQMPS(i1) + IQMPS(i2));

    CHECK_EQUAL(totalQN(iqpsi),QN(0,1));
    }

TEST(PositionTest)
    {
    Spinless model(10);

    InitState init(model,"Emp");
    init.set(2,"Occ");
    init.set(4,"Occ");
    init.set(6,"Occ");

    IQMPS psi(init);
    psi.cutoff(1E-8);
    psi.noise(1E-8);
    psi.Anc(1) *= Complex_i;

    psi.position(1);
    CHECK_EQUAL(findCenter(psi),1);

    psi.position(4);
    CHECK_EQUAL(findCenter(psi),4);
    }


BOOST_AUTO_TEST_SUITE_END()
