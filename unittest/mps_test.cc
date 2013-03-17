#include "test.h"
#include "mps.h"
#include "model/spinhalf.h"
#include "model/spinless.h"
#include <boost/test/unit_test.hpp>

struct MPSDefaults
    {
    static const int N = 10;
    SpinHalf shmodel;

    InitState shNeel, 
              shFerro;

    MPSDefaults() :
    shmodel(N),
    shNeel(shmodel),
    shFerro(shmodel,&SpinHalf::Up)
        {
        for(int j = 1; j <= N; ++j)
            {
            shNeel.set(j,j%2==1 ? &SpinHalf::Up : &SpinHalf::Dn);
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
    IQMPS psiNeel(shmodel,shNeel);
    CHECK(checkQNs(psiNeel));

    CHECK_EQUAL(totalQN(psiNeel),QN(0));

    IQMPS psiFerro(shmodel,shFerro);
    CHECK(checkQNs(psiFerro));

    CHECK_EQUAL(totalQN(psiFerro),QN(10));
    }

TEST(MPSAddition)
    {
    Spinless model(10);

    InitState i1(model,&Spinless::Emp),
              i2(model,&Spinless::Emp);

    i1.set(1,&Spinless::Occ);
    i2.set(2,&Spinless::Occ);

    //"Valence bond" between sites 1 and 2
    MPS psi = ISqrt2*(MPS(model,i1) + MPS(model,i2));

    CHECK_CLOSE(psi.norm(),1,1E-5);

    IQMPS iqpsi = ISqrt2*(IQMPS(model,i1) + IQMPS(model,i2));

    CHECK_EQUAL(totalQN(iqpsi),QN(0,1));
    }

TEST(PositionTest)
    {
    Spinless model(10);

    InitState init(model,&Spinless::Emp);
    init.set(2,&Spinless::Occ);
    init.set(4,&Spinless::Occ);
    init.set(6,&Spinless::Occ);

    IQMPS psi(model,init);
    psi.cutoff(1E-8);
    psi.noise(1E-8);
    psi.Anc(1) *= IQComplex_i();

    psi.position(1);
    CHECK_EQUAL(findCenter(psi),1);

    psi.position(4);
    CHECK_EQUAL(findCenter(psi),4);
    }


BOOST_AUTO_TEST_SUITE_END()
