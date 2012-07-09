#include "test.h"
#include "localmpo.h"
#include "model/spinhalf.h"
#include <boost/test/unit_test.hpp>

struct LocalMPODefaults
    {
    static const int N = 10;
    SpinHalf shmodel;

    InitState shNeel, shFerro;

    LocalMPODefaults() :
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

    ~LocalMPODefaults() { }

    };

BOOST_FIXTURE_TEST_SUITE(LocalMPOTest,LocalMPODefaults)

BOOST_AUTO_TEST_CASE(LocalMPOAsMPS)
    {
    IQMPS psiNeel(shmodel,shNeel),
          psiFerro(shmodel,shFerro);

    LocalMPO<IQTensor> lmps(psiNeel);
    lmps.position(3,psiFerro);
    }


BOOST_AUTO_TEST_SUITE_END()
