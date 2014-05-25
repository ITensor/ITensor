#include "test.h"
#include "localmpo.h"
#include "sites/spinhalf.h"

using namespace itensor;

TEST_CASE("LocalMPOAsMPS")
    {
    static const int N = 10;
    SpinHalf shsites(N);

    InitState shFerro(shsites,"Up");
    InitState shNeel(shsites);

    for(int j = 1; j <= N; ++j)
        {
        shNeel.set(j,j%2==1 ? "Up" : "Dn");
        }

    IQMPS psiNeel(shNeel),
          psiFerro(shFerro);

    LocalMPO<IQTensor> lmps(psiNeel);
    lmps.position(3,psiFerro);
    }


