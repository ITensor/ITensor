#include "test.h"
#include "localmpo.h"
#include "model/spinhalf.h"


TEST_CASE("LocalMPOAsMPS")
    {
    static const int N = 10;
    SpinHalf shmodel(N);

    InitState shFerro(shmodel,"Up");
    InitState shNeel(shmodel);

    for(int j = 1; j <= N; ++j)
        {
        shNeel.set(j,j%2==1 ? "Up" : "Dn");
        }

    IQMPS psiNeel(shNeel),
          psiFerro(shFerro);

    LocalMPO<IQTensor> lmps(psiNeel);
    lmps.position(3,psiFerro);
    }


