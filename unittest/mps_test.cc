#include "test.h"
#include "mps.h"
#include "model/spinhalf.h"
#include "model/spinless.h"

using namespace itensor;

TEST_CASE("MPSTest")
{

static const int N = 10;
SpinHalf shmodel(N);

InitState shFerro(shmodel,"Up");
InitState shNeel(shmodel);

for(int j = 1; j <= N; ++j)
    {
    shNeel.set(j,j%2==1 ? "Up" : "Dn");
    }

SECTION("Constructors")
    {
    }

SECTION("QNCheck")
    {
    IQMPS psiNeel(shNeel);
    CHECK(checkQNs(psiNeel));

    CHECK_EQUAL(totalQN(psiNeel),QN(0));

    IQMPS psiFerro(shFerro);
    CHECK(checkQNs(psiFerro));

    CHECK_EQUAL(totalQN(psiFerro),QN(10));
    }

SECTION("MPSAddition")
    {
    Spinless model(10);

    InitState i1(model,"Emp"),
              i2(model,"Emp");

    i1.set(1,"Occ");
    i2.set(2,"Occ");

    //"Valence bond" between sites 1 and 2
    MPS psi = ISqrt2*sum(MPS(i1),MPS(i2));

    CHECK_CLOSE(psi.norm(),1,1E-5);

    IQMPS iqpsi = ISqrt2*sum(IQMPS(i1),IQMPS(i2));

    CHECK_EQUAL(totalQN(iqpsi),QN(0,1));
    }

SECTION("PositionTest")
    {
    Spinless model(10);

    InitState init(model,"Emp");
    init.set(2,"Occ");
    init.set(4,"Occ");
    init.set(6,"Occ");

    IQMPS psi(init);
    psi.Anc(1) *= Complex_i;

    psi.position(1,"Cutoff=1E-8");
    CHECK_EQUAL(findCenter(psi),1);

    psi.position(4,"Cutoff=1E-8");
    CHECK_EQUAL(findCenter(psi),4);
    }


}
