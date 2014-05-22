#include "test.h"
#include "model/spinone.h"
#include "hams/heisenberg.h"

using namespace itensor;
using namespace std;

TEST_CASE("MPOTest")
{
const int N = 10;
SpinOne s1model(N);

SECTION("Position")
    {
    MPO H = Heisenberg(s1model);
    H.position(1);
    CHECK(H.isOrtho());
    CHECK_EQUAL(H.orthoCenter(),1);
    }

}
