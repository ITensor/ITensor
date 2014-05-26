#include "test.h"
#include "sites/spinone.h"
#include "hams/Heisenberg.h"

using namespace itensor;
using namespace std;

TEST_CASE("MPOTest")
{
const int N = 10;
SpinOne s1sites(N);

SECTION("Position")
    {
    MPO H = Heisenberg(s1sites);
    H.position(1);
    CHECK(H.isOrtho());
    CHECK_EQUAL(H.orthoCenter(),1);
    }

}
