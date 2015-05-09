#include "test.h"
#include "svdalgs.h"

using namespace itensor;
using namespace std;

TEST_CASE("SVDAlgsTest")
{

SECTION("Transpose SVD")
    {
    Index a("a",3),
          b("b",2);

    ITensor A(a,b);
    auto el = 1./sqrt(2);

    A.set(el,a(1),b(1));
    A.set(-el,a(2),b(1));
    A.set(el,a(3),b(1));
    A.set(el,a(1),b(2));
    A.set(el,a(2),b(2));
    A.set(el,a(3),b(2));

    ITensor U(b),D,V;
    svd(A,U,D,V,{"Truncate",false});

    CHECK(norm(A-U*D*V) < 1E-10);
    }
}
