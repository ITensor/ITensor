#include "test.h"
#include "itensor/mps/mpo.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/mps/sites/spinone.h"
#include "itensor/util/print_macro.h"

using namespace itensor;
using namespace std;

TEST_CASE("MPOTest")
{

SECTION("Orthogonalize")
    {
    auto N = 10;
    auto m = 4;
    auto sites = SpinHalf(10);
    auto W = MPO(sites);

    //Make a random MPS of bond dim. m
    auto links = vector<Index>(N+1);
    for(auto n : range1(N))
        {
        links.at(n) = Index(nameint("l",n),m);
        }
    W.Aref(1) = randomTensor(links.at(1),sites(1),prime(sites(1)));
    for(auto n : range1(2,N-1))
        {
        W.Aref(n) = randomTensor(links.at(n-1),sites(n),prime(sites(n)),links.at(n));
        }
    W.Aref(N) = randomTensor(links.at(N-1),sites(N),prime(sites(N)));

    //Normalize W
    auto n2 = overlap(W,W);
    W.Aref(1) /= sqrt(n2);

    auto oW = W;

    W.orthogonalize();

    CHECK_CLOSE(overlap(oW,W),1.0);

    for(int n = N; n > 1; --n)
        {
        auto li = commonIndex(W.A(n),W.A(n-1),Link);
        auto rho = W.A(n) * dag(prime(W.A(n),li));
        auto id = ITensor(li,prime(li));
        for(auto l : range1(li.m()))
            {
            id.set(li(l),prime(li)(l),1.0);
            }
        CHECK(norm(rho-id) < 1E-10);
        }
    }

}
