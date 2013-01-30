#include "test.h"
#include <boost/test/unit_test.hpp>
#include "global.h"
#include "option.h"

using namespace std;
using namespace boost;

BOOST_AUTO_TEST_SUITE(OptTest)

TEST(BasicUsage)
    {
    Opt o1 = Quiet();
    Opt o2 = Quiet();

    CHECK(o1 == o2);

    Opt o3 = Quiet(false);
    CHECK(o1 == o3);
    CHECK(o2 == o3);

    CHECK(Weight() == Weight());
    CHECK(Weight(1.1) == Weight(1.1));
    CHECK(Weight(0.3) == Weight(0.2));
    CHECK(Weight(0.3) == Weight(0.3+1E-20));

    //cout << Quiet() << endl;
    //cout << Pinning(-0.42) << endl;
    }

TEST(TestOptSet)
    {
    OptSet& gopts = OptSet::GlobalOpts();

    Opt o1 = Quiet();
    Opt o2 = Pinning(0.4);
    Opt o3 = Auto(false);
    Opt o4 = Opt();

    gopts.add(o1,o2);

    CHECK(gopts.getBool("Quiet") == true);

    //cout << "Global opts: " << endl;
    //cout << gopts << endl;

    OptSet opts1(Quiet(false));

    //cout << "opts1: " << endl;
    //cout << opts1 << endl;

    CHECK(opts1.defined("Quiet"));
    CHECK(opts1.defined(Quiet()));
    CHECK(opts1.get("Quiet").name() == "Quiet");
    CHECK(opts1.getBool("Quiet") == false);

    OptSet opts2(opts1);
    opts2.add(o3,o4);

    //cout << "opts2: " << endl;
    //cout << opts2 << endl;

    CHECK(opts2.defined("Pinning"));
    CHECK(opts2.getReal("Pinning") == 0.4);
    CHECK(opts2.get("Pinning").name() == "Pinning");

    CHECK(opts2.defined("Auto"));
    CHECK(opts2.getBool("Auto") == false);
    }

TEST(Operator)
    {
    OptSet oset1 = Quiet() & Auto();
    CHECK(oset1.defined("Quiet"));
    CHECK(oset1.defined("Auto"));

    OptSet oset2 = Quiet() & Auto() & Pinning(1);
    CHECK(oset1.defined("Quiet"));
    CHECK(oset1.defined("Auto"));
    CHECK(oset1.defined("Pinning"));
    }


BOOST_AUTO_TEST_SUITE_END()

