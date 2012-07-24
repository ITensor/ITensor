#include "test.h"
#include <boost/test/unit_test.hpp>
#include "option.h"

using namespace std;
using namespace boost;

BOOST_AUTO_TEST_SUITE(OptionTest)

TEST(BasicUsage)
    {
    Option o1 = Quiet();
    Option o2 = Quiet();

    CHECK(o1 == o2);

    Option o3 = Quiet(false);
    CHECK(o1 == o3);
    CHECK(o2 == o3);

    CHECK(Weight() == Weight());
    CHECK(Weight(1.1) == Weight(1.1));
    CHECK(Weight(0.3) == Weight(0.2));
    CHECK(Weight(0.3) == Weight(0.3+1E-20));

    //cout << Quiet() << endl;
    //cout << Pinning(-0.42) << endl;
    }

TEST(TestOptionSet)
    {
    Option o1 = Quiet();
    Option o2 = Pinning(0.4);
    Option o3 = Auto(false);
    Option o4 = Option();

    OptionSet oset(o1,o2,o3,o4);

    CHECK(oset.defined("Quiet"));
    CHECK(oset.defined(Quiet()));
    CHECK(oset.get("Quiet").name() == "Quiet");
    CHECK(oset.boolVal("Quiet") == true);

    CHECK(oset.defined("Pinning"));
    CHECK(oset.realVal("Pinning") == 0.4);
    CHECK(oset.get("Pinning").name() == "Pinning");

    CHECK(oset.defined("Auto"));
    CHECK(oset.boolVal("Auto") == false);

    //cout << oset;
    }


BOOST_AUTO_TEST_SUITE_END()

