#include "test.h"
#include <boost/test/unit_test.hpp>
#include "option.h"

BOOST_AUTO_TEST_SUITE(OptionTest)

TEST(BasicUsage)
    {
    Option o1 = Quiet();
    Option o2 = Quiet();

    CHECK(o1 == o2);

    Option o3 = Quiet(false);
    CHECK(o1 != o3);
    CHECK(o2 != o3);

    CHECK(Weight() == Weight());
    CHECK(Weight(1.1) == Weight(1.1));
    CHECK(Weight(0.3) != Weight(0.2));
    CHECK(Weight(0.3) == Weight(0.3+1E-20));
    }

TEST(TestOptionSet)
    {
    Option o1 = Quiet();
    Option o2 = DoPinning(0.4);
    Option o3 = Auto(false);
    Option o4 = Option();

    OptionSet oset(o1,o2,o3,o4);

    CHECK(oset.includes("Quiet"));
    CHECK(oset.includes(Quiet()));

    CHECK(oset.get("Quiet").name() == "Quiet");
    CHECK(oset.boolVal("Quiet") == true);

    CHECK(oset.includes("DoPinning"));
    CHECK(oset.includes(DoPinning(0.4)));
    CHECK(!oset.includes(DoPinning(0.3)));

    CHECK(oset.get("DoPinning").name() == "DoPinning");
    CHECK(oset.realVal("DoPinning") == 0.4);

    CHECK(oset.includes("Auto"));
    CHECK(oset.includes(Auto(false)));
    CHECK(!oset.includes(Auto(true)));

    CHECK(oset.get("Auto").name() == "Auto");
    CHECK(oset.boolVal("Auto") == false);
    }


BOOST_AUTO_TEST_SUITE_END()

