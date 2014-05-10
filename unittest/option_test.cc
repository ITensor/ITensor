#include "test.h"
#include <boost/test/unit_test.hpp>
#include "global.h"
#include "option.h"

using namespace std;
using namespace boost;

OptSet
string_test_function(const OptSet& opts)
    {
    return opts;
    }

BOOST_AUTO_TEST_SUITE(OptTest)

TEST(BasicUsage)
    {
    Opt o1("Quiet");
    Opt o2("Quiet");

    CHECK(o1 == o2);

    Opt o3("Quiet",false);
    CHECK(o1 == o3);
    CHECK(o2 == o3);

    Opt o4("Weight",0.3);
    CHECK_EQUAL(o4.name(),"Weight");
    CHECK_EQUAL(o4.type(),Opt::Numeric);
    CHECK_CLOSE(o4.realVal(),0.3,1E-5);

    Opt o5("UseSVD",true);
    CHECK_EQUAL(o5.name(),"UseSVD");
    CHECK_EQUAL(o5.type(),Opt::Boolean);
    CHECK_EQUAL(o5.boolVal(),true);
    CHECK_THROW(o5.intVal(),ITError);

    Opt o6("UseSVD",false);
    CHECK_EQUAL(o6.name(),"UseSVD");
    CHECK_EQUAL(o6.type(),Opt::Boolean);
    CHECK_EQUAL(o6.boolVal(),false);
    CHECK_THROW(o6.realVal(),ITError);

    Opt o7("Name","index");
    CHECK_EQUAL(o7.name(),"Name");
    CHECK_EQUAL(o7.type(),Opt::String);
    CHECK_EQUAL(o7.stringVal(),"index");
    CHECK_THROW(o7.boolVal(),ITError);
    }

TEST(TestOptSet)
    {
    OptSet& gopts = OptSet::GlobalOpts();

    Opt o1("Quiet");
    Opt o2 = Opt("Pinning",0.4);
    Opt o3("Auto",false);
    Opt o4 = Opt();

    gopts.add(o1,o2);

    CHECK(gopts.getBool("Quiet") == true);

    //cout << "Global opts: " << endl;
    //cout << gopts << endl;

    OptSet opts1(Opt("Quiet",false));

    //cout << "opts1: " << endl;
    //cout << opts1 << endl;

    CHECK(opts1.defined("Quiet"));
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

    OptSet opts3 = opts1 & opts2;

    CHECK(opts3.defined("Quiet"));
    CHECK(opts3.defined("Auto"));
    }

TEST(Operator)
    {
    OptSet oset1; 
    oset1 &= "Quiet";
    oset1 &= "Auto";
    CHECK(oset1.defined("Quiet"));
    CHECK(oset1.defined("Auto"));

    OptSet oset2 = Opt("Pinning",1) & "Quiet" & "Auto";
    CHECK(oset1.defined("Quiet"));
    CHECK(oset1.defined("Auto"));
    CHECK(oset1.defined("Pinning"));
    }

TEST(StringConstructor)
    {
    OptSet o1("Quiet,Auto,Pinning=-0.5");
    CHECK(o1.defined("Quiet"));
    CHECK(o1.defined("Auto"));
    CHECK(o1.defined("Pinning"));
    CHECK(o1.getReal("Pinning") == -0.5);

    //string_test_function is defined at the top of this file
    OptSet o2 = string_test_function("WriteM=500,UseSVD=false,");
    CHECK(o2.getInt("WriteM") == 500);
    CHECK(o2.getBool("UseSVD") == false);
    
    OptSet o3 = string_test_function(Opt("Quiet") & "WriteM=500,UseSVD=false,");
    CHECK(o3.getBool("Quiet") == true);
    CHECK(o3.getInt("WriteM") == 500);
    CHECK(o3.getBool("UseSVD") == false);

    OptSet o4("Name=new,Cutoff=5E-12,Val=-4.235235");
    CHECK(o4.getString("Name") == "new");
    CHECK(fabs(o4.getReal("Cutoff")-5E-12) < 1E-14);
    CHECK(fabs(o4.getReal("Val")-(-4.235235)) < 1E-12);
    }


BOOST_AUTO_TEST_SUITE_END()

