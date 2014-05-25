#include "test.h"
//#undef USE_CPP11 //For switching between C++11/C++98 safe bool idioms
#include "safebool.h"
#include "print.h"

using namespace itensor;
using namespace std;

class Testable : public safe_bool<Testable>
    {
    public:

    enum Initializer { Initialize };

    Testable() : initialized_(false) { }
    Testable(Initializer i) : initialized_(true) { }

    bool
    valid() const
        {
        return initialized_;
        }

    private:
    bool initialized_;
    };


TEST_CASE("SafeBool")
    {

    Testable t1,
             t2,
             t3(Testable::Initialize),
             t4(Testable::Initialize);

    SECTION("Basic")
        {
        if(t1)
            {
            //println("t1 not default constructed");
            CHECK(false);
            }
        else
            {
            //println("t1 is default constructed");
            CHECK(true);
            }

        if(!t2)
            {
            //println("t2 is default constructed");
            CHECK(true);
            }
        else
            {
            //println("t2 not default constructed");
            CHECK(false);
            }

        if(t3)
            {
            //println("t3 not default constructed");
            CHECK(true);
            }
        else
            {
            //println("t3 is default constructed");
            CHECK(false);
            }

        if(!t4)
            {
            //println("t4 is default constructed");
            CHECK(false);
            }
        else
            {
            //println("t4 not default constructed");
            CHECK(true);
            }
        }

    SECTION("Comparison")
        {
        //This should not compile:
        //if(t1 == t2) { }

        //This should not compile:
        //if(t1 < t2) { }
        }
    }

