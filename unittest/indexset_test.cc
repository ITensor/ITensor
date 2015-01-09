#include "test.h"
#include "indexset.h"

using namespace itensor;

TEST_CASE("IndexSetTest")
    {
    Index i1("i1",1),
          i2("i2",2),
          i3("i3",3),
          i4("i4",4),
          i5("i5",5),
          i6("i6",6),
          i7("i7",7),
          i8("i8",8),
          i9("i9",9),
          i10("i10",10),
          j1("j1",1),
          j2("j2",2),
          j3("j3",3),
          j4("j4",4),
          j5("j5",5),
          j6("j6",6),
          j7("j7",7),
          j8("j8",8),
          j9("j9",9),
          j10("j10",10);

    SECTION("Constructors")
        {
        SECTION("One")
            {
            auto is = IndexSet(i4);
            CHECK(is.r() == 1);
            CHECK(is.rn() == 1);
            CHECK(area(is) == 4);
            CHECK(is[0] == i4);
            CHECK(is[0] == is.index(1));

            is = IndexSet(i1);
            CHECK(is.r() == 1);
            CHECK(is.rn() == 0);
            CHECK(area(is) == 1);
            CHECK(is[0] == i1);
            CHECK(is[0] == is.index(1));
            }

        SECTION("Two")
            {
            auto is1 = IndexSet(i4,i3);
            CHECK(is1.r() == 2);
            CHECK(is1.rn() == 2);
            CHECK(area(is1) == 4*3);
            CHECK(is1[0] == i4);
            CHECK(is1[1] == i3);

            auto is2 = IndexSet(i3,i4);
            CHECK(is2.r() == 2);
            CHECK(is2.rn() == 2);
            CHECK(area(is2) == 4*3);
            //IndexSet automatically sorts the indices
            CHECK(is2[0] == i4);
            CHECK(is2[1] == i3);

            auto isA = IndexSet(i3,j3);
            auto isB = IndexSet(j3,i3);
            //IndexSet sorts equal-m indices
            //by their id #
            CHECK(isA[0] == isA[0]);
            CHECK(isB[1] == isB[1]);
            }

        SECTION("THREE")
            {
            auto is = IndexSet(i3,i1,i4);
            CHECK(is.r() == 3);
            CHECK(is.rn() == 2);
            CHECK(area(is) == 4*3);
            CHECK(is[0] == i4);
            CHECK(is[1] == i3);
            CHECK(is[2] == i1);
            }

        SECTION("TEN")
            {
            auto is = IndexSet(i3,i1,i4,i2,i5,i6,i7,i8,i10,i9);
            CHECK(is.r() == 10);
            CHECK(is.rn() == 9);
            CHECK(is[0] == i10);
            CHECK(is[1] == i9);
            CHECK(is[2] == i8);
            CHECK(is[8] == i2);
            CHECK(is[9] == i1);
            }
        }

    SECTION("PrimeLevelMethods")
        {
        SECTION("PrimeIndex")
            {
            }

        SECTION("NoPrimeIndex")
            {
            }

        SECTION("NoPrimeType")
            {
            }

        SECTION("AddIndex")
            {
            }
        }
    }
