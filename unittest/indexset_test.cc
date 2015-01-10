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
            CHECK(is2[0] == i3);
            CHECK(is2[1] == i4);
            }

        SECTION("THREE")
            {
            auto is = IndexSet(i3,i1,i4);
            CHECK(is.r() == 3);
            CHECK(is.rn() == 2);
            CHECK(area(is) == 4*3);
            CHECK(is[0] == i3);
            CHECK(is[1] == i4);
            CHECK(is[2] == i1);
            }

        SECTION("TEN")
            {
            auto is = IndexSet(i3,i1,i4,i2,i5,i6,i7,i8,i10,i9);
            CHECK(is.r() == 10);
            CHECK(is.rn() == 9);
            CHECK(is[0] == i3);
            CHECK(is[1] == i4);
            CHECK(is[2] == i2);
            CHECK(is[3] == i5);
            CHECK(is[4] == i6);
            CHECK(is[5] == i7);
            CHECK(is[6] == i8);
            CHECK(is[7] == i10);
            CHECK(is[8] == i9);
            CHECK(is[9] == i1);
            }

        SECTION("m=1 Indices")
            {
            Index a1("a1",1),
                  a2("a2",1),
                  a3("a3",1);
            IndexSet is(j8,a2,j9,a1,a3);
            CHECK(is[0] == j8);
            CHECK(is[1] == j9);
            CHECK(is[2] == a2);
            CHECK(is[3] == a1);
            CHECK(is[4] == a3);

            CHECK(is.dim(0) == 8);
            CHECK(is.stride(0) == 1);
            CHECK(is.dim(1) == 9);
            CHECK(is.stride(1) == 8);
            CHECK(is.dim(2) == 1);
            CHECK(is.stride(2) == 72);
            CHECK(is.dim(3) == 1);
            CHECK(is.stride(3) == 72);
            CHECK(is.dim(4) == 1);
            CHECK(is.stride(4) == 72);
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
