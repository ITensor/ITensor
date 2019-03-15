#include "test.h"
#include "itensor/indexset.h"
#include "itensor/util/set_scoped.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

TEST_CASE("IndexSetTest")
{
auto i1 = Index(1,"Link,i1");
auto i2 = Index(2,"Link,i2");
auto i3 = Index(3,"Link,i3");
auto i4 = Index(4,"Link,i4");
auto i5 = Index(5,"Link,i5");
auto i6 = Index(6,"Link,i6");
auto i7 = Index(7,"Link,i7");
auto i8 = Index(8,"Link,i8");
auto i9 = Index(9,"Link,i9");
auto i10 = Index(10,"Link,i10");
auto v1 = Index(2,"Vtype,v1");
auto w1 = Index(2,"Wtype,w1");

SECTION("Constructors")
    {
    SECTION("One")
        {
        auto is = IndexSet(i4);
        CHECK(order(is) == 1);
        CHECK(dim(is) == 4);
        CHECK(is[0] == i4);
        CHECK(is[0] == is.index(1));

        is = IndexSet(i1);
        CHECK(order(is) == 1);
        CHECK(dim(is) == 1);
        CHECK(is[0] == i1);
        CHECK(is[0] == is.index(1));
        }

    SECTION("Two")
        {
        auto is1 = IndexSet(i4,i3);
        CHECK(order(is1) == 2);
        CHECK(dim(is1) == 4*3);
        CHECK(is1[0] == i4);
        CHECK(is1[1] == i3);

        auto is2 = IndexSet(i3,i4);
        CHECK(is2.r() == 2);
        CHECK(dim(is2) == 4*3);
        CHECK(is2[0] == i3);
        CHECK(is2[1] == i4);
        }

    SECTION("THREE")
        {
        auto is = IndexSet(i3,i1,i4);
        CHECK(is.r() == 3);
        CHECK(dim(is) == 4*3);
        CHECK(is[0] == i3);
        CHECK(is[1] == i1);
        CHECK(is[2] == i4);
        }

    SECTION("TEN")
        {
        auto is = IndexSet(i3,i1,i4,i2,i5,i6,i7,i8,i10,i9);
        CHECK(is.r() == 10);
        CHECK(is[0] == i3);
        CHECK(is[1] == i1);
        CHECK(is[2] == i4);
        CHECK(is[3] == i2);
        CHECK(is[4] == i5);
        CHECK(is[5] == i6);
        CHECK(is[6] == i7);
        CHECK(is[7] == i8);
        CHECK(is[8] == i10);
        CHECK(is[9] == i9);
        }

    }

SECTION("PrimeLevelMethods")
{
SECTION("Prime All")
    {
    SECTION("Case 1")
        {
        IndexSet is(i2,i3,prime(i2),i4);
        is.prime();
        CHECK(is[0] == prime(i2));
        CHECK(is[1] == prime(i3));
        CHECK(is[2] == prime(i2,2));
        CHECK(is[3] == prime(i4));
        }
    SECTION("Case 2")
        {
        IndexSet is(i2,i3,prime(i2),i4);
        is.prime(4);
        CHECK(is[0] == prime(i2,4));
        CHECK(is[1] == prime(i3,4));
        CHECK(is[2] == prime(i2,5));
        CHECK(is[3] == prime(i4,4));
        }
    }

SECTION("Prime IndexType")
    {
    SECTION("Case 1")
        {
        IndexSet is(i2,w1,prime(i2),v1);
        is.prime("Link");
        CHECK(is[0] == prime(i2));
        CHECK(is[1] == w1);
        CHECK(is[2] == prime(i2,2));
        CHECK(is[3] == v1);
        }
    SECTION("Case 2")
        {
        IndexSet is(i2,w1,prime(i2),v1);
        is.prime("Wtype");
        CHECK(is[0] == i2);
        CHECK(is[1] == prime(w1));
        CHECK(is[2] == prime(i2));
        CHECK(is[3] == v1);
        }
    SECTION("Case 3")
        {
        IndexSet is(i2,w1,prime(i2),v1);
        is.prime("Wtype");
        is.prime("Vtype");
        CHECK(is[0] == i2);
        CHECK(is[1] == prime(w1));
        CHECK(is[2] == prime(i2));
        CHECK(is[3] == prime(v1));
        }
    SECTION("Case 4")
        {
        IndexSet is(i2,w1,prime(i2),v1);
        is.prime(4,"Wtype");
        is.prime(4,"Vtype");
        CHECK(is[0] == i2);
        CHECK(is[1] == prime(w1,4));
        CHECK(is[2] == prime(i2));
        CHECK(is[3] == prime(v1,4));
        }
    }

SECTION("Prime Indices")
    {
    SECTION("Case 1")
        {
        IndexSet is(i5,i2,prime(i2),i4);
        is.prime(i2,prime(i2));
        CHECK(is[0] == i5);
        CHECK(is[1] == prime(i2));
        CHECK(is[2] == prime(i2,2));
        CHECK(is[3] == i4);
        }
    SECTION("Case 2")
        {
        IndexSet is(i5,i2,prime(i2),i4);
        is.prime(prime(i2),i2);
        CHECK(is[0] == i5);
        CHECK(is[1] == prime(i2));
        CHECK(is[2] == prime(i2,2));
        CHECK(is[3] == i4);
        }

    SECTION("Case 3 - include inc")
        {
        IndexSet is(i5,i2,prime(i2),i4);
        is.prime(4,prime(i2),i2);
        CHECK(is[0] == i5);
        CHECK(is[1] == prime(i2,4));
        CHECK(is[2] == prime(i2,5));
        CHECK(is[3] == i4);
        }
    //TODO: add this check back in
    //SECTION("Check Error Condition")
    //    {
    //    IndexSet is(i5,i2,i3,i4);
    //    CHECK_THROWS_AS(is.prime(prime(i2)),ITError);
    //    }
    }

SECTION("Prime Mix of IndexType and Index")
    {
    SECTION("Case 1 Basic")
        {
        auto is = IndexSet(i5,i2,prime(i2,2),v1,prime(w1));
        is.prime(i2);
        is.prime("Vtype");
        CHECK(is[0] == i5);
        CHECK(is[1] == prime(i2));
        CHECK(is[2] == prime(i2,2));
        CHECK(is[3] == prime(v1));
        CHECK(is[4] == prime(w1));
        }
    SECTION("Case 2 Increment")
        {
        auto is = IndexSet(i5,i2,prime(i2,2),v1,prime(w1));
        is.prime(3,i2);
        is.prime(3,"Vtype");
        CHECK(is[0] == i5);
        CHECK(is[1] == prime(i2,3));
        CHECK(is[2] == prime(i2,2));
        CHECK(is[3] == prime(v1,3));
        CHECK(is[4] == prime(w1));
        }
    SECTION("Case 3 Other Order")
        {
        auto is = IndexSet(i5,i2,prime(i2,2),v1,prime(w1));
        is.prime(3,"Vtype");
        is.prime(3,i2);
        CHECK(is[0] == i5);
        CHECK(is[1] == prime(i2,3));
        CHECK(is[2] == prime(i2,2));
        CHECK(is[3] == prime(v1,3));
        CHECK(is[4] == prime(w1));
        }
    SECTION("Check Error: No Matching Index")
        {
        auto is = IndexSet(i5,i2,prime(i2,2),v1,prime(w1));
        //TODO: add debug check that index is in IndexSet
        //CHECK_THROWS_AS(prime(is,Vtype,i3),ITError);
        }
    SECTION("Check Error: Invalid Prime Levels")
        {
        auto is = IndexSet(i5,i2,prime(i2,2),v1,prime(w1));
        //TODO: add debug check that index is in IndexSet
        //CHECK_THROWS_AS(prime(is,Vtype,i2,2),ITError);
        }
    }

// For now, primeExcept has been removed
// TODO: consider putting it back?
//SECTION("Prime Except")
//    {
//    SECTION("Case 1")
//        {
//        IndexSet is(i5,i2,i3,i4);
//        primeExcept(is,i5);
//        CHECK(is[0] == i5);
//        CHECK(is[1] == prime(i2));
//        CHECK(is[2] == prime(i3));
//        CHECK(is[3] == prime(i4));
//        }
//    SECTION("Case 2")
//        {
//        IndexSet is(i5,i2,i3,i4);
//        primeExcept(is,i2,2);
//        CHECK(is[0] == prime(i5,2));
//        CHECK(is[1] == i2);
//        CHECK(is[2] == prime(i3,2));
//        CHECK(is[3] == prime(i4,2));
//        }
//    SECTION("Case 3")
//        {
//        IndexSet is(i5,i2,i3,i4);
//        primeExcept(is,i2,i4,2);
//        CHECK(is[0] == prime(i5,2));
//        CHECK(is[1] == i2);
//        CHECK(is[2] == prime(i3,2));
//        CHECK(is[3] == i4);
//        }
//    SECTION("Case 4")
//        {
//        IndexSet is(i5,i2,i3,i4);
//        primeExcept(is,i4,i2,2);
//        CHECK(is[0] == prime(i5,2));
//        CHECK(is[1] == i2);
//        CHECK(is[2] == prime(i3,2));
//        CHECK(is[3] == i4);
//        }
//    SECTION("Regression Test 1")
//        {
//        Index x("x",2,Xtype),
//              z("z",2,Ztype),
//              y("y",2,Ytype),
//              v("v",2,Vtype);
//        IndexSet is(prime(x,3),prime(y,1),prime(z,3),y,v);
//        primeExcept(is,Vtype,Ytype,-2);
//        CHECK(is[0]==prime(x,1));
//        CHECK(is[1]==prime(y,1));
//        CHECK(is[2]==prime(z,1));
//        CHECK(is[3]==      y   );
//        CHECK(is[4]==      v   );
//        }
//    }

SECTION("NoPrime Index")
    {
    SECTION("Case 1")
        {
        IndexSet is(i5,i2,prime(i2),prime(i4));
        is.noPrime(prime(i4));
        CHECK(is[0] == i5);
        CHECK(is[1] == i2);
        CHECK(is[2] == prime(i2));
        CHECK(is[3] == i4);
        }
    SECTION("Case 2")
        {
        IndexSet is(i5,i2,prime(i2),prime(i4,4));
        is.noPrime(prime(i4,4));
        CHECK(is[0] == i5);
        CHECK(is[1] == i2);
        CHECK(is[2] == prime(i2));
        CHECK(is[3] == i4);
        }
    }

SECTION("NoPrimeType")
    {
    SECTION("Case 0")
        {
        IndexSet is(i2,v1,w1,i4);
        is.prime(2);
        is.noPrime();
        CHECK(is[0] == i2);
        CHECK(is[1] == v1);
        CHECK(is[2] == w1);
        CHECK(is[3] == i4);
        }
    SECTION("Case 1")
        {
        IndexSet is(i2,v1,w1,i4);
        is.prime(2);
        is.noPrime("Vtype");
        CHECK(is[0] == prime(i2,2));
        CHECK(is[1] == v1);
        CHECK(is[2] == prime(w1,2));
        CHECK(is[3] == prime(i4,2));
        }
    SECTION("Case 2")
        {
        IndexSet is(i2,v1,w1,i4);
        is.prime(2);
        is.noPrime("Wtype");
        CHECK(is[0] == prime(i2,2));
        CHECK(is[1] == prime(v1,2));
        CHECK(is[2] == w1);
        CHECK(is[3] == prime(i4,2));
        }
    SECTION("Case 3")
        {
        IndexSet is(i2,v1,w1,i4);
        is.prime(2);
        // TODO: allow listing multiple TagSets?
        is.noPrime("Wtype");
        is.noPrime("Vtype");
        CHECK(is[0] == prime(i2,2));
        CHECK(is[1] == v1);
        CHECK(is[2] == w1);
        CHECK(is[3] == prime(i4,2));
        }
    SECTION("Case 4")
        {
        IndexSet is(i2,v1,w1,i4);
        is.prime(2);
        // TODO: allow listing multiple TagSets?
        is.noPrime("Vtype");
        is.noPrime("Wtype");
        CHECK(is[0] == prime(i2,2));
        CHECK(is[1] == v1);
        CHECK(is[2] == w1);
        CHECK(is[3] == prime(i4,2));
        }
    }

SECTION("Map Prime")
    {
    SECTION("Case 1")
        {
        auto is = IndexSet(i1,prime(i1));
        is.setPrime(2,i1);
        CHECK(is[0] == prime(i1,2));
        CHECK(is[1] == prime(i1,1));
        }
    SECTION("Case 2")
        {
        auto is = IndexSet(i1,prime(i1));
        is.setPrime(2,prime(i1));
        CHECK(is[0] == i1);
        CHECK(is[1] == prime(i1,2));
        }
    SECTION("Case 3")
        {
        auto is = IndexSet(i1,prime(i1),v1);
        is.replaceTags("0","2","Link");
        is.setPrime(4,v1);
        CHECK(is[0] == prime(i1,2));
        CHECK(is[1] == prime(i1,1));
        CHECK(is[2] == prime(v1,4));
        }
    }

} //PrimeLevelMethods

SECTION("IndexSet Iterator")
    {
    auto s1 = Index(2,"s1");
    auto s2 = Index(3,"s2"); 
    auto s3 = Index(4,"s3"); 
    auto s4 = Index(5,"s4"); 
    auto i = IndexSet(s1,s2,s3,s4); 
    auto c = i.begin(); 
    CHECK(*c == s1);
    CHECK(*(c+1) == s2);
    CHECK(*(c+2) == s3);
    CHECK(*(c+3) == s4);
    CHECK((c+4) == i.end());
    }

} //IndexSetTest
