#include "test.h"

#include "itensor/global.h"
#include "itensor/util/infarray.h"
#include "itensor/util/stats.h"

using namespace itensor;
using namespace std;

TEST_CASE("InfArray")
{

SECTION("Constructors")
    {
    SECTION("Default")
        {
        InfArray<Real,10> ia;
        CHECK(ia.arr_size()==10);
        CHECK(ia.size()==0);
        CHECK(ia.empty());
        }

    SECTION("Small")
        {
        InfArray<Real,10> ia(8);
        CHECK(ia.size()==8);
        CHECK(ia.vec_size()==0);
        CHECK(!ia.empty());
        }

    SECTION("Large")
        {
        InfArray<Real,10> ia(20);
        CHECK(ia.size()==20);
        CHECK(ia.vec_size()==20);
        CHECK(!ia.empty());
        }

    SECTION("Value")
        {
        InfArray<int,10> ia(8,2);
        CHECK(ia.size()==8);
        CHECK(ia.vec_size()==0);
        CHECK(!ia.empty());
        for(size_t j = 0; j < ia.size(); ++j)
            {
            CHECK(ia[j]==2);
            }
        }
    SECTION("Initializer List")
        {
        InfArray<int,10> ia = {1,2,3,4,5};
        CHECK(ia.size()==5);
        for(int j = 0; j < int(ia.size()); ++j)
            CHECK(ia[j]==j+1);

        InfArray<int,5> ib = {1,2,3,4,5,6,7};
        CHECK(ib.size()==7);
        for(int j = 0; j < int(ib.size()); ++j)
            CHECK(ib[j]==j+1);
        }
    }

SECTION("push_back")
    {
    InfArray<int,10> ib;
    for(int j = 0; j < int(ib.arr_size()); ++j)
        {
        ib.push_back(j);
        }
    CHECK(ib.vec_size()==0);

    InfArray<int,4> ia;
    int final_size = 80;
    CHECK(ia.size()==0);
    for(int j = 0; j < final_size; ++j)
        {
        ia.push_back(j);
        CHECK(ia.size()==j+1);
        }
    CHECK(ia.vec_size()==final_size);
    for(int j = 0; j < final_size; ++j)
        {
        CHECK(ia[j]==j);
        }
    }

SECTION("Iteration")
    {
    int size = 8;
    InfArray<int,10> ia(size);
    for(int j = 0; j < size; ++j)
        {
        ia[j] = j;
        }
    int count = 0;
    for(auto& el : ia)
        {
        CHECK(el == count);
        ++count;
        }

    size = 20;
    ia.resize(size);
    for(int j = 0; j < size; ++j)
        {
        ia[j] = j;
        }
    count = 0;
    for(auto& el : ia)
        {
        CHECK(el == count);
        ++count;
        }
    }

//SECTION("Fill")
//    {
//    int size = 8;
//    InfArray<int,10> ia(size);
//    for(int j = 0; j < size; ++j)
//        {
//        ia[j] = j;
//        }
//
//    ia.fill(1);
//    for(auto& el : ia)
//        {
//        CHECK(el == 1);
//        }
//
//    size = 20;
//    ia.resize(size);
//    for(int j = 0; j < size; ++j)
//        {
//        ia[j] = j;
//        }
//    ia.fill(7);
//    for(auto& el : ia)
//        {
//        CHECK(el == 7);
//        }
//    }

SECTION("Resize")
    {
    int size = 8;
    InfArray<int,10> ia(size);
    for(int j = 0; j < size; ++j)
        {
        ia[j] = j;
        }
    ia.resize(4);
    CHECK(ia.size()==4);
    int count = 0;
    for(auto& el : ia)
        {
        CHECK(el == count);
        ++count;
        }
    ia.resize(10);
    CHECK(ia.size()==10);

    size = 20;
    ia.resize(size);
    for(int j = 0; j < size; ++j)
        {
        ia[j] = j;
        }
    ia.resize(18);
    CHECK(ia.size()==18);
    count = 0;
    for(auto& el : ia)
        {
        CHECK(el == count);
        ++count;
        }
    ia.resize(30);
    CHECK(ia.size()==30);
    CHECK(ia.vec_size()==30);
    }

SECTION("Erase")
    {
    int size = 9;
    InfArray<int,10> ia(size);
    for(int j = 0; j < size; ++j)
        {
        ia[j] = j;
        }
    auto vals_to_erase = {5,8,1};
    for(auto val : vals_to_erase)
        {
        auto it = ia.begin();
        for(; it != ia.end(); ++it)
            if(*it == val) break;
        ia.erase(it);
        }
    CHECK(ia.size()==(size-vals_to_erase.size()));

    auto erased = [&vals_to_erase](int elt)
                  {
                  return stdx::count(vals_to_erase,elt) != 0;
                  };

    //print("New els:");
    int count = 0;
    for(auto& el : ia)
        {
        //print(" ",el);
        while(erased(count)) ++count;
        CHECK(el == count);
        ++count;
        }
    //println();
    }
}

TEST_CASE("Stats")
{
SECTION("Average")
    {
    auto st = Stats();
    auto v = std::vector<Real>{{1.,2.,3.,4.}};
    auto avg = 0.;
    for(auto& el : v) 
        {
        st.putin(el);
        avg += el;
        }
    avg /= v.size();

    CHECK_CLOSE(st.avg(), avg);
    }

SECTION("Err")
    {
    auto st = Stats();
    auto v = std::vector<Real>{{1.,2.,3.,4.}};
    auto avg = 0.;
    auto avg2 = 0.;
    for(auto& el : v) 
        {
        st.putin(el);
        avg += el;
        avg2 += el*el;
        }
    avg /= v.size();
    avg2 /= v.size();
    auto err = std::sqrt((avg2-avg*avg)/(v.size()-1));

    CHECK_CLOSE(st.err(),err);
    }

SECTION("BinErr")
    {
    auto st = Stats();
    int N = 100;
    for(int n = 1; n <= N; ++n) st.putin(Global::random());
    CHECK(st.binerr(10) < 0.51/std::sqrt(N-1));
    }
}

