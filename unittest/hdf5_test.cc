#include "test.h"

#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace std;
using namespace itensor;


TEST_CASE("HDF5Test")
{

SECTION("TagSet")
    {
    auto to = TagSet("Red,Blue,n=3");
    auto fo = h5_open("test.h5",'w');
    h5_write(fo,"tagset",to);
    close(fo);

    auto fi = h5_open("test.h5",'r');
    auto ti = h5_read<TagSet>(fi,"tagset");
    //close(fi);

    CHECK(ti == to);
    }

SECTION("Index")
    {
    auto io = Index(3,"Link,n=1");
    auto fo = h5_open("test.h5",'w');
    h5_write(fo,"index",io);
    close(fo);

    auto fi = h5_open("test.h5",'r');
    auto ii = h5_read<Index>(fi,"index");
    //close(fi);

    CHECK(ii == io);
    }

SECTION("IndexSet")
    {
    auto i1 = Index(1,"Link,n=1");
    auto i2 = Index(2,"n=2,Blue");
    auto i3 = Index(3,"n=2,Red");
    auto iso = IndexSet(i2,i1,i3);
    auto fo = h5_open("test.h5",'w');
    h5_write(fo,"indexset",iso);
    close(fo);

    auto fi = h5_open("test.h5",'r');
    auto isi = h5_read<IndexSet>(fi,"indexset");
    //close(fi);

    for(auto n : range1(iso.length()))
        {
        CHECK(isi(n) == iso(n));
        }
    }

SECTION("ITensor")
    {
    auto i2 = Index(2,"n=2,Blue");
    auto i3 = Index(3,"n=2,Red");

    auto R = ITensor(i2,i3);
    auto n = 1;
    for(auto n2 : range1(i2.dim()))
    for(auto n3 : range1(i3.dim()))
        {
        R.set(i2=n2,i3=n3,n);
        n += 1;
        }

    auto fo = h5_open("test.h5",'w');
    h5_write(fo,"itensor_R",R);
    close(fo);

    auto fi = h5_open("test.h5",'r');
    auto read_R = h5_read<ITensor>(fi,"itensor_R");

    CHECK(norm(R-read_R) < 1E-10);
    }

}

