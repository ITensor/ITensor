#include "test.h"

#include "itensor/all.h"

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
    close(fi);

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
    close(fi);

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
    close(fi);

    for(auto n : range1(iso.length()))
        {
        CHECK(isi(n) == iso(n));
        }
    }

}

