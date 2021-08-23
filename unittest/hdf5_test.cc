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

SECTION("Dense ITensor")
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

    auto C = ITensor(i2,i3);
    n = 1;
    for(auto n2 : range1(i2.dim()))
    for(auto n3 : range1(i3.dim()))
        {
        C.set(i2=n2,i3=n3,n-n*Cplx_i);
        n += 1;
        }

    auto fo = h5_open("test.h5",'w');
    h5_write(fo,"itensor_R",R);
    h5_write(fo,"itensor_C",C);
    close(fo);

    auto fi = h5_open("test.h5",'r');
    auto read_R = h5_read<ITensor>(fi,"itensor_R");
    auto read_C = h5_read<ITensor>(fi,"itensor_C");

    CHECK(norm(R-read_R) < 1E-10);
    CHECK(norm(C-read_C) < 1E-10);
    }

SECTION("QN")
    {
    auto q = QN({"Sz",-1},{"Nf",1,-1});
    auto fo = h5_open("test.h5",'w');
    h5_write(fo,"qn_q",q);
    close(fo);

    auto fi = h5_open("test.h5",'r');
    auto read_q = h5_read<QN>(fi,"qn_q");

    CHECK(q == read_q);
    }

SECTION("QN Index")
    {
    auto N = 4;
    auto sites = SpinHalf(N);
    auto s1 = sites(1);
    auto fo = h5_open("test.h5",'w');
    h5_write(fo,"qnindex_s1",s1);
    close(fo);

    auto fi = h5_open("test.h5",'r');
    auto read_s1 = h5_read<Index>(fi,"qnindex_s1");
    CHECK(s1 == read_s1);
    }

SECTION("QN ITensor")
    {
    auto N = 4;
    auto s = SpinHalf(N);

    auto R = randomITensor(QN({"Sz",0}),prime(s(1)),prime(s(2)),dag(s(1)),dag(s(2)));

    auto I = randomITensor(QN({"Sz",0}),prime(s(1)),prime(s(2)),dag(s(1)),dag(s(2)));
    auto C = R + Cplx_i*I;

    auto fo = h5_open("test.h5",'w');
    h5_write(fo,"qnitensor_R",R);
    h5_write(fo,"qnitensor_C",C);
    close(fo);

    auto fi = h5_open("test.h5",'r');
    auto read_R = h5_read<ITensor>(fi,"qnitensor_R");
    auto read_C = h5_read<ITensor>(fi,"qnitensor_C");
    CHECK(norm(R - read_R) < 1E-10);
    CHECK(norm(C - read_C) < 1E-10);
    }

SECTION("MPS")
    {
    auto N = 4;
    auto s = SpinHalf(N,{"ConserveQNs=",false});

    auto M = randomMPS(s,4);

    auto ampo = AutoMPO(s);
    for(auto j = 1; j < N; ++j)
        {
        ampo += "Sx",j,"Sx",j+1;
        }
    for(auto j = 1; j <= N; ++j)
        {
        ampo += 1.0,"Sz",j;
        }
    auto H = toMPO(ampo);

    auto fo = h5_open("test.h5",'w');
    h5_write(fo,"mps_M",M);
    h5_write(fo,"mpo_H",H);
    close(fo);

    auto fi = h5_open("test.h5",'r');
    auto read_M = h5_read<MPS>(fi,"mps_M");
    auto read_H = h5_read<MPO>(fi,"mpo_H");

    CHECK(abs(inner(M,read_M)-1.0) < 1E-8);

    CHECK(abs(inner(M,H,M)-inner(M,read_H,M)) < 1E-8);
    }

}

