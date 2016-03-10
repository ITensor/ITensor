#include "test.h"
#include "itensor/decomp.h"

using namespace itensor;


TEST_CASE("Spectrum Test")
    {
    SECTION("ITensor")
        {   
        Index l("l",4),
              s1("s1",2,Site),
              s2("s2",2,Site),
              r("r",4);

        ITensor t(l,s1,s2,r);
        t.randomize();
        t *= 1./t.norm();

        ITensor a(l,s1),d,b;
        Spectrum itspec = svd(t,a,d,b);

        Spectrum spec(d);
        for(int n = 1; n <= spec.numEigsKept(); ++n)
            {
            CHECK(fabs(spec.eig(n)-itspec.eig(n)) < 1E-5);
            }

        CHECK(!spec.hasQNs());
        }

    SECTION("IQTensor")
        {   
        Index s1u("Site1 Up",1,Site);
        Index s1d("Site1 Dn",1,Site);
        Index s2u("Site2 Up",1,Site);
        Index s2d("Site2 Dn",1,Site);
        Index l1u("Link1 Up",4,Link);
        Index l10("Link1 Z0",2,Link);
        Index l1d("Link1 Dn",4,Link);
        Index l2uu("Link2 UU",2,Link);
        Index l20("Link2 Z0",2,Link);
        Index l2dd("Link2 DD",2,Link);

        IQIndex S1("S1",
                   s1u,QN(+1),
                   s1d,QN(-1),Out);
        IQIndex S2("S2",
                   s2u,QN(+1),
                   s2d,QN(-1),Out);
        IQIndex L1("L1",
                   l1u,QN(+1),
                   l10,QN( 0),
                   l1d,QN(-1),
                   Out);
        IQIndex L2("L2",
                   l2uu,QN(+2),
                   l20,QN( 0),
                   l2dd,QN(-2),
                   Out);

        IQTensor phi(L1(5),S1(1),S2(2),L2(3));
        phi.randomize();
        phi *= 1./phi.norm();

        IQTensor A(L1,S1),D,B;
        Spectrum iqtspec = svd(phi,A,D,B);

        Spectrum spec(D);
        //PrintData(D);
        for(int n = 1; n <= spec.numEigsKept(); ++n)
            {
            CHECK(fabs(spec.eig(n)-iqtspec.eig(n)) < 1E-5);
            //printfln("%s %.10f",spec.qn(n),spec.eig(n));
            }
        CHECK(spec.hasQNs());
        }

    }
