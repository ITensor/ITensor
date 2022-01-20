#include "test.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/mps/sites/spinone.h"
#include "itensor/mps/sites/spintwo.h"
#include "itensor/mps/sites/electron.h"
#include "itensor/mps/sites/fermion.h"
#include "itensor/mps/sites/tj.h"
#include "itensor/util/print_macro.h"
#include "itensor/util/str.h"

using namespace itensor;

TEST_CASE("SiteSetTest")
{

    const int N = 10;
    std::string file_name("SiteSetIO.dat");


SECTION("Generic SiteSet")
    {
    auto sites = SiteSet(N,3);
    for(auto i : range1(N))
        {
        CHECK(dim(sites(i)) == 3);
        CHECK(hasTags(sites(i),"Site,n="+str(i)));
        }
    }

SECTION("SpinHalf (QNs)")
    {
    auto sites = SpinHalf(N,{"ConserveQNs=",true});
    for(auto i : range1(N))
        {
        CHECK(dim(sites(i)) == 2);
        CHECK(hasTags(sites(i),"Site,S=1/2,n="+str(i)));
        }

    op(sites,"Sz",2); 
    op(sites,"S+",2); 
    op(sites,"S-",2); 
    op(sites,"Sp",2); 
    op(sites,"Sm",2); 
    //TODO: test these throw the correct error
    //op(sites,"Sx",2); 
    //op(sites,"Sy",2); 
    //op(sites,"ISy",2); 
    }

SECTION("hasQNs")
    {
    auto sites = SpinHalf(N,{"ConserveQNs=",true});
    CHECK(hasQNs(sites));
    auto sitesNoQNs = SpinHalf(N,{"ConserveQNs=",false});
    CHECK(not hasQNs(sitesNoQNs));
    }

SECTION("SpinHalf (no QNs)")
    {
    auto sites = SpinHalf(N,{"ConserveQNs=",false});
    for(auto i : range1(N))
        {
        CHECK(dim(sites(i)) == 2);
        CHECK(hasTags(sites(i),"Site,S=1/2,n="+str(i)));
        }

    op(sites,"Sz",2); 
    op(sites,"S+",2); 
    op(sites,"S-",2); 
    op(sites,"Sp",2); 
    op(sites,"Sm",2); 
    op(sites,"Sx",2); 
    op(sites,"Sy",2); 
    op(sites,"ISy",2); 
    }

SECTION("SpinOne")
    {
    auto sites = SpinOne(N,{"ConserveQNs=",false});
    for(auto i : range1(N)) 
        {
        CHECK(dim(sites(i)) == 3);
        CHECK(hasTags(sites(i),"Site,S=1,n="+str(i)));
        }

    op(sites,"Sz",2); 
    op(sites,"S+",2); 
    op(sites,"S-",2); 
    op(sites,"Sp",2); 
    op(sites,"Sm",2); 
    op(sites,"Sx",2); 
    op(sites,"Sy",2); 
    op(sites,"ISy",2); 
    }

SECTION("SpinTwo")
    {
    auto sites = SpinTwo(N,{"ConserveQNs=",false});
    for(auto i : range1(N))
        {
        CHECK(dim(sites(i)) == 5);
        CHECK(hasTags(sites(i),"Site,S=2,n="+str(i)));
        }

    op(sites,"Sz",2);
    op(sites,"S+",2);
    op(sites,"S-",2);
    op(sites,"Sp",2);
    op(sites,"Sm",2);
    op(sites,"Sx",2);
    op(sites,"Sy",2);
    op(sites,"ISy",2);
    }

SECTION("Electron")
    {
    auto sites = Electron(N,{"ConserveQNs=",true});
    for(auto i : range1(N))
        {
        CHECK(dim(sites(i)) == 4);
        CHECK(hasTags(sites(i),"Site"));
        }

    op(sites,"Nup",2); 
    op(sites,"Ndn",2); 
    op(sites,"Nupdn",2); 
    op(sites,"Ntot",2); 
    op(sites,"Sz",2); 
    op(sites,"Cup",2); 
    op(sites,"Cdn",2); 
    op(sites,"Aup",2); 
    op(sites,"Adn",2); 
    op(sites,"F",2); 
    }

SECTION("Fermion")
    {
    auto sites = Fermion(N,{"ConserveQNs=",true});
    for(auto i : range1(N))
        {
        CHECK(dim(sites(i)) == 2);
        CHECK(hasTags(sites(i),"Site,Fermion,n="+str(i)));
        }

    op(sites,"N",2); 
    op(sites,"A",2); 
    op(sites,"Adag",2); 
    op(sites,"F",2); 
    }

SECTION("tJ")
    {
    auto sites = tJ(N,{"ConserveQNs=",true});
    for(auto i : range1(N))
        {
        CHECK(dim(sites(i)) == 3);
        CHECK(hasTags(sites(i),"Site,tJ,n="+str(i)));
        }

    op(sites,"Nup",2); 
    op(sites,"Ndn",2); 
    op(sites,"Aup",2); 
    op(sites,"Adn",2); 
    op(sites,"F",2); 
    }


    SECTION("Spin 1/2 legacy IO")
    {
        // Dump sites to disk
        {
            SiteSet sites=SpinHalf(N,{"ConserveQNs=",false});
            std::ofstream file(file_name);
            sites.write(file);
        }
        {
            //Restore from Disk (sort of)
            SiteSet sites;
            std::ifstream file(file_name);
            sites.read(file); //Creates a set of GenericSite objects.
            CHECK(sites.length()==N);//This is the easy part.

            // If we try and use the sites for something it fails with:
            //       'op' method not defined for generic site
    //        sites.op("Sz",1); //auto MPO will fail for the same reason.            
        }
        {
            //Restore from Disk with hind sight about the SiteSet type.
            SpinHalf sites;
            std::ifstream file(file_name);
            sites.read(file); //Creates a set of GenericSite objects.
            CHECK(sites.length()==N);//This is the easy part.
            sites.op("Sz",1); //auto MPO will fail for the same reason.            
        }
        {
            //Restore from Disk with *incorrect* hind sight about the SiteSet type.
            //SpinOne sites; this will actually work!! because of spin 1/2 the edge site handling code.
            Electron sites;
            std::ifstream file(file_name);
            sites.read(file); //Creates a set of GenericSite objects.
            CHECK(sites.length()==N);//This is the easy part.
            //auto sz=sites.op("Sz",1); Out of range: IndexVal at position 1 has val 3 > (dim=2|id=646|"n=1,Site,S=1/2")      
            //auto sx=sites.op("Sx",1); Not defined for Electron .. but it could be        
            //auto nup=sites.op("Nup",1); Out of range: IndexVal at position 1 has val 4 > (dim=2|id=413|"n=1,Site,S=1/2")          
        }

        std::remove(file_name.c_str()); //cleanup (CATCH2 should provide some sort of teardown function for this?)
    }

    SECTION("Spin 1 no QNs legacy IO")
    {
        // Dump sites to disk
        {
            SiteSet sites=SpinOne(N,{"ConserveQNs=",false,"SHalfEdge",true});
            std::ofstream file(file_name);
            sites.write(file);
        }
        //Restore from Disk (sort of)
        SiteSet sites;
        std::ifstream file(file_name);
        sites.read(file); //Creates a set of GenericSite objects.
        CHECK(sites.length()==N);
        CHECK(hasTags(sites(1  ),"S=1/2")); //Make sure we the edge spins correct after readback.
        CHECK(hasTags(sites(2  ),"S=1"  ));
        CHECK(hasTags(sites(N-1),"S=1"  ));
        CHECK(hasTags(sites(N  ),"S=1/2"));
        std::remove(file_name.c_str()); //cleanup (CATCH2 should provide some sort of teardown function for this?)
    }
    
    SECTION("Spin 2 no QNs legacy IO")
    {
        // Dump sites to disk
        {
            SiteSet sites=SpinTwo(N,{"ConserveQNs=",false,"SHalfEdge",true});
            std::ofstream file(file_name);
            sites.write(file);
        }
        //Restore from Disk (sort of)
        SiteSet sites;
        std::ifstream file(file_name);
        sites.read(file); //Creates a set of GenericSite objects.
        CHECK(sites.length()==N);
        CHECK(hasTags(sites(1  ),"S=1/2")); //Make sure we the edge spins correct after readback.
        CHECK(hasTags(sites(2  ),"S=2"  ));
        CHECK(hasTags(sites(N-1),"S=2"  ));
        CHECK(hasTags(sites(N  ),"S=1/2"));
        std::remove(file_name.c_str()); //cleanup (CATCH2 should provide some sort of teardown function for this?)
    }

    SECTION("Spin 1 with QNs legacy IO")
    {
        // Dump sites to disk
        {
            SiteSet sites=SpinOne(N,{"ConserveQNs=",true,"SHalfEdge",true});
            std::ofstream file(file_name);
            sites.write(file);
        }
        //Restore from Disk (sort of)
        SiteSet sites;
        std::ifstream file(file_name);
        sites.read(file); //Creates a set of GenericSite objects.
        CHECK(sites.length()==N);
        
        //Make sure we the edge spins correct after readback.
        CHECK(hasTags(sites(1  ),"S=1/2")); 
        CHECK(hasTags(sites(2  ),"S=1"  ));
        CHECK(hasTags(sites(N-1),"S=1"  ));
        CHECK(hasTags(sites(N  ),"S=1/2"));
        
        //Make sure the edge spins correct QNs after readback.
        CHECK(hasQNs(sites));
        CHECK(sites(1).qn(1).store()[0].name()=="Sz");
        CHECK(sites(1).qn(1).store()[0].val ()==1);
        CHECK(sites(1).qn(2).store()[0].name()=="Sz");
        CHECK(sites(1).qn(2).store()[0].val ()==-1);
        CHECK(sites(2).qn(1).store()[0].name()=="Sz");
        CHECK(sites(2).qn(1).store()[0].val ()==2);
        CHECK(sites(2).qn(2).store()[0].name()=="Sz");
        CHECK(sites(2).qn(2).store()[0].val ()==0);
        CHECK(sites(2).qn(3).store()[0].name()=="Sz");
        CHECK(sites(2).qn(3).store()[0].val ()==-2);
        std::remove(file_name.c_str()); //cleanup (CATCH2 should provide some sort of teardown function for this?)
    }

}



